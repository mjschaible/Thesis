import numpy as np
import re
import csv
import math

import useful_funcs
from operator import add

class LogData (object):

    def __init__(self, label, energy, fluence, Flux, totYld):
        self.label=label # Simulation details (single string)
        self.energy=energy # Simulation energy (single integer)
        self.fluence=fluence # Array of fluence (n) steps in a simulation 
        self.Flux=Flux # Flux of each species  (m) at a given fluence step (n by m array)
        self.totYld=totYld # Total yield at each fluence step
    def __repr__(seld):
        return repr((self.label, self.energy, self.fluence, self,Flux, self.totYld))

    
elemlist = ['C','Na','Mg','Al','Si','S','K','Ca','Ti','Fe','Ni','O']
ifrac = 0.1
relcorr = [0.01, 23.319, 2.886, 3.637, 1.000, 0.01, 81.080, 3.950, 3.582, 2.767, 2.886, 0.01]
navg=10

def read_exptfile(fn):
    with open(fn,'r') as logfile:
        contents = logfile.read()
    logfile.close()
    simname = fn.split('_')[-1].split('.srim')[0]
    descrip=[]
    exptName=[]
    ion=[]
    target=[]
    author=[]
    SBEO=[]
    SBEX=[]
    energies = []
    totYld = []
    fluence=0
    exptFlux=0
    
    lines=contents.split('\n')
    for line in lines:
        lin = line.strip('\r')
        col = lin.split(',')
        if 'Energies' in line:
            Header = col
            for i in range(len(Header)-1):
                ion_tar=col[i+1].split('->')
                ion.append(ion_tar[0])
                target.append(ion_tar[1])
                exptName.append(col[i+1])
            #print exptName
            totYld= [[] for i in range(len(exptName))]
        elif 'Authors' in line:
            for i in range(len(Header)-1):
                author.append(col[i+1])
        elif 'SBE_O' in line:
            for i in range(len(Header)-1):
                SBEO.append(col[i+1])
        elif 'SBE_X' in line:
            for i in range(len(Header)-1):
                SBEX.append(col[i+1])
        else:
            energies.append(float(col[0]))
            for i in range(1,len(Header)):
                if col[i]=='':
                    totYld[i-1].append(float('nan'))
                else:
                    totYld[i-1].append(float(col[i]))
    
    descrip.append(exptName)
    descrip.append(ion)
    descrip.append(target)
    descrip.append(author)
    descrip.append(SBEO)
    descrip.append(SBEX)
    
    expt_data=LogData(descrip, energies, fluence, exptFlux, totYld)
    return expt_data

def read_outfile(fn):
    counter = 0

    filename=[]
    simeng=[]
    incion=[]
    target=[]
    isbv=[]
    
    flulines=[]
    flu_start=[]
    flu_end=[]

    nelem=[]
    descrip=[]
    sumlines=[]
    histlines=[]

    sim_data=[]
    
    r1 = re.compile('.*:.*:.*:.*: .*')
    with open(fn,'r') as logfile:
        contents = logfile.read()
    logfile.close()

    lines=contents.split('\n')

    for line in lines:
        counter+=1
        columns = line.split()
        line = re.sub('([:])', r'\1 ', line)
        line = re.sub('(dt)', r'\1 ', line)
        
        if 'ERROR' in line:
            continue
        elif 'filename' in line:
            
            filename.append(line.split('_')[-1])
            #print filename[-1]
        elif '->' in line and 'eV' in line:
            #print line
            ncp=0
            simeng.append(columns[0])
            incion.append(columns[1])
            target.append(columns[3])
        elif 'isbv model' in line:
            isbv.append(columns[5])
        elif 'use inelastic' in line:
            ncp+=1
        elif 'outgasing' in line:
            diff = line.split(':')
            outgass = '{} diffusion ON'.format(diff[1].rstrip())
            print line
        elif 'INFO' in line and 'flc' in line:
            #print line
            totflu=columns[3]
        elif 'Fluence step' in line:
            flu_step=columns[3] 
        elif r1.match(line) is not None:
            flulines.append(line)
            if float(columns[0])==1:
                flu_start.append(len(flulines)-1)
                iYsum=columns.index('Ysum:')
                iqumax=columns.index('qumax')
                nYsum=iqumax-iYsum-1
                nqumax=len(columns)-iqumax-1
                if nYsum == nqumax and nYsum == ncp:
                    nelem.append(nYsum)
                    #print 'The number of elements is {}'.format(nYsum)
                else:
                    print "The number of elements does not match"
                    print 'nYsum = {}, nqumax = {}, ncp = {}'.format(nYsum, nqumax, ncp)
                    
            if float(columns[2])==float(totflu):
                flu_end.append(len(flulines))
        elif 'sum' in line:
            sumlines.append(line)
        elif 'ihist' in line:
            histlines.append(line)
        elif 'isbv' in line and 'move' not in line and 'WARNING' not in line:
            print line
        elif counter == len(lines):
            if len(flu_end) < len(flu_start):
                flu_end.append(counter-1)

        if len(isbv)<len(flu_start):
            isbv.append('isbv1')

    try:
        outgass
    except:
        outgass = 'Diffusion OFF'

    num_cur=len(flu_end)
    data_vals=['']*num_cur
    for n in range(num_cur):
        data_vals[n]=flulines[flu_start[n]:flu_end[n]]

        descrip=[]
        flustep=[]
        energy=[]
        dz=[]
        dzdt=[]
        flu_yld=[]
        flu_qumax=[]
        tot_yld=[]
        ylds=[[] for i in range(ncp)]
        tyld=0
        
        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()
            for i in range(10,10+ncp):
                ylds[i-10].append(float(vals[i]))
                if i>10:
                    tyld+=float(vals[i])
                                
            qumax=np.array(vals[10+ncp+1:len(vals)]).astype(np.float)
            flustep.append(vals[2])
            energy.append(float(vals[4]))
            dz.append(float(vals[6]))
            dzdt.append(float(vals[8]))

            flu_qumax.append(qumax)
            tot_yld.append(tyld)
            tyld=0

        descrip.append(filename[n])
        descrip.append('{}->{}'.format(incion[n],target[n]))
        descrip.append(outgass)
        descrip.append(isbv[n])
        arr_step=np.array(flustep,dtype='float')
        arr_eng=np.array(energy,dtype='float')
        arr_yld=np.array(ylds,dtype='float')
        arr_qumax=np.array(flu_qumax,dtype='float')
        arr_tyld=np.array(tot_yld,dtype='float')
        sim_data.append(LogData(descrip,arr_eng,arr_step,arr_yld,arr_tyld))
            
    return sim_data
    
def read_logfile(fn):
    with open(fn,'r') as logfile:
        contents = logfile.read()
    logfile.close()

    regex = re.compile(r'\d')
    stuff = fn.split('_')
    for i in range(len(stuff)):
        if 'eV' in stuff[i]:
            numbers = [int(s) for s in regex.findall(stuff[i])]
            energy = ''.join(map(str,numbers[:]))
    
    numbers = [int(s) for s in regex.findall(fn)]
    isbv = numbers[len(numbers)-1]
    descrip = "The energy={} and isbv={}".format(energy, isbv)
    #print descrip
    lines = contents.split('\n')
    count = 1
    nE = 0

    e_set = 0
    readData = 0
    energline = 0
    potline = 0
    fluline = 0
    sysLine = 0
    tarLine = 0
    targetline = 0
    sbeline=0
    elemName = []
    Amass = []
    DNS0 = []
    RHO = []
    Q0=[]
    Qbeam=[]
    Qmax=[]
    Charge=[]
    INEL0 = []
    Ecutoff=[]
    Edispl=[]
    Ebulkb=[]
    Esurfb=[]
    Delta_H_d = []
    fluenz = 0
    tt = 0
    ttdyn = 0
    dsf = 0
    ttemp = 0 # Sample temperature during irradiation
    nh = 0
    nr_proj = 0
    idrel = 0 # Mode of simulation (0 = fully dynamic)
    sfin = 0 # Inelastic energy loss outside surface (1=yes)
    ipot = 0 # Potential used during the simulation
    target = 't' # Multielement target name
    mdens = 0 # atomic density of target (g/cm^3)
    adens = 0 # atomix density of targer (atom/ang^3)
    deltahf = 0 # target formation enthalpy
    SBV = []
    elementNum=[]
    incidentF=[]
    reflectedF=[]
    reemittedF=[]
    sputteredF=[]
    totYld = []
    depositedF=[]

    simName = []

    for line in lines:
        columns = line.split()
        if 'SDTrimSP' in line:
            version = line # Read in the version of SDTrimSP being used
        if '->' in line and len(columns)>3:
            target = columns[3]
        if line.strip()!='' and sysLine > 0:
            elemName.append(columns[1])
            Amass.append(columns[3])
            DNS0.append(columns[4])
            RHO.append(columns[5])
            nE+=1
        else:
            sysLine=0
        if 'SYMBOL A-Z  A-MASS' in line:
            sysLine = count
        if 'Q-0      Q-BEAM' in line:
            tarLine=count
        if tarLine > 0 and count in range(tarLine+1, tarLine+1+nE):   
            columns=line.split()
            Q0.append(float(columns[1]))
            Qbeam.append(float(columns[2]))
            Qmax.append(float(columns[3]))
            Charge.append(float(columns[4]))
        if ' CPT          E0      AlPHA0       INEL0' in line:
            energline = count
        if energline > 0 and count in range(energline+1, energline+1+nE):
            if count == energline+1:
                Esim=float(columns[1]) # Incident ion energy
                Asim=float(columns[2]) # Incident ion angle
            INEL0.append(columns[3])        
        if 'CPT    E_CUTOFF     E_DISPL     E_BULKB     E_SURFB' in line:
            e_set = count
        if e_set > 0 and count in range(e_set+1, e_set+1+nE):   
            columns=line.split()
            if line.strip()=='':
                for i in range(3):
                    Ecutoff.append(0)
                    Edispl.append(0)
                    Ebulkb.append(0)
                    Esurfb.append(0)
                    Delta_H_d.append(0)
                e_set=0
            else:
                Ecutoff.append(float(columns[1]))
                Edispl.append(float(columns[2]))
                Ebulkb.append(float(columns[3]))
                Esurfb.append(float(columns[4]))
                if len(columns)>6:
                    Delta_H_d.append(float(columns[6]))
                else:
                    Delta_H_d.append('')        
        if 'FLUENZ          TT       TTDYN         DSF       TTEMP' in line:
            fluline = count
        if fluline > 0 and count == fluline + 1:
            fluenz = float(columns[0]) # x10^16 ions/cm^2/s
            tt = columns[1]
            ttdyn = columns[2]
            dsf = columns[3]
            ttemp = columns[4]            
        if 'NH  NR-PPROJ       NCP     IDREL      SFIN      IPOT' in line:
            potline = count
        if potline > 0 and count == potline + 1:
            nh = columns[0]
            nr_proj = columns[1]
            idrel = columns[3]
            sfin = columns[4]
            ipot = columns[5]
        if 'composit.   mol.mass    rhom  atoms/A**3   deltahf' in line:
            targetline = count
        if targetline > 0 and count == targetline + 1:
            target = columns[0]
            mdens = columns[2]
            adens = columns[3]
            deltahf = columns[4]
        if 'SBV(I,J)' in line:
            sbeline = count
        if sbeline > 0 and count in range(sbeline+2, sbeline+2+nE*nE):
            numCol = len(columns)
            SBV.append(float(columns[numCol-4]))
        if 'incident   reflected   reemitted   sputtered   deposited' in line:
            readData = count
        if readData > 0 and count in range(readData+1, readData+1+nE):
            elementNum.append(int(columns[0]))
            incidentF.append(float(columns[1]))
            reflectedF.append(float(columns[2]))
            reemittedF.append(float(columns[3]))
            sputteredF.append(float(columns[4])/fluenz)
            depositedF.append(float(columns[5]))

        count +=1
    print DNS0
    print Q0
    tardensity=np.sum()
    totYld = np.sum(sputteredF[1:])
    simName.append('{}->{}'.format(elemName[0],target))
    simName.append('isbv={}'.format(isbv))
    simName.append('ipot={}'.format(ipot))
    simName.append(elemName)
    simName.append(target)
    simName.append(Esurfb)
#    simName.append('inel{}={}'.format(elemName[1],INEL0[1]))
#    simName.append('Edispl{}={}'.format(elemName[1],Edispl[1]))
#    simName.append('SBE{}-{}={}'.format(elemName[nE-2],elemName[nE-2],SBV[nE*nE-5]))
#    simName.append('SBE{}-{}={}'.format(elemName[nE-2],elemName[nE-1],SBV[nE*nE-4]))
#    simName.append('SBE{}-{}={}'.format(elemName[nE-1],elemName[nE-2],SBV[nE*nE-2]))
#    simName.append('SBE{}-{}={}'.format(elemName[nE-1],elemName[nE-1],SBV[nE*nE-1]))
#    simName.append('Fluence={}'.format(fluenz))    
#    simName.append('Thickness={}'.format(ttdyn))
#    simName.append('Hist={}'.format(nh))

    log_data = LogData(simName, Esim, fluenz, sputteredF, totYld)
    return log_data

def read_sputfile(fn):
    with open(fn,'r') as sputfile:
        contents = sputfile.read()
    sputfile.close()

    regex = re.compile(r'\d')
    stuff = fn.split('_')
    for i in range(len(stuff)):
        if 'eV' in stuff[i]:
            numbers = [int(s) for s in regex.findall(stuff[i])]
            energy = ''.join(map(str,numbers[:]))
    
    numbers = [int(s) for s in regex.findall(fn)]
    isbv = numbers[len(numbers)-1]
    descrip = "The energy={} and isbv={}".format(energy, isbv)
    #print descrip
    
    lines = contents.split('\n')
    num=1
    pn = 0 
    fn = 0
    numSpec = 0
    numcur=0
    totYld_cur=0
 
    # Initialize a list to hold the header steps
    simName = []
    totYld=[]
    fsteps=[]
    for line in lines:
        line = line.strip()
        columns = line.split()
        if 'SDTrimSP' in line:
            version = line # Read in the version of SDTrimSP being used
        if num==2:
            short_title = line
        if 'nh     idout       ncp   fluence' in line:
            pn = num
        if pn > 1 and num == pn+1:
            nh = int(columns[0]) # number of simulated projectiles
            ncp = int(columns[2]) # number of species in system
            fluenz = float(columns[3]) # simulated total fluence (x10^16 ion/cm^2)
            specFlu= [[] for i in range(ncp)]
            specFlu_cur=[0]*ncp
        if pn > 1 and num == pn+2:
            species = columns
        if 'flu step' in line:
            fluLabel = line
            fn = num 
        if fn > 0 and num == fn+1:
            fsteps.append(float(columns[0]))
        if 'by' in line:
            numSpec+=1
        if numSpec > 0:
            num_cols = len(columns)
            if 'by' in line or num_cols == 0:
                continue
            elif num_cols == 1:
                fsteps.append(float(columns[0]))
                numcur=0
            elif numcur<ncp:
                for i in range(ncp):
                    specFlu_cur[i]+=float(columns[i])
                numcur+=1
            if numcur == ncp:
                for i in range(ncp):
                    specFlu[i].append(specFlu_cur[i])
                    #totYld_cur+=specFlu_cur[i]
                totYld_cur = np.sum(specFlu_cur[1:])
                totYld.append(totYld_cur)
#                print totYld, numcur
                specFlu_cur = [0]*ncp
                totYld_cur=0
        num += 1

    if len(specFlu[0]) != len(fsteps):
        print "There is a problem"

    simName.append('{}'.format(short_title))
    simName.append('isbv={}'.format(isbv))
    simName.append('Fluence={}'.format(fluenz))    
    simName.append('Hist={}'.format(nh))
    simName.append(species)

    sput_data = LogData(simName, float(energy), fsteps, specFlu, totYld)
    return sput_data

def average(sput):
    sput_data=[]
    for j in range(len(sput)):
        iyld=[]
        fs=sput[j].fluence[1:]
        span = math.ceil(len(fs)/2)
        if span % 2 ==0:
            span += 1
        for k in range(len(sput[j].Flux)):
            smflux = useful_funcs.savitzky_golay(np.array(sput[j].Flux[k][1:]), span, 3)
            iyld.append(smflux)
        smooth=useful_funcs.savitzky_golay(np.array(sput[j].totYld[1:]), span, 3)
        totYld=smooth.tolist()
        sput_data.append(LogData(sput[j].label,sput[j].energy,fs,iyld,totYld))

    return sput_data

def find_sputvar(sput, i):
    sput_data=[]
    for j in range(len(sput)):
        iyld=[]
        fs=sput[j].fluence[i]
        print sput[0].fluence
        print sput[0].Flux

        for k in range(len(sput[j].Flux)):
            iyld.append(sput[j].Flux[k][i])
        totYld=sput[j].totYld[i]
#        print fs, totYld

        sput_data.append(LogData(sput[j].label, sput[j].energy, fs, iyld, totYld))

    return sput_data

def comp_yield(out_yld, tar):
    hyld=[[0]*len(elemlist) for i in range(len(tar))]
    h_iyld=[[0]*len(elemlist) for i in range(len(tar))]
    heyld=[[0]*len(elemlist) for i in range(len(tar))]
    he_iyld=[[0]*len(elemlist) for i in range(len(tar))]
    sw_yld=[]
    sw_iyld=[]
    Elem_Ratio=[[0]*4 for i in range(len(tar))]
    nE=[2, 3, 7, 9]
    #print out_yld[i].label[4]
    #print out_yld[i].Flux
    print tar
    for i in range(len(out_yld)):
        for j in range(len(out_yld[i])):
            for k in range(len(tar)):
                ion = out_yld[i][j].label[1].split('->')[0]
                target = out_yld[i][j].label[1].split('->')[1]
                tyld=np.mean(out_yld[i][j].totYld[len(out_yld[i][j].totYld)-navg:])
                #print ion, target, tar[k]
                if ion == 'H' and target in tar[k]:
                    print 'The {} yield is {}'.format(out_yld[i][j].label[1], tyld)
                    for y, elem in enumerate(elemlist):
                        match = [l for l, x in enumerate(out_yld[i][j].label[4]) if x == elem]
                        if match:
                            yld=np.mean(0.95*out_yld[i][j].Flux[match[0]][len(out_yld[i][j].totYld)-navg:])
                            hyld[k][y]=yld
                            h_iyld[k][y]=yld*relcorr[y]*ifrac

                elif ion == 'He' and target == tar[k]:
                    print 'The {} yield is {}'.format(out_yld[i][j].label[1], tyld)
                    for y, elem in enumerate(elemlist):
                        match = [l for l, x in enumerate(out_yld[i][j].label[4]) if x == elem]
                        if match:
                            yld=np.mean(0.05*out_yld[i][j].Flux[match[0]][len(out_yld[i][j].totYld)-navg:])
                            heyld[k][y]=yld
                            he_iyld[k][y]=yld*relcorr[y]*ifrac

    for k in range(len(tar)):
        sw_yld.append(map(add,hyld[k],heyld[k]))
        sw_iyld.append(map(add,h_iyld[k],he_iyld[k]))
        for i in range(len(nE)):
            Elem_Ratio[k][i]=sw_iyld[k][nE[i]]/sw_iyld[k][4]
        print 'The {} total SW yield is {:.3f}'.format(tar[k], np.sum(sw_yld[k]))
        print 'The {} total SW ion yield is {:.5f}'.format(tar[k], np.sum(sw_iyld[k]))

    #print Elem_Ratio
    
    return Elem_Ratio
