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
        self.Flux=Flux # Flux of each species (m) at a given fluence step (n by m array)
        self.totYld=totYld # Total yield at each fluence step
    def __repr__(seld):
        return repr((self.label, self.energy, self.fluence, self,Flux, self.totYld))

    
elemlist = ['H', 'He', 'C_g','Na','Mg','Al','Si','S','K','Ca','Ti','Mn','Fe','Ni','O']
ifrac = 0.01
relcorr = [0.00,0.00,0.01,42.13,8.80,3.85,1.000,0.01,229.488,7.00,5.71,8.81,1.78,2.886,0.01]
navg=100

def f4(seq): 
   # order preserving
   noDupes = []
   [noDupes.append(i) for i in seq if not noDupes.count(i)]
   return noDupes

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
    sN=[]
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
        elif 'filname' in line:
            ln=line.split()
            ln=ln[-1].split('_')[-4]+'->'+line.split('_')[-3]
            filename.append(ln)
        elif '->' in line:
            sN.append(line) # The full simulation description (energy + ion -> target)
            simeng.append(columns[0])
            incion.append(columns[1])
            target.append(columns[3])
        elif 'isbv model' in line:
            isbv.append(columns[5])
        elif 'outgasing' in line:
            diff = line.split(':')
            outgass = '{} diffusion ON'.format(diff[1].rstrip())
        elif 'INFO' in line and 'flc' in line:
            totflu=columns[3]
        elif 'Fluence step' in line:
            flu_step=columns[3] 
        elif r1.match(line) is not None:
            flulines.append(line)
            if float(columns[0])==1:
                flu_start.append(len(flulines)-1)
                iYsum=columns.index('Ysum:')
                iqumax=columns.index('qumax')
                ncp=iqumax-iYsum-1
                nqumax=len(columns)-iqumax-1
                nelem.append(ncp)
                #print 'The number of elements is {}'.format(nYsum)
                    
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
        if len(isbv)==0:
            isbv.append('isbv1')
        elif len(isbv)<len(flu_start):
            isbv.append(isbv[-1])

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
        descrip.append(sN[n])
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
            incident = columns[:2] # The incident ion [1] and energy [0]
            target = columns[3] # The simulation target
            sN=line # The full description (incident -> target)
            sN=sN.strip()
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

    #print DNS0
    #print Q0
    #tardensity=np.sum()
    totYld = np.sum(sputteredF[1:])
    simName.append(sN)
    simName.append('isbv={}'.format(1))
    simName.append('ipot={}'.format(ipot))
    simName.append(target)
    #print elemName
    simName.append(elemName)
    simName.append(Amass)
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
            smflux = useful_funcs.savitzky_golay(np.array(sput[j].Flux[k][1:]), span, 5)
            iyld.append(smflux)
        smooth=useful_funcs.savitzky_golay(np.array(sput[j].totYld[1:]), span, 5)
        totYld=smooth.tolist()
        sput_data.append(LogData(sput[j].label,sput[j].energy,fs,iyld,totYld))

    return sput_data

def find_sputvar(sput, i):
    sput_data=[]
    for j in range(len(sput)):
        iyld=[]
        fs=sput[j].fluence[i]

        for k in range(len(sput[j].Flux)):
            iyld.append(sput[j].Flux[k][i])
        totYld=sput[j].totYld[i]

        sput_data.append(LogData(sput[j].label, sput[j].energy, fs, iyld, totYld))

    return sput_data

def yld_comp(out_yld, tar):
    # Array of neut and ion yields vs. fluence for all elements
    nfyld=[[0]*len(elemlist) for i in range(len(tar))] 
    ifyld=[[0]*len(elemlist) for i in range(len(tar))]
    # Array of average neut and ion yields at steady state for all elements
    avg_nyld=[[0]*len(elemlist) for i in range(len(tar))]
    avg_iyld=[[0]*len(elemlist) for i in range(len(tar))]
    # Arrays for returned LogData objects
    tot_nyld=[] 
    tot_iyld=[]
    dd = [] 
      
    k=-1
    # out_yld = files
    # cycle over ion/target combinations
    for i in range(len(out_yld)):
        print out_yld[i].label[0]
        ion = out_yld[i].label[0].split('->')[0]
        target = out_yld[i].label[0].split('->')[1]
        tyld=np.mean(out_yld[i].totYld[len(out_yld[i].totYld)-navg:])
        dd.append(len(out_yld[i].totYld)) # number of fluence steps

        k+=1
        for y, elem in enumerate(elemlist):
            match = [l for l, x in enumerate(out_yld[i].label[3]) if x == elem]
            if match:
                yld=out_yld[i].Flux[match[0]]
                #print k, y
                nfyld[k][y]=yld
                ifyld[k][y]=yld*relcorr[y]*ifrac

        for y, elem in enumerate(elemlist):
            if isinstance(nfyld[k][y], (int,long)):
                pass
            else:
                avg=np.mean(nfyld[k][y][dd[k]-navg:])
                if avg>1: avg-=1
                avg_nyld[k][y]=avg
                avg_iyld[k][y]=avg_nyld[k][y]*relcorr[y]*ifrac
                print '{}: Y^tot={:.5f}, Y^i={:.6f}'.format(elem,avg_nyld[k][y],avg_iyld[k][y])
        ttot=np.sum(avg_nyld[k])
        itot=np.sum(avg_iyld[k])
        print '{}: Y^tot={:.4f}, Y^i={:.6f}'.format(tar[k],ttot,itot)

        tyld=np.zeros(len(nfyld[k][-1]))
        iyld=np.zeros(len(ifyld[k][-1]))
        # Determine the total SW yield as a function of fluence
        # len(swtyld) = len(fluence)
        for y, elem in enumerate(elemlist):
            if not isinstance(nfyld[k][y],(int,long)) and y>1:
                tyld=map(add,nfyld[k][y],tyld) 
                iyld=map(add,ifyld[k][y],iyld)

        lbl=out_yld[i].label[:]
        lbl[3]=elemlist
        tot_nyld.append(LogData(lbl,out_yld[i].energy,out_yld[i].fluence,nfyld[k],tyld))
        tot_iyld.append(LogData(lbl,out_yld[i].energy,out_yld[i].fluence,ifyld[k],iyld))
        
    return tot_nyld, tot_iyld
        
def comp_yield(out_yld, tar):
    hyld=[[0]*len(elemlist) for i in range(len(tar))]
    h_iyld=[[0]*len(elemlist) for i in range(len(tar))]
    heyld=[[0]*len(elemlist) for i in range(len(tar))]
    he_iyld=[[0]*len(elemlist) for i in range(len(tar))]
    sw_avg=[[0]*len(elemlist) for i in range(len(tar))]
    sw_iavg=[[0]*len(elemlist) for i in range(len(tar))]
    sw_yld=[]
    sw_iyld=[]
    totyld=[]
    totiyld=[]
    swtot_yld = []
    swion_yld = []
    dd = []
    Elem_Ratio=[[0]*4 for i in range(len(tar))]
    #elemlist = ['H', 'He', 'C_g','Na','Mg','Al','Si','S','K','Ca','Ti','Mn','Fe','Ni','O']
    nE=[4, 5, 9, 12]

    for i in range(len(out_yld)): # cycle over files
        for j in range(len(out_yld[i])): # cycle over ion/target combinations
            ion = out_yld[i][j].label[0].split()[1]
            target = out_yld[i][j].label[0].split()[3]
            #print ion, target
            tyld=np.mean(out_yld[i][j].totYld[len(out_yld[i][j].totYld)-navg:])
            dd.append(len(out_yld[i][j].totYld))
            
            for k in range(len(tar)):
                if ion == 'H' and target in tar[k]:
                    #print 'The {} total yield is {}'.format(out_yld[i][j].label[0], tyld)
                    for y, elem in enumerate(elemlist):
                        match = [l for l, x in enumerate(out_yld[i][j].label[3]) if x == elem]
                        if match:
                            yld=0.95*out_yld[i][j].Flux[match[0]]
                            hyld[k][y]=yld
                            h_iyld[k][y]=yld*relcorr[y]*ifrac
    
                elif ion == 'He' and target in tar[k]:
                    #print 'The {} total yield is {}'.format(out_yld[i][j].label[0], tyld)
                    for y, elem in enumerate(elemlist):
                        match = [l for l, x in enumerate(out_yld[i][j].label[3]) if x == elem]
                        if match:
                            yld=0.05*out_yld[i][j].Flux[match[0]]
                            heyld[k][y]=yld
                            he_iyld[k][y]=yld*relcorr[y]*ifrac

    for k in range(len(tar)):
        print tar[k]
        for y, elem in enumerate(elemlist):
            if isinstance(hyld[k][y], (int,long)) or isinstance(heyld[k][y], (int,long)):
                if y<2:
                    havg=np.mean(hyld[k][0][dd[k]-navg:])
                    heavg=np.mean(heyld[k][1][dd[k]-navg:])
                else:
                    havg=0
                    heavg=0
                #print 'Y{}={} and Y{}={}'.format(elem, havg, elem, heavg)
            else:
                #print elemlist[match[0]], elem
                #print match[0]
                #print hyld[k][y]
                havg=np.mean(hyld[k][y][dd[k]-navg:])
                heavg=np.mean(heyld[k][y][dd[k]-navg:])
                sw_avg[k][y]=havg+heavg
                sw_iavg[k][y]=sw_avg[k][y]*relcorr[y]*ifrac
                print 'Y^tot({})={:.5f}, Y^i({})={:.6f}'.format(elem,sw_avg[k][y],elem,sw_iavg[k][y])
        ttot=np.sum(sw_avg[k][2:])
        itot=np.sum(sw_iavg[k][2:])
        print 'Y^tot({})={:.4f}, Y^i_({})={:.6f}'.format(tar[k],ttot ,tar[k],itot)
        
        for l in range(len(nE)):
            Elem_Ratio[k][l]=sw_iavg[k][nE[l]]/sw_iavg[k][6]

        sw_yld.append(map(add,hyld[k],heyld[k]))
        sw_iyld.append(map(add,h_iyld[k],he_iyld[k]))
        swtyld=np.zeros(len(sw_iyld[k][0]))
        swiyld=np.zeros(len(sw_iyld[k][0]))
        # Determine the total SW yield as a function of fluence
        # len(swtyld) = len(fluence)
        for y, elem in enumerate(elemlist):
            if not isinstance(sw_iyld[k][y], (int,long)) and y>1:
                swtyld=map(add,sw_yld[k][y],swtyld) 
                swiyld=map(add,sw_iyld[k][y],swiyld)
        totyld.append(swtyld)
        totiyld.append(swiyld)
        
    for k, ttar in enumerate(tar):
        tot_iyld=0
        tot_yld=0
        sw_iavg_c=0
        for y,elem in enumerate(elemlist):
            if not isinstance(sw_iyld[k][y], (int,long)) and y>1:
                havg=np.mean(hyld[k][y][dd[k]-navg:])
                heavg=np.mean(heyld[k][y][dd[k]-navg:])
                sw_iavg_c+=(havg+heavg)*relcorr[y]*ifrac

                elem_yld= np.mean(sw_yld[k][y][dd[k]-navg:])
                elem_iyld = np.mean(sw_iyld[k][y][dd[k]-navg:])
                
                tot_iyld+=elem_iyld
                tot_yld+=elem_yld
        #print ttar
        #print tot_yld, tot_iyld, sw_iavg_c
        sw_yld_lbl=out_yld[i][k].label[:]
        sw_yld_lbl[4]=elemlist
        asdf=sw_yld_lbl[1].replace('He','SW')
        sw_yld_lbl[1]=asdf
        asdf=sw_yld_lbl[1].replace('H','SW')
        sw_yld_lbl[1]=asdf
        swtot_yld.append(LogData(sw_yld_lbl,out_yld[i][k].energy,out_yld[i][k].fluence,sw_yld[k],totyld[k]))
        swion_yld.append(LogData(sw_yld_lbl,out_yld[i][k].energy,out_yld[i][k].fluence,sw_iyld[k],totiyld[k]))
        
    return Elem_Ratio, swtot_yld, swion_yld

def log_comp(yld, tar):
    hyld=[[] for i in range(len(tar))]
    h_iyld=[[] for i in range(len(tar))]
    heyld=[[] for i in range(len(tar))]
    he_iyld=[[] for i in range(len(tar))]
    sw_yld=[]
    sw_iyld=[]
    Elem_Ratio=[[0]*4 for i in range(len(tar))]
    nE=[4, 5, 9, 12]
    #print out_yld[i].label[4]
    #print out_yld[i].Flux
    #print tar
    neut_data=[]
    ion_data=[]

    for i in range(len(yld)):
        for k in range(len(tar)):
            ion = yld[i].label[0].split('->')[0]
            target = yld[i].label[0].split('->')[1]
            tyld=yld[i].totYld
            if ion == 'H' and target == tar[k]:
                #print 'The {} yield is {}'.format(yld[i].label[0], tyld)
                #print yld[i].label[4]
                for y, elem in enumerate(elemlist):
                    match = [l for l, x in enumerate(yld[i].label[4]) if x == elem]
                    if match:
                        ny=yld[i].Flux[match[0]]
                        hyld[k].append(ny)
                        #print elem, yld[i].label[4][match[0]], relcorr[y], ny
                        h_iyld[k].append(ny*relcorr[y]*ifrac)

            elif ion == 'He' and target == tar[k]:
                #print 'The {} yield is {}'.format(yld[i].label[0], tyld)
                for y, elem in enumerate(elemlist):
                    match = [l for l, x in enumerate(yld[i].label[4]) if x == elem]
                    if match:
                        cy=yld[i].Flux[match[0]]
                        heyld[k].append(cy)
                        he_iyld[k].append(cy*relcorr[y]*ifrac)

    for k in range(len(tar)):
        sw_yld.append(map(add,hyld[k],heyld[k]))
        sw_iyld.append(map(add,h_iyld[k],he_iyld[k]))
        neut_data.append(LogData(yld[k].label,yld[k].energy,yld[k].fluence,sw_yld[k],np.sum(sw_yld[k])))
        ion_data.append(LogData(yld[k].label,yld[k].energy,yld[k].fluence,sw_iyld[k],np.sum(sw_iyld[k])))
        npR=0
        if npR==1:
            print 'H Yield'
            for l in range(len(yld[0].label[4])):
                #print '{}, {:.4f}, {:.6f}'.format(out_yld[i][j].label[4][l],sw_yld[k][l],sw_iyld[k][l])
                print '{} {:.2e}'.format(yld[0].label[4][l], hyld[k][l])
            print 'He Yield'
            for l in range(len(yld[i].label[4])):
                print '{} {:.2e}'.format(yld[0].label[4][l], heyld[k][l])
        #print 'The {} total SW yield is {:.3f}'.format(tar[k], np.sum(sw_yld[k]))
        #print 'The {} total SW ion yield is {:.5f}'.format(tar[k], np.sum(sw_iyld[k]))

    return neut_data, ion_data
                        
