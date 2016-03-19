import numpy as np
import re
import csv

import useful_funcs

class LogData (object):

    def __init__(self, label, energy, fluence, Flux, totYld):
        self.label=label # Simulation details (single string)
        self.energy=energy # Simulation energy (single integer)
        self.fluence=fluence # Array of fluence (n) steps in a simulation 
        self.Flux=Flux # Flux of each species  (m) at a given fluence step (n by m array)
        self.totYld=totYld # Total yield at each fluence step

def read_exptfile(fn):
    with open(fn,'r') as logfile:
        contents = logfile.read()
    logfile.close()

    exptName=[]
    energies = []
    totYld = []
    fluence=0
    exptFlux=0
    
    lines=contents.split('\n')
    for line in lines:
        col = line.split(',')
        if 'Energies' in line:
            Header = col
            for i in range(len(Header)-1):
                exptName.append(col[i+1])
            print len(exptName), exptName
            totYld= [[] for i in range(len(exptName))]
        else:
            print line
            energies.append(float(col[0]))
            for i in range(1,len(Header)):
#                print line[0], line[i]
                if col[i]=='':
                    totYld[i-1].append(None)
                else:
                    totYld[i-1].append(col[i])

    expt_data=LogData(exptName, energies, fluence, exptFlux, totYld)

    return expt_data
            
def read_logfile(fn):
    with open(fn,'r') as logfile:
        contents = logfile.read()
    logfile.close()

    regex = re.compile(r'\d')
    numbers = [int(s) for s in regex.findall(fn)]
    energy = ''.join(map(str,numbers[1:(len(numbers)-2)]))
    isbv = numbers[len(numbers)-1]
    descrip = "The energy={} and isbv={}".format(energy, isbv)

    lines = contents.split('\n')
    count = 1
    nE = 0

    e_set = 0
    readData = 0
    energline = 0
    potline = 0
    fluline = 0
    sysLine = 0
    targetline = 0
    sbeline=0
    elemName = []
    Amass = []
    DNS0 = []
    RHO = []
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
    isbv = 0
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
            if len(SBV) == nE*nE:
                if SBV[nE*nE-2]==0:
                    isbv = 1
                else: 
                    isbv = 3
        if 'incident   reflected   reemitted   sputtered   deposited' in line:
            readData = count
        if readData > 0 and count in range(readData+1, readData+1+nE):
            elementNum.append(int(columns[0]))
            incidentF.append(float(columns[1]))
            reflectedF.append(float(columns[2]))
            reemittedF.append(float(columns[3]))
            sputteredF.append(float(columns[4]))
            depositedF.append(float(columns[5]))

        count +=1

    totYld = np.sum(sputteredF)/fluenz
    simName.append('{}{}->{}'.format(Esim,elemName[0],target))
    simName.append('isbv={}'.format(isbv))
    simName.append('ipot={}'.format(ipot))
    simName.append('inel{}={}'.format(elemName[1],INEL0[1]))
    simName.append('Edispl{}={}'.format(elemName[1],Edispl[1]))
    simName.append('SBE{}-{}={}'.format(elemName[nE-2],elemName[nE-2],SBV[nE*nE-5]))
    simName.append('SBE{}-{}={}'.format(elemName[nE-2],elemName[nE-1],SBV[nE*nE-4]))
    simName.append('SBE{}-{}={}'.format(elemName[nE-1],elemName[nE-2],SBV[nE*nE-2]))
    simName.append('SBE{}-{}={}'.format(elemName[nE-1],elemName[nE-1],SBV[nE*nE-1]))
    simName.append('Fluence={}'.format(fluenz))    
    simName.append('Thickness={}'.format(ttdyn))
    simName.append('Hist={}'.format(nh))

    log_data = LogData(simName, Esim, fluenz, sputteredF, totYld)
    return log_data

def read_sputfile(fn):
    with open(fn,'r') as sputfile:
        contents = sputfile.read()
    sputfile.close()

    regex = re.compile(r'\d')
    numbers = [int(s) for s in regex.findall(fn)]
    energy = ''.join(map(str,numbers[4:(len(numbers)-2)]))
    isbv = numbers[len(numbers)-1]
    descrip = "The energy={} and isbv={}".format(energy, isbv)


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
        # First 7 lines should ALWAYS be the same. Read in values as 'line'.
        # Strip removes whitespace/carrage return from EOL
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
#                print specFlu_cur
                numcur+=1
            if numcur == ncp:
                for i in range(ncp):
                    specFlu[i].append(specFlu_cur[i])
                    totYld_cur+=specFlu_cur[i]
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

    sput_data = LogData(simName, float(energy), fsteps, specFlu, totYld)
    return sput_data

def average(sput):
    sput_data=[]
    for j in range(len(sput)):
        iyld=[]
        fs=sput[j].fluence[1::]
        for k in range(len(sput[j].Flux)):
            smflux = useful_funcs.savitzky_golay(np.array(sput[j].Flux[k][1::]), 41, 3)
            iyld.append(smflux)
        smooth=useful_funcs.savitzky_golay(np.array(sput[j].totYld[1::]), 41, 3)
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
#        print fs, totYld
        sput_data.append(LogData(sput[j].label, sput[j].energy, fs, iyld, totYld))

    return sput_data

