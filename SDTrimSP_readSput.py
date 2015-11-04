def read_exptfile(fn):
    import numpy as np

    energies = []
    exptFlux_HeSiO2 = []

    for lines in open(fn):
        line = lines.split(',')
        if 'Energies' in lines:
            Header = line
            exptName = Header[2]
#            print Header
        elif line[2] == '':
            continue
        else:
            energies.append(float(line[0]))
            exptFlux_HeSiO2.append(float(line[2]))
#            print energies

    return exptName, energies, exptFlux_HeSiO2
            
def read_logfile(fn):
    import numpy as np
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
    targetline = 0
    sbeline=0
    elemName = []
    Amass = []
    DNS0 = []
    RHO = []
    Einit = []
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
            Einit.append(float(columns[1]))
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
            fluenz = columns[0] # x10^16 ions/cm^2/s
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
    
    simName.append('{}->{}'.format(elemName[0],target))
    simName.append('ipot={}'.format(ipot))
#    simName.append('isbv={}'.format(isbv))
    simName.append('inel{}={}'.format(elemName[1],INEL0[1]))
    simName.append('Edispl{}={}'.format(elemName[1],Edispl[1]))
    simName.append('SBE{}-{}={}'.format(elemName[nE-1],elemName[nE-1],SBV[nE*nE-1]))
    simName.append('Fluence={}'.format(fluenz))    
    simName.append('Thickness={}'.format(ttdyn))
    simName.append('Hist={}'.format(nh))
    
    return simName, fluenz, Einit, elemName, sputteredF

def read_sputfile(fn):
    import numpy as np
    with open(fn,'r') as sputfile:
        contents = sputfile.read()
    sputfile.close()

    lines = contents.split('\n')
    num=1
    pn = 0 
    fn = 0
    numSpec = 0
 
    # Initialize a list to hold the header steps
    fsteps = []
    sputdat1 = []
    sputdat2 = []
    sputdat3 = []

    for line in lines:
        # First 7 lines should ALWAYS be the same. Read in values as 'line'.
        # Strip removes whitespace/carrage return from EOL
 #       line = line.strip()
        if 'SDTrimSP' in line:
            version = line # Read in the version of SDTrimSP being used
        if '->' in line:
            descrip = line.split()
            target = line[4]
        if 'nh     idout       ncp   fluence' in line:
            paramLabel = line.split(None)
            pn = num
        if pn > 1 and num == pn+1:
            params = line.split(None)
            ncp = int(params[2]) # Determine # of particles in system
#            print "The number of species is", ncp
        if pn > 1 and num == pn+2:
            species = line.split(None)
        if 'flu step' in line:
            fluLabel = line
            fn = num 
        if fn > 0 and num == fn+1:
            fluence = line.split()
            fsteps.append(fluence[0])
        if 'by' in line:
            numSpec = numSpec+1
        if numSpec > 1:
            row = line.split(None)
            num_cols = len(row)
            if line.startswith(' by') or num_cols == 0:
                continue
            elif num_cols == 1:
                fsteps.append(row[0])
            else:
# I am cheating here. Should be able to account for an arbitrary number of species (not just sputdat1-3)
#                print row[0]
                sputdat1.append(row[0])
                sputdat2.append(row[1])
                sputdat3.append(row[2])
        num += 1

    numSteps = len(fsteps)
    specFlu1 = []
    specFlu2 = []
    specFlu3 = []

    # This part also needs updated to handle an arbitrary number of sample species
    for j in range(0,numSteps):
#        data_vals[j]=sputdat1[j]
        cur = j*3
        specTot1 = float(sputdat1[cur])+float(sputdat1[cur+1])+float(sputdat1[cur+2])
        specTot2 = float(sputdat2[cur])+float(sputdat2[cur+1])+float(sputdat2[cur+2])
        specTot3 = float(sputdat3[cur])+float(sputdat3[cur+1])+float(sputdat3[cur+2])  
        specFlu1.append(specTot1)
        specFlu2.append(specTot2)
        specFlu3.append(specTot3)
#        if j == numSteps-1:
#            print sputdat1[cur], sputdat1[cur+1], sputdat1[cur+2]
#            print specTot1, specTot2, specTot3

 #   print len(fsteps), len(specFlu1), len(specFlu2), len(specFlu3), numSteps
 #   print fsteps[numSteps-1], specFlu1[numSteps-1], specFlu2[numSteps-1], specFlu3[numSteps-1]

    return species, fsteps, specFlu1, specFlu2, specFlu3


def movingaverage(interval, window_size):
    import numpy as np
    aveFlux = []
#    window = np.zeros(int(window_size))/float(window_size)
    begtot=0
    endtot =0
    count=0
    for i in range(1,window_size):
        begtot += interval[i]
        begavg = begtot/i
        aveFlux.append(begavg)
 
#    print len(aveFlux)
    curRange = []
    for i in range(window_size,len(interval)-window_size):
        curRange = interval[i-window_size/2:i+window_size/2]
        aveFlux.append(sum(curRange)/window_size)
#    print len(aveFlux)    
#    aveFlux.append(np.convolve(interval, window, 'same'))
 
    for i in range(len(interval),len(interval)-window_size-1,-1):
        count +=1
        endtot += interval[i-1]
        endavg = endtot/count
        aveFlux.append(endavg)

    count=0
    endtot=0
    for i in range(len(interval),len(interval)-window_size,-1):
        count +=1
        endtot += interval[i-1]
    avYld = endtot/count

#    print avYld
#    print len(aveFlux)
#    print len(aveFlux), aveFlux
    return avYld, aveFlux

