def read_logfile(fn):
    import numpy as np
    with open(fn,'r') as logfile:
        contents = logfile.read()

    lines = contents.split('\n')
    count = 1
    readData = 0
    energline = 0
    sysLine = 0
    elemName = []
    elementNum=[]
    incidentF=[]
    reflectedF=[]
    reemittedF=[]
    sputteredF=[]
    depositedF=[]

    for line in lines:
        if 'SDTrimSP' in line:
            version = line # Read in the version of SDTrimSP being used
        if '->' in line:
            descrip = line
        if ' CPT          E0      AlPHA0       INEL0' in line:
            energline = count
        if energline > 0 and count == energline + 1:
            columns = line.split()
            energy = columns[1]
        if 'CPT  SYMBOL A-Z  A-MASS        DNS0         RHO' in line:
            sysLine = count
        if sysLine > 0 and count in range(sysLine+1, sysLine+4):
            columns = line.split()
            elemName.append(columns[1])
        if 'cpt    incident   reflected   reemitted   sputtered   deposited' in line:
            readData = count
#            print "The sputter data starts on line ", readData+1
        if readData > 0 and count in range(readData+1, readData+4):
            columns = line.split()
            elementNum.append(int(columns[0]))
            incidentF.append(float(columns[1]))
            reflectedF.append(float(columns[2]))
            reemittedF.append(float(columns[3]))
            sputteredF.append(float(columns[4]))
            depositedF.append(float(columns[5]))
#            print "Line ", count, ", Sputter Flux = ", columns[4]
        
        count +=1
    
#    print "The incident ion energy is " + name

    return energy, elemName, sputteredF

def read_sputfile(fn):
    import numpy as np
    with open(fn,'r') as sputfile:
        contents = sputfile.read()

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
            descrip = line
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

    sputfile.close()

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

    return descrip, species, fsteps, specFlu1, specFlu2, specFlu3


def movingaverage(interval, window_size):
    import numpy as np
    aveFlux = []
#    window = np.zeros(int(window_size))/float(window_size)
    begtot=0
    endtot =0
    count=0
    for i in range(1,window_size):
        begtot += interval[i-1]
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

#    print len(aveFlux)
#    print len(aveFlux), aveFlux
    return aveFlux
