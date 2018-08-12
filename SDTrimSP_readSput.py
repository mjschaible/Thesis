import numpy as np
import re
import useful_funcs
import math

class LogData (object):
# Flux is misnamed and actually corresponds to elemYld
# The equation for the sputtered secondary ion flux at the surface of a small body of radius $r_{sb}$ 
# as a function of solar incidence angle is
#\begin{linenomath*}
#\begin{equation}\label{eq:surfFlux}
#	\Phi^+_i(r_{sb},\theta') = \Phi_{SW}\:Y^+_i P \cos^{-1.6}\theta'
#\end{equation}
#\end{linenomath*}
# where $P$ is a reduction factor due to regolith porosity and $\theta'$ is the incident ion angle.

    def __init__(self, label, energy, fluence, Flux, totYld):
        self.label=label # Simulation details (single string)
        self.energy=energy # Simulation energy (single integer)
        self.fluence=fluence # Array of fluence (n) steps in a simulation 
        self.Flux=Flux # Yield of each species (m) at a given fluence step (n by m array)
        self.totYld=totYld # Total yield at each fluence step
    def __repr__(seld):
        return repr((self.label, self.energy, self.fluence, self,Flux, self.totYld))

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
    totYld = np.sum(sputteredF[:])
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

def read_outfile(fn):
    counter = 0
    # Define arrays to contain the sim description [filename,sN]
    # as well as the fluence value at each simulation step, and the beginning/ending line
    # numbers for each specific meterorite in the file [flu]
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
            flulines.append(line) # Array containing all fluence values in sequence 

            # this array does not differentiate between meteorites and therefore
            # pointers to array positions are needed so that the array can be parsed
            # into the individual meteorites
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

    # Function to parse the individual meteorites from the [flulines] array containing all values
    # num_cur is the number of meteorites in the file. e.g. for the Aubrites class there are 
    # eight individual meteorite compositions -> num_cur=8
    # *Note that at this point each [fluline] contains an unbroken line of sputtering data for 
    #  all elements in the meteroite composition. These must be parsed into induvidual elements
    num_met=len(flu_end)
    data_vals=['']*num_met
    for n in range(num_met):
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
        
        # parse the sputtering line by line to form an array of sputtering yields at each fluence step
        # [ylds] contains an nElem x nFluence array of sputtering yields
        # while [energy], [dz], and [dz/dt] contain energy, erosion depth and rate 1 x nFluence arrays
        # [tot_yld] is the summation of the elemental yields at each fluence setp
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
        
        # Define a 'description' containing the simulation (energy + ion -> target)
        descrip.append(sN[n])
        #descrip.append(outgass)
        #descrip.append(isbv[n])
        # Define floating point arrays of all data values 
        arr_step=np.array(flustep,dtype='float')
        arr_eng=np.array(energy,dtype='float')
        arr_yld=np.array(ylds,dtype='float')
        arr_qumax=np.array(flu_qumax,dtype='float')
        arr_tyld=np.array(tot_yld,dtype='float')
        # Append a log data object containing the above value to a single array
        # of nMet x nFluence objects
        sim_data.append(LogData(descrip,arr_eng,arr_step,arr_yld,arr_tyld))
    
    return sim_data

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

def SGsmooth(sput):
# Function to smooth the yield data. 
# This seems largely unneccessary provided sufficient incident particles are tested
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
