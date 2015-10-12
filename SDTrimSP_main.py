import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from Tkinter import Tk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename

from SDTrimSP_readSput import read_sputfile
from SDTrimSP_readSput import movingaverage
from SDTrimSP_readSput import read_logfile
from SDTrimSP_readSput import read_exptfile
from SDTrimSP_plotSput import plot_sputFile
from SDTrimSP_plotSput import plot_sputYld
from SDTrimSP_plotSput import plot_sputExpt

class LogData (object):

    def __init__(self, simName, energ, Flux1, Flux2, Flux3):
        self.simName=simName
        self.energ=energ
        self.Flux1=Flux1
        self.Flux2=Flux2
        self.Flux3=Flux3

class SputEqData (object):

    def __init__(self, simNameeq, energeq, Flux1eq, Flux2eq, Flux3eq):
        self.simNameeq=simNameeq
        self.energeq=energeq
        self.Flux1eq=Flux1eq
        self.Flux2eq=Flux2eq
        self.Flux3eq=Flux3eq

class SputData (object):

    def __init__(self,species,fsteps,Flu1,Flu2,Flu3):
        self.species=species
        self.fsteps=fsteps
        self.Flu1=Flu1
        self.Flu2=Flu2
        self.Flu3=Flu3

# ------------ Begin main program ----------
cont = 1
num_runs=0
runs=[]
numFlu=2

# These variables take the class type
SSyld=[]
LogYld=[]
EqYld=[]


while cont != 0:
#    name = raw_input('Enter filename: ')
#    filename.append(name)

    simCount = 0

    if num_runs==0:
        print "Please identify the experimental data file."
        Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
        filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

        exptName, exptEng, exptFlux_HeSiO2 = read_exptfile(filename)
        b = plot_sputExpt(exptName, exptEng, exptFlux_HeSiO2, simCount)
        simCount +=1
    else:
        Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
        path = askdirectory() # show an "Open" dialog box and return the path to the selected file

        scanfiles = path+"/*.dat"                 
        print scanfiles
        # Define log file arrays
        simName=[]
        energies = []
        Flux1 = []
        Flux2 = []
        Flux3 = []
        
        #Define sputter data file arrays
        species=[]
        fsteps=[]
        SY1=[]
        SY2=[]
        SY3=[]

        energeq = []
        Flux1eq = []
        Flux2eq = []
        Flux3eq = []    

        for filename in glob.glob(scanfiles):
        #    if os.path.isfile(fn):
            print 'The file is {0}'.format(filename)
            
            # Read the log file
            if filename.find('E0') == -1:
                sN, fluence, energy, elemName, sputteredF = read_logfile(filename)
                name = energy + 'eV'

                simName.append(sN)
                energies.append(float(energy))
                Flux1.append(float(sputteredF[0])/float(fluence))
                Flux2.append(float(sputteredF[1])/float(fluence))
                Flux3.append(float(sputteredF[2])/float(fluence))
            
            #Read the sputter data file
            elif 'E0_33' in filename:
                regex = re.compile(r'\d')
                numbers = [int(s) for s in regex.findall(filename)]
                energy = ''.join(map(str,numbers[5:10]))
 
                species, fsteps, SY1, SY2, SY3 = read_sputfile(filename)

                specarr=np.array(species,dtype='string')
                fstepsarr=np.array(fsteps,dtype='float')
                specFlu1arr=np.array(SY1,dtype='float')
                specFlu2arr=np.array(SY2,dtype='float')
                specFlu3arr=np.array(SY3,dtype='float')

#                SSyld.append(SputData(specarr, fstepsarr, specFlu1arr, specFlu2arr, specFlu3arr))
                
                numave = 20 # specify number of points over which to average equilibrium yield
                avYld0, av0 = movingaverage(SY1,numave)
                avYld1, av1 = movingaverage(SY2,numave)
                avYld2, av2 = movingaverage(SY3,numave)       

                energeq.append(energy)
                Flux1eq.append(avYld0)
                Flux2eq.append(avYld1)
                Flux3eq.append(avYld2)

#                a = plot_sputFile(energy, species, fsteps, SY1, SY2, SY3, av0, av1, av2, numFlu)
            else:
                continue

        # After all files in a folder are read, sort the yields from the log and sputter data files
        # by energy
        simName_sorted = [x for (y,x) in sorted(zip(energies,Flux1), key=lambda pair:pair[0])]
        Flux1_sorted = [x for (y,x) in sorted(zip(energies,Flux1), key=lambda pair:pair[0])]
        Flux2_sorted = [x for (y,x) in sorted(zip(energies,Flux2), key=lambda pair:pair[0])]
        Flux3_sorted = [x for (y,x) in sorted(zip(energies,Flux3), key=lambda pair:pair[0])]
        energies_sorted=energies.sort()

        Flux1eq_sorted = [x for (y,x) in sorted(zip(energeq,Flux1eq), key=lambda pair:pair[0])]
        Flux2eq_sorted = [x for (y,x) in sorted(zip(energeq,Flux2eq), key=lambda pair:pair[0])]
        Flux3eq_sorted = [x for (y,x) in sorted(zip(energeq,Flux3eq), key=lambda pair:pair[0])]
        energeq_sorted=energeq.sort()
        
        # Convert log file data to arrays
        simNamearr =np.array(simName_sorted,dtype='string')
        energiesarr =np.array(energies_sorted,dtype='float')
        Flux1arr =np.array(Flux1_sorted,dtype='float')
        Flux2arr =np.array(Flux2_sorted,dtype='float')
        Flux3arr =np.array(Flux3_sorted,dtype='float')
        # Append log data to the class and plot
        LogYld.append(SputData(simNamearr, energiesarr, Flux1arr, Flux2arr, Flux3arr))

        plog_yld = plot_sputYld(simName, energies, Flux1_sorted, Flux2_sorted, Flux3_sorted, simCount)
        simCount +=1 

        # Convert equilibrium sputter file data to arrays
        energeqarr =np.array(energeq_sorted,dtype='float')
        Flux1eqarr =np.array(Flux1eq_sorted,dtype='float')
        Flux2eqarr =np.array(Flux2eq_sorted,dtype='float')
        Flux3eqarr =np.array(Flux3eq_sorted,dtype='float')
        # Append equilibrium sputter data to class and plot
        EqYld.append(SputEqData(simNamearr, energeqarr, Flux1eqarr, Flux2eqarr, Flux3eqarr))
        peq_yld = plot_sputYld(simName, energeq, Flux1eq_sorted, Flux2eq_sorted, Flux3eq_sorted, simCount)
        simCount +=1

        numFlu +=1

    #    print EqYld[num_runs]
    num_runs+=1
    cont = input('Enter 0 to end: ')

plt.show()





