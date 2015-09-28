import numpy as np
import matplotlib.pyplot as plt
import glob
from readSput import read_sputfile
from readSput import movingaverage
from readSput import read_logfile
from plotSput import plot_sputFile

## First open and scan the file to determine the total number of lines, where the header lines end, and where the footer lines begin                
# To read in from command line, use:
#filename = raw_input('Enter filename: ')                                       

# To input a specific file
#filename = 'E0_33_1000eVHe_SiO2.dat'
path = "./*.dat"

energies = []
Flux1 = []
Flux2 = []
Flux3 = []
simCount = 0

# open, read only ('r') file. Pull out sim version, description, and parameters
#with open(filename, 'r') as logfile:
for filename in glob.glob(path):
#    if os.path.isfile(fn):
#    print('The file is ' + filename)

    if 'E0_33' in filename:
        descrip, species, fsteps, specFlu1, specFlu2, specFlu3 = read_sputfile(filename)

        #print specFlu1
        numave = 20

        av0 = movingaverage(specFlu1,numave)
        av1 = movingaverage(specFlu2,numave)
        av2 = movingaverage(specFlu3,numave)

        a = plot_sputFile(descrip, species, fsteps, specFlu1, specFlu2, specFlu3, av0, av1, av2)

    if 'log.dat' in filename:
        energy, elemName, sputteredF = read_logfile(filename)

        name = energy + 'eV'
 #       print sputteredF
        energies.append(energy)
        Flux1.append(sputteredF[0])
        Flux2.append(sputteredF[1])
        Flux3.append(sputteredF[2])
 

        fig = plt.figure(2)
        ax1 = fig.add_subplot(111)
        #ax1.set_title("Sputtered Fluence")    
        ax1.set_xlabel('Energy')
        ax1.set_ylabel('Sputtered Fluence')

        ax1.scatter(energies,Flux2, c='r', label="{0}".format(elemName[1]))
        ax1.scatter(energies,Flux3, c='b', label="{0}".format(elemName[2]))    
        if simCount == 0:
            leg = ax1.legend()
        
        simCount +=1

plt.show()

