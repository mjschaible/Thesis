import numpy as np
import matplotlib.pyplot as plt
import glob
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
'''
class LogData (object):

    def __init__(self,potential,step,temp,tempave,pe,peave,ke,keave,etot):
        self.potential=potential
        self.step=step
        self.temp=temp
        self.tempave=tempave
        self.pe=pe
        self.peave=peave
        self.ke=ke
        self.keave=keave
        self.etot=etot
'''

# ------------ Begin main program ----------
cont = 1
num_runs=0
runs=[]


while cont != 0:
#    name = raw_input('Enter filename: ')
#    filename.append(name)

    energeq = []
    Flux1eq = []
    Flux2eq = []
    Flux3eq = []    
    energies = []
    Flux1 = []
    Flux2 = []
    Flux3 = []
    simCount = 0

    if num_runs==0:
        print "Please identify the experimental data file."
        Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
        filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

        exptName, exptEng, exptFlux_HeSiO2 = read_exptfile(filename)
        b = plot_sputExpt(exptName, exptEng, exptFlux_HeSiO2, simCount)
        simCount +=1
    # open, read only ('r') file. Pull out sim version, description, and parameters
    #with open(filename, 'r') as logfile:
    else:
        Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
        path = askdirectory() # show an "Open" dialog box and return the path to the selected file

    ## First open and scan the file to determine the total number of lines, where the header lines end, and where the footer lines begin                
        scanfiles = path+"/*.dat"
        print scanfiles
        for filename in glob.glob(scanfiles):
        #    if os.path.isfile(fn):
            print 'The file is {0}'.format(filename)

            if 'E0_33' in filename:
                descrip, species, fsteps, specFlu1, specFlu2, specFlu3 = read_sputfile(filename)

                params = descrip.split()
                simEng = float(params[0])
                #print specFlu1
                numave = 20

                avYld0, av0 = movingaverage(specFlu1,numave)
                avYld1, av1 = movingaverage(specFlu2,numave)
                avYld2, av2 = movingaverage(specFlu3,numave)

                energeq.append(simEng)
                Flux1eq.append(avYld0)
                Flux2eq.append(avYld1)
                Flux3eq.append(avYld2)

                a = plot_sputFile(descrip, species, fsteps, specFlu1, specFlu2, specFlu3, av0, av1, av2)

            elif filename.find('E0') == -1:
                fluence, energy, elemName, sputteredF = read_logfile(filename)

                name = energy + 'eV'
                #       print sputteredF
                energies.append(float(energy))
                Flux1.append(float(sputteredF[0])/float(fluence))
                Flux2.append(float(sputteredF[1])/float(fluence))
                Flux3.append(float(sputteredF[2])/float(fluence))
                #print Flux1, energies

                #print Flux1_sorted, energies_sorted

            else:
                continue
    
        Flux1eq_sorted = [x for (y,x) in sorted(zip(energeq,Flux1eq), key=lambda pair:pair[0])]
        Flux2eq_sorted = [x for (y,x) in sorted(zip(energeq,Flux2eq), key=lambda pair:pair[0])]
        Flux3eq_sorted = [x for (y,x) in sorted(zip(energeq,Flux3eq), key=lambda pair:pair[0])]
        energies_sortedeq=energeq.sort()

        Flux1_sorted = [x for (y,x) in sorted(zip(energies,Flux1), key=lambda pair:pair[0])]
        Flux2_sorted = [x for (y,x) in sorted(zip(energies,Flux2), key=lambda pair:pair[0])]
        Flux3_sorted = [x for (y,x) in sorted(zip(energies,Flux3), key=lambda pair:pair[0])]
        energies_sorted=energies.sort()

        b = plot_sputYld(elemName, energies, Flux1_sorted, Flux2_sorted, Flux3_sorted, simCount)
        simCount +=1 
        c = plot_sputYld(species, energeq, Flux1eq_sorted, Flux2eq_sorted, Flux3eq_sorted, simCount)

    num_runs+=1
    cont = input('Enter 0 to end: ')

plt.show()





