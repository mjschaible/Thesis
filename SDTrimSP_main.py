import numpy as np
import matplotlib.pyplot as plt
plt.interactive(True)
import glob

from Tkinter import Tk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename

import SDTrimSP_readSput
import SDTrimSP_plotSput

# ------------ Begin main program ----------
cont = 1
num_runs=0


# ----- Define a default array of simulation energies -----
defEng = np.linspace(100, 10000, num=100)
np.savetxt(
    'output.dat',
    defEng,
    fmt='%.0f',
    delimiter=',',
    newline='\n',
    header='energies')

expt_yld=[]

loglist=[]
log_yld=[]
sput_yld=[]
datlist=[]
trfile=[]

# read experimental data file
print "Please identify the experimental data file."
# show an "Open" dialog box and return the path to the selected file
#Tk().withdraw() 
#efn = askopenfilename()
efn='Expt_Data.dat'
expt_yld.append(SDTrimSP_readSput.read_exptfile(efn))
ion_target_pairs=expt_yld[0].label
for expt in range(len(ion_target_pairs)):
    plot_expt = SDTrimSP_plotSput.plot_sputExpt(expt_yld, expt, 'x')
#plt.show()

while cont != 0:
    Tk().withdraw() 
    path = askdirectory() 
    logfiles = path+"/*.log"
    datfiles = path+"/*.dat"                 

# first read log data files to get full simulation description
    for filename in glob.glob(logfiles):
#        print 'The file is {0}'.format(filename)
        loglist.append(filename)
        log_yld.append(SDTrimSP_readSput.read_logfile(filename)) 

# Plot total yields derived from the log files
    plot1_yld = SDTrimSP_plotSput.plot_log(log_yld, num_runs, '-', 'isbv')
    plt.show()

    # After all files in a folder are read, sort the yields from the log and sputter data files
    # by energy. 
#    sputteredF = [x for (y,x) in sorted(zip(Esim,sputteredF), key=lambda pair:pair[0])]
#    totYld = [x for (y,x) in sorted(zip(Esim,totYld), key=lambda pair:pair[0])]
#    Esim[j].sort()
#    print Esim
#    print fluenz
#    print sputteredF
#    print totYld

    # Read in sputter data and trajectory files
    # and arrange according to energy

    for filename in glob.glob(datfiles):
#        print 'The file is {0}'.format(filename)
        datlist.append(filename)        
    print
#    print datlist

# Call read functions for each data file type 
    for filename in datlist:
        if 'E0_31'in filename:
#            print 'The file is {0}'.format(filename)
            continue
           # File with variation in elemental concentration vs. timestep(+?)
        elif 'E0_33' in filename:
            # File with variation in elemetal sputter yield vs. timestep
            sput_yld.append(SDTrimSP_readSput.read_sputfile(filename))
        elif 'E0_34'in filename:
#            print 'The file is {0}'.format(filename)
            continue
            # File with energy loss moments?
        elif 'EngAn' in filename:
#            print 'The file is {0}'.format(filename)
            continue
            # File with analysis of how energy of projectile and target atoms is lost/partitioned
        elif 'layer' in filename:
#            print 'The file is {0}'.format(filename)
            continue
            # File with final(?) elemental variation as a function of depth
        elif 'tr_all' in filename:
            trfile.append(filename)
#            projectr.append(read the trajectory file)
#            print 'The file is {0}'.format(filename)
            # Trajectories for a number of projectiles and excited particles
        elif 'tr' in filename:
#            print 'The file is {0}'.format(filename)
            continue
            # Trajectories for primary and secondary knock-on atoms
        else:
            print 'The file {} does not have an analysis function written for it'.format(filename)

    index = [i for i, s in enumerate(ion_target_pairs) if ion_target_pairs[s] in sput_yld[num_runs].label[3]]
    print index

    # Determine average yields
    sput_yld_avg=SDTrimSP_readSput.average(sput_yld)
    # Extract the sputter yield from the first time step
    sput_yld_init=SDTrimSP_readSput.find_sputvar(sput_yld_avg, 1)
    # Extract the sputter yield from the last time step
    last= len(sput_yld_avg[0].fluence)-1
    sput_yld_final=SDTrimSP_readSput.find_sputvar(sput_yld_avg, last)
    #----- Plot total yields derived from the sputter data files -----
    plot_init_yld=SDTrimSP_plotSput.plot_log(sput_yld_init, nf, '--', None)
    plot_final_yld=SDTrimSP_plotSput.plot_log(sput_yld_final, nf, '-.', None)
#    plot_avg_yld=SDTrimSP_plotSput.plot_sput(sput_yld_init, num_runs)

    plot3_sput=SDTrimSP_plotSput.plot_flu(sput_yld_avg, ns+1)
    plt.show()

# Arrange energy in the default list np.linspace(100,10000)    
    # Define log file arrays
    sim_name=[None]*len(defEng)
    e_sort=[None]*len(defEng)
    yld_sort=[None]*len(defEng)

    ymax=0
    for j in range(len(log_yld)):
        for i in range(len(defEng)):
            if log_yld[j].energy==defEng[i]:
                if log_yld[j].label[2]=='isbv=1':
                    sim_name[i]=log_yld[j].label
                    e_sort[i]=log_yld[j].energy
                    yld_sort[i]=log_yld[j].totYld
                    if yld_sort[i]>ymax:
                        ymax=yld_sort[i]
                        maxi=i
            elif e_sort[i] != None:
                continue
            else:
                e_sort[i]=defEng[i]
                yld_sort[i]=0

#    for i in range(len(sim_name)):
#        print e_sort[i], yld_sort[i]

#       

    num_runs+=1
    cont = input('Enter 0 to end: ')

'''

    row0=[]
    allr=[]
    with open('output.dat', 'r') as csvinput:
        r = csv.reader(csvinput)
        row0 = next(r)
        row0.append([' '.join(simName[0])])
        allr.append(row0)
        i=0
        for row in r:
            row.append(yld_sort[i])
            allr.append(row)
            i+=1            
    with open('output.dat', 'w') as csvoutput:
        w = csv.writer(csvoutput, lineterminator='\n')
        w.writerows(allr)


'''





