import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
plt.interactive(True)
import glob
import itertools

from Tkinter import Tk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename

import SDTrimSP_readSput
import SDTrimSP_plotSput

def find_between( s, first, last=None ):
    try:
        start = s.index( first ) + len( first )
        return s[start:]
    except ValueError:
        return ""

def find_name(yld_label, color): 
    for i in range(len(ion_target_pairs)):
        if ion_target_pairs[i] in yld_label:
            nf = i
            c=next(color[i])
            return nf, c
        else:
            return None
    
# ------------ Begin main program ----------
cont = 1
nr=0

# ----- Define a default array of simulation energies -----
defEng = np.linspace(100, 10000, num=100)

expt_yld=[]
# read experimental data file
print "Please identify the experimental data file."
# show an "Open" dialog box and return the path to the selected file
root = Tk()
#root.withdraw() 
#efn = askopenfilename()
efn='ExpData.dat'
expt_yld.append(SDTrimSP_readSput.read_exptfile(efn))
ion_target_pairs=expt_yld[0].label

for expt in range(len(ion_target_pairs)-3):
    plot_expt = SDTrimSP_plotSput.plot_sputExpt(expt_yld, expt, 'o', 'Experimental')

color=[iter(plt.cm.rainbow(np.linspace(0,1,10))) for i in ion_target_pairs]
while cont != 0:
    loglist=[]
    log_yld=[]
    srim_yld=[]
    sput_yld=[]
    srimfn=[]
    trfile=[]
    n_eng=0
    print 'go go go'
    root.withdraw() 
    path = askdirectory() 
    logfiles = path+"/*.log"
    datfiles = path+"/*.dat"
    srimfiles = path+"/*.srim"
    path_name=find_between(path, 'data')
    print path_name
# first read log data files to get full simulation description
    for filename in glob.glob(logfiles):
#        print 'The file is {0}'.format(filename)
        loglist.append(filename)
        log_yld.append(SDTrimSP_readSput.read_logfile(filename))
        #print log_yld[n_eng].label[0], log_yld[n_eng].label[1]   
        n_eng+=1

    for filename in glob.glob(srimfiles):
        if 'Thiel' not in filename:
            srimfn.append(filename)
            #print filename
            srim_yld.append(SDTrimSP_readSput.read_exptfile(filename))

# Call read functions for each data file type 
    for filename in glob.glob(datfiles):
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
            #print 'The file {} does not have an analysis function written for it'.format(filename)
            continue
        
    #----- Plot total yeilds derived from SRIM simulations -----
    if not srim_yld:
        print "No srim yield files present"
    else:
        for i in range(len(srim_yld)):
            print len(srim_yld[i].label)
            for j in range(len(srim_yld[i].label)):
                nsim = srim_yld[i].label[j]
                #print nsim
                for k in range(len(ion_target_pairs)-3):
                    if ion_target_pairs[k] in nsim:
                        nf = k
                        c=next(color[k])
                #nf, c = find_name(nsim, color)
                plot_expt = SDTrimSP_plotSput.plot_srim(srim_yld,nf,'-.',c,nsim)

    #----- Plot total yields derived from the log files -----
    if not log_yld:
        print "No log files present"
    else:
        nsim = log_yld[0].label[0]
        for i in range(len(ion_target_pairs)):
            if ion_target_pairs[i] in nsim:
                nf = i
                c=next(color[i])
        #nf, c = find_name(nsim, color)
        #plot_yld = SDTrimSP_plotSput.plot_log(log_yld,nf,'-',c,path_name)

    if not sput_yld:
        print "No sput yield files present"
    else:
        nsim = sput_yld[0].label[0]
        for i in range(len(ion_target_pairs)):
            if ion_target_pairs[i] in nsim:
                nf = i
                c=next(color[i])
        # Determine average yields
        sput_yld_avg=SDTrimSP_readSput.average(sput_yld)
        # Extract the sputter yield from the first time step
        sput_yld_init=SDTrimSP_readSput.find_sputvar(sput_yld_avg, 1)
        # Extract the sputter yield from the last time step
        last= len(sput_yld_avg[0].fluence)-1
        sput_yld_final=SDTrimSP_readSput.find_sputvar(sput_yld_avg, last)
        
        #----- Plot total yields derived from the sputter data files -----
        plot_yld = SDTrimSP_plotSput.plot_log(log_yld,nf,'-',c,path_name)
        plot_init_yld=SDTrimSP_plotSput.plot_log(sput_yld_init,nf,'--',c,None)
        #    plot_final_yld=SDTrimSP_plotSput.plot_log(sput_yld_final, nf, '-.', None)
        #    plot_avg_yld=SDTrimSP_plotSput.plot_sput(sput_yld_init)

        plot3_sput=SDTrimSP_plotSput.plot_flu(sput_yld_avg, len(ion_target_pairs)+nf)
        #    plt.show()
    nr+=1
    cont = input('Enter 0 to end: ')

'''
np.savetxt(
    'output.dat',
    defEng,
    fmt='%.0f',
    delimiter=',',
    newline='\n',
    header='energies')

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

    # After all files in a folder are read, sort the yields from the log and sputter data files
    # by energy. 
#    sputteredF = [x for (y,x) in sorted(zip(Esim,sputteredF), key=lambda pair:pair[0])]
#    totYld = [x for (y,x) in sorted(zip(Esim,totYld), key=lambda pair:pair[0])]
#    Esim[j].sort()
#    print Esim
#    print fluenz
#    print sputteredF
#    print totYld
'''





