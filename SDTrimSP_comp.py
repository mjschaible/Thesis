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
from operator import add

import SDTrimSP_readSput
import SDTrimSP_plotSput

def find_between( s, first, last=None ):
    try:
        start = s.index( first ) + len( first )
        return s[start:]
    except ValueError:
        return ""
    
# ------------ Begin main program ----------
cont = 1
nr=0

# ----- Define a default array of simulation energies -----
defEng = np.linspace(100, 10000, num=100)
root = Tk()
root.withdraw()

elemlist = ['Na','Mg','Al','Si','K','Ca','Ti','Fe','Ni','O']
ifrac = 0.05
relcorr = [23.319, 2.886, 3.637, 1.000, 81.080, 3.950, 3.582, 2.767, 2.886, 0.01]

while cont != 0:
    tarlist=[]
    log_yld=[]
    sput_yld=[]
    esput=[]
    trfile=[]
    n_eng=0
 
    path = askdirectory() 
    logfiles = path+"/*.log"
    datfiles = path+"/*.dat"

    path_name=find_between(path, 'data')
    #print path_name
    # first read log data files to get full simulation description
    for filename in glob.glob(logfiles):
        #print 'The file is {0}'.format(filename)
        log_yld.append(SDTrimSP_readSput.read_logfile(filename))
        if log_yld[n_eng].label[4] in tarlist:
            continue
        else:
            tarlist.append(log_yld[n_eng].label[4])
        n_eng+=1
    print tarlist
    
# Call read functions for each data file type 
    for filename in glob.glob(datfiles):
        if 'E0_31'in filename:
            continue
           # File with variation in elemental concentration vs. timestep(+?)
        elif 'E0_33' in filename:
            # File with variation in elemetal sputter yield vs. timestep
            #print filename
            sput_yld.append(SDTrimSP_readSput.read_sputfile(filename))
        elif 'E0_34'in filename:
            continue
            # File with energy loss moments?
        elif 'EngAn' in filename:
            continue
            # File with analysis of how energy of projectile and target atoms is lost/partitioned
        elif 'layer' in filename:
            continue
            # File with final(?) elemental variation as a function of depth
        elif 'tr_all' in filename:
            trfile.append(filename)
#            projectr.append(read the trajectory file)
            # Trajectories for a number of projectiles and excited particles
        elif 'tr' in filename:
            continue
            # Trajectories for primary and secondary knock-on atoms
        else:
            continue

    for k in range(len(tarlist)):
        hyld=[]
        h_iyld=[]
        heyld=[]
        he_iyld=[]
        for i in range(len(log_yld)):
            print log_yld[i].label[3]
            print log_yld[i].Flux

            if log_yld[i].label[3][0] == 'H' and log_yld[i].label[4] == tarlist[k]:
                print 'The {} yield is {:.3f}'.format(log_yld[i].label[0], log_yld[i].totYld)
                for y in range(len(elemlist)):
                    match = [j for j, x in enumerate(log_yld[i].label[3]) if x == elemlist[y]]
                    if match:
                        yld=0.95*log_yld[i].Flux[match[0]]
                        hyld.append(yld)
                        h_iyld.append(yld*relcorr[y]*ifrac)
                print hyld
            elif log_yld[i].label[3][0] == 'He' and log_yld[i].label[4] == tarlist[k]:
                print 'The {} yield is {:.3f}'.format(log_yld[i].label[0], log_yld[i].totYld)
                for y in range(len(elemlist)):
                    match = [j for j, x in enumerate(log_yld[i].label[3]) if x == elemlist[y]]
                    if match:
                        yld=0.05*log_yld[i].Flux[match[0]]
                        heyld.append(yld)
                        he_iyld.append(yld*relcorr[y]*ifrac)
                print heyld
            else:
                print 'What ion are you using dumbass?'

        sw_yld= map(add,hyld,heyld)
        sw_iyld= map(add,h_iyld,he_iyld)
        print sw_iyld
            
    nr+=1
    cont = input('Enter 0 to end: ')

color=iter(plt.cm.rainbow(np.linspace(0,1,nr)))

'''
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
        #nf = 1
        #c=next(color[nf])
        plot_yld = SDTrimSP_plotSput.plot_log(log_yld,nf,'-',c,path_name)

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
        #plot_yld = SDTrimSP_plotSput.plot_log(log_yld,nf,'-',c,path_name)
        #plot_init_yld=SDTrimSP_plotSput.plot_log(sput_yld_init,nf,'--',c,None)
        #    plot_final_yld=SDTrimSP_plotSput.plot_log(sput_yld_final, nf, '-.', None)
        #    plot_avg_yld=SDTrimSP_plotSput.plot_sput(sput_yld_init)
        #plot3_sput=SDTrimSP_plotSput.plot_flu(sput_yld,len(ion_target_pairs)+nf,'-',1)
        #plot3_sput=SDTrimSP_plotSput.plot_flu(sput_yld,'-',1)
        #plot3_sput=SDTrimSP_plotSput.plot_flu(sput_yld_avg,nr, ':')
        #    plt.show()
    
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





