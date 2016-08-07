import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
plt.interactive(True)
import glob
import itertools

from Tkinter import Tk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename

import SDTrimSP_readSput
import SDTrimSP_plotSput

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.titlesize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['axes.labelpad'] = 2.5
mpl.rcParams['xtick.labelsize']='large'
mpl.rcParams['ytick.labelsize']='large'
mpl.rcParams['legend.handlelength']=2.5

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
nr=1
incl_expt=0
incl_srim=0
elem_comp=2
expt_comp=0

# Variable to determine whether or not to plot yields vs. fluence for specified energies
fyld=1

root = Tk()
# ----- Import and plot experimental data -----
if incl_expt==1:
    #print "Please identify the experimental data file."
    # show an "Open" dialog box and return the path to the selected file
    #root.withdraw() 
    #efn = askopenfilename()
    efn='ExpData.csv'
    expt_yld=SDTrimSP_readSput.read_exptfile(efn)
    ion_target_pairs=list(set(expt_yld.label[0]))
    target= list(set(expt_yld.label[2]))
    #print target

    plot_expt = SDTrimSP_plotSput.plot_sputExpt(expt_yld, ion_target_pairs)

# ----- Import and plot SRIM data -----
if incl_srim==1:
    srim_yld=[]
    srimfn=[]
    srimfiles = "./SRIM/*.srim"
    for filename in glob.glob(srimfiles):
        if 'Thiel' not in filename:
            srimfn.append(filename)
            #print filename
            srim_yld.append(SDTrimSP_readSput.read_exptfile(filename))
            met_class='srim'
    #----- Plot total yields derived from SRIM simulations ----
    for i in range(len(srim_yld)):
        plot_srim = SDTrimSP_plotSput.plot_srim(srim_yld[i],ion_target_pairs)

if elem_comp==1:
    color=[iter(plt.cm.viridis(np.linspace(0,1,6))) for i in ion_target_pairs]
if elem_comp==2:
    MetClass = ["Lunar", "MetStds", "Mars", "HEDs", "Aubrites", "Urelites", "CCs", "OCs"] #"LunarAnalog",
    color=iter(['orange','blue','red','green','cyan','purple','grey','black'])
    #color=iter(plt.cm.Set1(np.linspace(0,1,len(MetClass))))

# ---- Begin procedure to read in SDTrimSP output files ----
while cont != 0:
    loglist=[]
    log_yld=[]
    out_yld=[]
    ER=[]
    sput_yld=[]
    trfile=[]
    n_eng=0

    root.withdraw() 
    path = askdirectory() 
    logfiles = path+"/*.log"
    outfiles = path+"/*.out"
    datfiles = path+"/*.dat"
    path_name=find_between(path, 'data')
    #print path_name

    # first read log data files to get full simulation description
    for filename in glob.glob(logfiles):
        #print 'The file is {0}'.format(filename)
        loglist.append(filename)
        log_yld.append(SDTrimSP_readSput.read_logfile(filename))
        #print log_yld[-1].label[4]
        n_eng+=1
    log_yld=sorted(log_yld, key=lambda x: x.label[0])#, reverse=True)
    tar=[i.label[0] for i in log_yld]
    targets=list(set([i.label[0].split('->')[1] for i in log_yld]))

    # Print the list of targets used in simulations for the current directory
    # List is taken from the *.log data files
    #print targets

    # read in the output files which contain the sputtering yield vs. fluence data
    for filename in glob.glob(outfiles):
        #print 'The file is {0}'.format(filename)
        out_yld.append(SDTrimSP_readSput.read_outfile(filename))
    #out_yld=sorted(log_yld, key=lambda x: x.label[0])
    rat=[]

    for i in range(len(out_yld)):
        out_yld[i]=sorted(out_yld[i], key=lambda x: x.label[1])
        for j in range(len(out_yld[i])):
            ctar = out_yld[i][j].label[1]
            #print ctar
            for n, t in enumerate(tar):
                if t == ctar:
                    index=n
            rat.append(out_yld[i][j].label[1])
            #print len(out_yld[i][j].Flux)
            #print log_yld[index].label[4]
            if len(out_yld[i][j].Flux) == len(log_yld[index].label[4]):
                out_yld[i][j].label.append(log_yld[index].label[4]) # append element names
                out_yld[i][j].label.append(log_yld[index].label[5]) # append atomic masses
                out_yld[i][j].label.append(log_yld[index].label[6]) # append SBE
            else:
                print "the number of elements does not match"
                print len(out_yld[i][j].Flux)
                print len(log_yld[nlog].label[4])
            met_class=out_yld[i][j].label[0]

    if rat!=tar:
        print "ERROR: Out and Log file targets don't match"
        #print out_yld[0][-1].energy
        #print log_yld[0].energy
        
    #----- Plot total yields derived from the *.out files -----
    if out_yld:
        for out in out_yld:
            # Determine average yield
            out_yld_avg=SDTrimSP_readSput.average(out)
        nsim = out_yld[0][0].label[0]
        #print nsim
        
        if expt_comp==1:
            for j, k  in enumerate(ion_target_pairs):
                #print j, k
                if k in nsim:
                    nf = j
                    c=next(color[j])
            plot3_sput=SDTrimSP_plotSput.plot_out1(out_yld,nf,'-',c)

    if expt_comp==1:
        for i in range(len(ion_target_pairs)):
            if ion_target_pairs[i] in nsim:
                nr = i+1
                c=next(color[i])
                lbl=path_name
            else:
                nr=nf+10
        # plot_out1 plots the total yield vs. energy for SDTrimSP simulatioXns
        plot_log = SDTrimSP_plotSput.plot_log(log_yld,nr,'-',c,lbl)
        
    # ---- Plot relative elemental yield comparisons ----
    if elem_comp == 2:
        c=next(color)
        nf=0
        ER=SDTrimSP_readSput.comp_yield(out_yld,targets)
        ER_plot=SDTrimSP_plotSput.plot_ER(ER, c, nf, met_class, targets)
        #OvC_plot=SDTrimpSP_plotSput.plot_YvO()

    # ---- Plot the species yields vs. fluence ----
    if fyld==1:
        if elem_comp==2:
            tag='Met'
            nf=1
        elif elem_comp==1:
            tag='SOx'
        for out in out_yld:
            plot_out=SDTrimSP_plotSput.plot_out2(out,nf,c,met_class)
            #ploty_sput=SDTrimSP_plotSput.plot_out3(out,met_class,tag)

    #----- Plot total yields derived from the log files -----
    if not log_yld:
        print "No log files present"
    else:
        if expt_comp==2:
            neut_yld, ion_yld=SDTrimSP_readSput.log_comp(log_yld,targets)
            nsim = log_yld[0].label[0]
            nf=nr+1
            lbl='bar'
            plot_log = SDTrimSP_plotSput.plot_log(neut_yld,nf,'-',c,lbl)
            plot_log = SDTrimSP_plotSput.plot_log(ion_yld,nf,'-',c,lbl,'//',met_class)

    nr+=1
    cont = input('Enter 0 to end: ')

    '''
# Call read functions for each data file type 
    for filename in glob.glob(datfiles):
        if 'E0_31'in filename:
#            print 'The file is {0}'.format(filename)
            continue
           # File with variation in elemental concentration vs. timestep(+?)
        elif 'E0_33' in filename:
            # File with variation in elemetal sputter yield vs. timestep
            #print filename
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





