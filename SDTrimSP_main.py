import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
plt.interactive(True)
import glob
import itertools
import os.path

from Tkinter import Tk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename

import SDTrimSP_readSput
import SDTrimSP_plotSput

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.titlesize'] = 'x-large'
mpl.rcParams['axes.labelsize'] = 'x-large'
#mpl.rcParams['axes.labelpad'] = 2.5
mpl.rcParams['xtick.labelsize']='large'
mpl.rcParams['ytick.labelsize']='large'
mpl.rcParams['legend.fontsize']='large'
mpl.rcParams['legend.handlelength']=2.5
mpl.rcParams['lines.markersize']=8

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
# Specify which type of target is being analyzed
# 0 == single system; 1 == simple oxides; 2 == meteorites
elem_comp=2

if elem_comp==1:
    genClass = ["Mars","Aubrites"]
    tarClass = ["Mars","Aubrites"]
    color=iter(plt.cm.Set1(np.linspace(0,1,len(MetClass))))
elif elem_comp==1:
    genClass = ["H_SiO2", "H_Al2O3", "He_SiO2", "He_Al2O3"]
    color=[iter(plt.cm.viridis(np.linspace(.2,.8,3))) for i in genClass]
    tarClass = [ "SBEO1eV_nodiff", "SBEO2eV_nodiff", "SBEO3eV_nodiff"]
    marker=[iter(['v','^','<','>','v','^','<','>']) for i in genClass]
    mfac=[iter([None,None,None,None,'None','None','None','None']) for i in genClass]
elif elem_comp==2:
#    MetClass = ["Lunar","HEDs","Mars","Aubrites","Urelites","CCs","OCsECs"] #"LunarAnalog"
    genClass = ["Lunar","HEDs","Mars","Aubrites","Urelites","OCsECs","CCs"]
    tarClass = ["Lunar","HEDs","Mars","Aubrites","Urelites","OCsECs","CCs"]
    color=iter(['grey','blue','red','green','purple','orange','black'])
    marker=iter(['o', 's', 'D', '^', 'v', '<', '>'])
    #color=iter(plt.cm.Set1(np.linspace(0,1,len(MetClass))))

# Specify whether or not to include experimental and SRIM results
# 1 == yes; 0 == no
incl_expt=0
incl_srim=0

# Specify how to compare the experiments (?)
expt_comp=2 # 1= SOx, 2=Met
shift=-0.3

# Variable to determine whether or not to plot yields vs. fluence for specified energies
# no == 0, yes == 1
fyld=0

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
    plot_expt = SDTrimSP_plotSput.plot_sputExpt(expt_yld, genClass)

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
        plot_srim = SDTrimSP_plotSput.plot_srim(srim_yld[i],genClass)

# ---- Begin procedure to read in SDTrimSP output files ----
root = Tk()
root.withdraw() 
path = askdirectory()
dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

for cdir in genClass:
    print cdir
    ddirs=[d for d in os.listdir(cdir) if os.path.isdir(os.path.join(cdir, d))]
    ddirs.sort()
    for adir in ddirs:
        if cdir in genClass and adir in tarClass:
            loglist=[]
            log_yld=[]
            out_yld=[]
            ER=[]
            sput_yld=[]
            trfile=[]
            n_eng=0

            logfiles = path+'/'+cdir+'/'+adir+"/*.log"
            outfiles = path+'/'+cdir+'/'+adir+"/*.out"
            datfiles = path+'/'+cdir+'/'+adir+"/*.dat"
            path_name=find_between(path, 'data')
            print path_name

            # first read log data files to get full simulation description
            for filename in glob.glob(logfiles):
                print 'The file is {0}'.format(filename)
                loglist.append(filename)
                log_yld.append(SDTrimSP_readSput.read_logfile(filename))
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
                    for n, t in enumerate(tar):
                        if t == ctar:
                            index=n
                    rat.append(out_yld[i][j].label[1])
                    if len(out_yld[i][j].Flux) == len(log_yld[index].label[4]):
                        out_yld[i][j].label.append(log_yld[index].label[4]) # append element names
                        out_yld[i][j].label.append(log_yld[index].label[5]) # append atomic masses
                        out_yld[i][j].label.append(log_yld[index].label[6]) # append SBE
                        #if out_yld[i][j].label[6][-1]==2.0:
                        #print out_yld[i][j].label[6][-1]
                    else:
                        print "the number of elements does not match"
                        print len(out_yld[i][j].Flux)
                        print len(log_yld[nlog].label[4])
                    met_class=out_yld[i][j].label[0]
                    #print met_class

#----- Plot total yields derived from the *.out files -----
            if out_yld:
                for out in out_yld:
                    # Determine average yield
                    out_yld_avg=SDTrimSP_readSput.average(out)
                nsim = out_yld[0][0].label[0]

                if expt_comp==1:
                    for j, k  in enumerate(genClass):
                        k=k.replace('_', '->', 1)
                        if k in nsim:
                            nf = j+1
                            c=next(color[j])
                            mk=next(marker[j])
                            mfc=next(mfac[j])
                            lbl=path_name
                            # plot_out1 is the total yield vs. energy for SDTrimSP simulations
                            plot_sput=SDTrimSP_plotSput.plot_out1(out_yld,nf,c,mk,mfc)
                        else:
                            nf=-1
                else:
                    nf=-1

                # ---- Plot relative elemental yield comparisons ----
                if elem_comp==2:
                    c=next(color)
                    mk=next(marker)
                    tag='Met'
                    ER,sw_out,ion_out=SDTrimSP_readSput.comp_yield(out_yld,targets)
                    shift+=0.1
                    nf=1
                    #ploty_sput=SDTrimSP_plotSput.plot_iavg(sw_out,nf,c,shift,mk,adir)
                    nf=2
                    mean, std=SDTrimSP_plotSput.plot_iavg(ion_out,nf,c,shift,mk,adir)
                    nf=3
                    #ER_plot=SDTrimSP_plotSput.plot_ER(ER,c,met_class,adir,nf,mk)
                    nf=4
                    flux_plot=SDTrimSP_plotSput.plot_flux(mean,std,adir,nf)
                elif elem_comp==1:
                    tag='SOx'
                    nf=nr+10
                    #ploty_sput=SDTrimSP_plotSput.plot_iavg(out_yld,nf,c,shift,met_class)

            # ---- Plot the species yields vs. fluence ----
            if fyld==1:
                if elem_comp==1:
                    for out in out_yld:
                        ploty_sput=SDTrimSP_plotSput.plot_vflu(out,tag)
                elif elem_comp==2:
                    plot_sput=SDTrimSP_plotSput.plot_vflu(ion_out,tag)
                    for out in out_yld:
                        plot_out=SDTrimSP_plotSput.plot_vflu(sw_out,tag)

            #----- Plot total yields derived from the log files -----
            if log_yld:
                if expt_comp==2:
                    neut_yld, ion_yld=SDTrimSP_readSput.log_comp(log_yld,targets)
                    nsim = log_yld[0].label[0]
                    nf=nr+1
                    lbl='bar'
                    #plot_log = SDTrimSP_plotSput.plot_log(neut_yld,nf,'-',c,lbl)
                    #plot_log = SDTrimSP_plotSput.plot_log(ion_yld,nf,'-',c,lbl,'//',met_class)
                    #ploty_sput=SDTrimSP_plotSput.plot_out3(ion_yld,met_class,tag)

            plt.show()
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





