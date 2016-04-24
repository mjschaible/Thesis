import numpy as np
import scipy
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy import interpolate
import itertools
import re

def log_plots(runs, param):
    kcalToeV = 0.0433641 # convert kcal/mole to eV/atom

    a = len(runs)
    b = 3
#    print "a\t=\t%d\nb\t=\t%d\na*b\t=\t%d\nn\t=\t%d" % (a,b,a*b,n)
    fig, axs = plt.subplots(b,a,figsize=(a*b,b*b), sharex='col', sharey='row')
    axs.shape = (b,a)

    pemin=0
    pemax=0
    kemin=0
    kemax=0

    for i in range(len(runs)):
        time = (runs[i].step-runs[0].step[0])*param[i].timesteps/1000
        num_molec = re.findall('\d+\.\d+',runs[i].descrip[1])
        num_molec = float(num_molec[0])
        axs[0,i].plot(time,runs[i].pe*kcalToeV,label = runs[i].thermoCol[2])
        axs[0,i].plot(time,runs[i].peave*kcalToeV*num_molec)

        plt.setp(axs[0,i].get_xticklabels(), visible=False)
        axs[0,i].locator_params(axis='x', tight=True, nbins=2)

        axs[1,i].plot(time,runs[i].temp,label = runs[i].thermoCol[1])
        axs[1,i].plot(time,runs[i].tempave)
#        axs[1,i].plot(time,runs[i].ke*kcalToeV,label = runs[i].thermoCol[5])
#        axs[1,i].plot(time,runs[i].keave*kcalToeV*param[i].num_molec)
        plt.setp(axs[1,i].get_xticklabels(), visible=False)
        axs[1,i].locator_params(axis='x', tight=True, nbins=2)
 
        pemin_t=min(runs[i].peave)*kcalToeV*num_molec
        pemax_t=max(runs[i].peave)*kcalToeV*num_molec
        kemin_t=min(runs[i].keave)*kcalToeV*num_molec
        kemax_t=max(runs[i].keave)*kcalToeV*num_molec
        if i == 1:
            pemax = pemax_t
            pemin = pemin_t
            axs[0,i].set_ylim([pemin-0.5,pemax+5.5])
        elif pemin_t > pemin:
            pemax = pemin_t
            axs[0,i].set_ylim([pemin-0.5,pemax+5.5])
        if kemax_t > kemax:
            kemax = kemax_t
#            axs[1,i].set_ylim([kemax-5.5,kemax+0.5])

        axs[2,i].plot(time,runs[i].etot*kcalToeV,label = runs[i].thermoCol[7])
        axs[2,i].set_xlabel('(ss={})'.format(param[i].timesteps))
        axs[2,i].ticklabel_format(axis='x', style = 'sci')
        axs[2,i].locator_params(axis='x', tight=True, nbins=2)

        if i == 0:
            axs[0,i].set_ylabel('Sys. {} (eV)'.format(runs[i].thermoCol[3]))
            axs[1,i].set_ylabel('Sys. {} (K)'.format(runs[i].thermoCol[1]))
            axs[2,i].set_ylabel('{} (eV)'.format(runs[i].thermoCol[7]))
            plt.figtext(0.5,0.95,'The sim is {}'.format(runs[i].descrip),ha='center')
            plt.suptitle('Time [ps]',y=0.05,fontsize=16)
        if i > 0:
            plt.setp(axs[0,i].get_yticklabels(), visible=False)
            plt.setp(axs[1,i].get_yticklabels(), visible=False)
            plt.setp(axs[2,i].get_yticklabels(), visible=False)
        plt.subplots_adjust(wspace=0.001)
        #axs[2,i].set_ylim(-1700,-1600)

    return

def msd_plots(axrow, runs, param, label):
    for i in range(len(runs)):
        time = (runs[i].step-runs[i].step[0])*param[i].timesteps/1000

        axrow[i].plot(time,runs[i].data,label = runs[i].Head[2])
        axrow[i].set_ylabel(label)
        plt.figtext(0.5,0.95,'The sim is {}'.format(runs[i].descrip),ha='center')
        if i > 0:
            plt.setp(axrow[i].get_yticklabels(), visible=False)
            axrow[i].ticklabel_format(axis='x', style = 'sci')
            axrow[i].locator_params(axis='x', tight=True, nbins=2)
        if label=='PKA msd (A^2)':
            axrow[i].set_xlabel('Time [ps] (step size = {} fs)'.format(param[i].timesteps))
        
        return

def rdf_plots(runs, param):
    for i in range(len(runs)):
        
        fig, axs = plt.subplots(2,1, sharex='col')

        axs[0].plot(runs[i].pos,runs[i].RDF)
        axs[0].set_ylabel('$g_{OO}$(r)')
        axs[1].plot(runs[i].pos,runs[i].coordN)
        axs[1].set_ylabel('Coordination Number')
        
        axs[1].set_xlabel('Distance [A]')
        
        return

def msd2_plots(runs, msd, param):
    
    steparr=[runs[x].timestep for x in range(len(runs))]
#    print param[1].timesteps
#    time = (steparr[:]-steparr[0])*param[1].timesteps/1000
    plt.figure()
    #print range(len(msd))
    for j in range(len(msd)):
        plt.plot(steparr,msd[j])
        
    return
