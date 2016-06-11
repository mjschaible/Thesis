import numpy as np
import scipy
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy import interpolate
import itertools
import re

def log_plots(runs, b):
    kcalToeV = 0.0433641 # convert kcal/mole to eV/atom

    a = len(runs)-1 #num columns
    #b =  num rows
#    print "a\t=\t%d\nb\t=\t%d\na*b\t=\t%d\nn\t=\t%d" % (a,b,a*b,n)
#    fig, axs = plt.subplots(b,a,figsize=(a*b,b*b), sharex='col', sharey='row')
#    axs.shape = (b,a)

    pemin=0
    pemax=0
    kemin=0
    kemax=0
    fig=plt.figure(0)
    ax=[]
    ax2=[]
    for i in range(1,len(runs)):
        time = (runs[i].step-runs[0].step[0])*runs[i].descrip[0]/1000
        num_molec = runs[i].descrip[2]
        ax.append(plt.subplot(b,a,i))
        ax2.append(ax[-1].twinx())
        #ax.plot(time,runs[i].pe*kcalToeV,label = runs[i].thermoCol[2])
        ax[-1].plot(time,runs[i].peave*kcalToeV*num_molec,color='blue',label='Sys. {} (eV)'.format(runs[i].thermoCol[3]))
        ax[-1].locator_params(axis='x', tight=True, nbins=2) 
        plt.setp(ax[-1].get_xticklabels(), visible=True)
        plt.setp(ax[-1].get_yticklabels(), visible=False)
        
        #ax2.plot(time,runs[i].temp,label = runs[i].thermoCol[1])
        ax2[-1].plot(time,runs[i].tempave,color='red',label='Sys. {} (K)'.format(runs[i].thermoCol[1]))
        ax2[-1].locator_params(axis='x', tight=True, nbins=2) 
        plt.setp(ax2[-1].get_yticklabels(), visible=False)
        #axs[1,i].plot(time,runs[i].ke*kcalToeV,label = runs[i].thermoCol[5])
        #axs[1,i].plot(time,runs[i].keave*kcalToeV*param[i].num_molec)
 
        pemin_t=min(runs[i].peave)*kcalToeV*num_molec
        pemax_t=max(runs[i].peave)*kcalToeV*num_molec
        kemin_t=min(runs[i].keave)*kcalToeV*num_molec
        kemax_t=max(runs[i].keave)*kcalToeV*num_molec
        if i == 1:
            pemax = pemax_t
            pemin = pemin_t
            ax[-1].set_ylim([pemin,pemax+10.5])
        elif pemin_t > pemin:
            pemax = pemin_t
            ax[-1].set_ylim([pemin,pemax+10.5])
        if i == len(runs)-1:
            kemax=kemax_t
            kemin=kemin_t
            ax2[-1].set_ylim([kemax-10.5,kemax+0.75])
            
        for ac in ax:
            ac.get_shared_y_axes().join(ac,ax[-1])
        for ac in ax2:
            ac.get_shared_y_axes().join(ac,ax2[-1])
        #axs[2,i-1].plot(time,runs[i].etot*kcalToeV,label = runs[i].thermoCol[7])
        #axs[2,i].set_xlabel('(ss={})'.format(runs[i].descrip[0]))
        #axs[2,i-1].ticklabel_format(axis='x', style = 'sci')
        #axs[2,i-1].locator_params(axis='x', tight=True, nbins=2)

        if i == 1:
            h1,l1=ax[-1].get_legend_handles_labels()
            h2,l2=ax2[-1].get_legend_handles_labels()
            ax[-1].legend(h1+h2, l1+l2, loc=2, fontsize=10)
            ax[-1].set_ylabel('Energy ($eV$)')
            #axs[2,i-1].set_ylabel('{} (eV)'.format(runs[i].thermoCol[7]))
            #plt.figtext(0.5,0.95,'The sim is {}'.format(runs[i].descrip),ha='center')
        plt.suptitle('Time [ps]',y=0.5,fontsize=16)
        plt.subplots_adjust(wspace=0.001)
        #axs[2,i].set_ylim(-1700,-1600)

    
    return

def msd_plots(runs, label, c):
    fig=plt.figure(1)
    for i in range(len(runs)):
        #time = (runs[i].step-runs[i].step[0])*param[i].timesteps/1000
        if len(runs[i].step) == len(runs[i].data):
            continue
            #axrow[i].plot(runs[i].step,runs[i].data,label = runs[i].Head[2])
        elif 'Shell Avg' in label:
            ax=plt.subplot(2,2,c)
            for j in range(len(runs[i].data)):
                lbl = '{} to {} $\AA$'.format(runs[i].Head[0]*j,runs[i].Head[0]*(j+1))
                ax.plot(runs[i].step,runs[i].data[j],label = lbl)
            if 'MSD' in label:
                ax.set_ylabel('{} ($\AA^2$)'.format(label))
            elif 'KE' in label:
                ax.legend(fontsize=10, loc=2)
                ax.set_ylabel('{} ($eV$)'.format(label))
                #ax.set_ylim([0,100])

        #plt.figtext(0.5,0.95,'The sim is {}'.format(runs[i].descrip),ha='center')
        ax.set_xlabel('Time [ps]')
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.ticklabel_format(axis='x', style = 'sci')
        #axrow[i].locator_params(axis='x', tight=True, nbins=2)
        
        return

def rdf_plots(runs):
    for i in range(len(runs)):
        fig, axs = plt.subplots(2,1, sharex='col')
        axs[0].plot(runs[i].pos,runs[i].RDF)
        axs[0].set_ylabel('$g_{OO}$(r)')
        axs[1].plot(runs[i].pos,runs[i].coordN)
        axs[1].set_ylabel('Coordination Number')
        axs[1].set_xlabel('Distance [A]')
        
        return

def msd2_plots(com, msd):
    plt.figure()
    #print range(len(msd))
    #print len(steparr), len(msd[0])
    for j in range(len(msd)):
        plt.plot(steparr,msd[j])
        
    return
