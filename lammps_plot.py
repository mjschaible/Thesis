import numpy as np
import scipy
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy import interpolate
import itertools
import re

def log_plots(runs, b, ft):
    kcalToeV = 0.0433641 # convert kcal/mole to eV/atom

    a = len(runs) #num columns

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
    if b>1:
        ax3=[]
        
    for i in range(ft,len(runs)):
        time = (runs[i].step-runs[0].step[0])*runs[i].descrip[0]/1000
        num_molec = runs[i].descrip[2]

        if i==0:
            ax.append(plt.subplot(b,a,i))
        else:
            ax.append(plt.subplot(b,a,i, sharey=ax[-1]))
        #ax[-1].plot(time,runs[i].pe*kcalToeV,label = runs[i].thermoCol[2])
        ax[-1].plot(time,runs[i].peave*kcalToeV*num_molec,color='blue',label='Sys. {} (eV)'.format(runs[i].thermoCol[3]))
        ax[-1].locator_params(axis='x', tight=True, nbins=2) 
        plt.setp(ax[-1].get_xticklabels(), visible=False)
             
        ax2.append(ax[-1].twinx())
        ax2[-1].plot(time,runs[i].tempave,color='red',label='Sys. {} (K)'.format(runs[i].thermoCol[1]))
        ax2[-1].locator_params(axis='x', tight=True, nbins=2)
        plt.setp(ax2[-1].get_xticklabels(), visible=False)
                
            
        if i>0:
            plt.setp(ax2[-1].get_yticklabels(), visible=False)
            plt.setp(ax[-1].get_yticklabels(), visible=False)

        try:
            ax3
            if i==0:
                ax3.append(plt.subplot(b,a,i+1+a))
            else:
                ax3.append(plt.subplot(b,a,i+1+a, sharey=ax3[-1]))
            ax3[-1].plot(time,runs[i].etot*kcalToeV,label = runs[i].thermoCol[7])
            #ax3[-2].set_xlabel('(ss={})'.format(runs[i].descrip[0]))
            ax3[-1].ticklabel_format(axis='x', style = 'sci')
            ax3[-1].locator_params(axis='x', tight=True, nbins=2)
            if i>0:
                plt.setp(ax3[-1].get_yticklabels(), visible=False)
        except NameError:
            print 'Only one axis defned'
                        
        pemin_t=min(runs[i].peave)*kcalToeV*num_molec
        pemax_t=max(runs[i].peave)*kcalToeV*num_molec
        #kemin_t=min(runs[i].keave)*kcalToeV*num_molec
        #kemax_t=max(runs[i].keave)*kcalToeV*num_molec
        kemin_t=min(runs[i].tempave)
        kemax_t=max(runs[i].tempave)
        if i == 1:
            pemax = pemax_t
            pemin = pemin_t
            ax[-1].set_ylim([pemin-5,pemax+5])
        elif pemin_t > pemin:
            pemax = pemin_t
            ax[-1].set_ylim([pemin-5,pemax+5])
        if i == len(runs)-1:
        #if kemax_t>kemax:
            kemax=kemax_t
            #kemin=kemin_t
            #ax2[-1].set_ylim([kemax-10.5,kemax+0.75])
            ax2[-1].set_ylim([kemin-5,kemax+5])

        #for ac in ax:
        #    ac.get_shared_y_axes().join(ac,ax[-1])
        for ac in ax2:
            ac.get_shared_y_axes().join(ac,ax2[-1])

        if i == 0:
            h1,l1=ax[-1].get_legend_handles_labels()
            h2,l2=ax2[-1].get_legend_handles_labels()
            #ax[-1].legend(h1+h2, l1+l2, loc=2, fontsize=10)
            ax[-1].set_ylabel('Energy ($eV$)')
            ax2[-1].set_ylabel('Temperature (K)')
            ax3[-1].set_ylabel('{} (eV)'.format(runs[i].thermoCol[7]))
            #plt.figtext(0.5,0.95,'The sim is {}'.format(runs[i].descrip),ha='center')
        plt.suptitle('Time [ps]',fontsize=16)
        plt.subplots_adjust(wspace=0.001)
        #axs[2,i].set_ylim(-1700,-1600)
    
    return

def msd_plots(runs, label, nr, c):
    fig = plt.figure(1)
    for i in range(len(runs)):
        #time = (runs[i].step-runs[i].step[0])*param[i].timesteps/1000
        if 'Shell Avg' in label:
            ax=fig.add_subplot(nr,2,c)
            for j in range(len(runs[i].data)):
                #print len(runs[i].step), len(runs[i].data[j])
                lbl = '{} to {} $\AA$'.format(runs[i].Head[0]*j,runs[i].Head[0]*(j+1))
                ax.plot(runs[i].step,runs[i].data[j],label = lbl)
            if 'MSD' in label:
                ax.set_ylabel('{} ($\AA^2$)'.format(label))
            elif 'KE' in label:
                ax.legend(fontsize=10, loc=1)
                ax.set_ylabel('{} ($eV$)'.format(label))
                ax.locator_params(axis='x', tight=True, nbins=10)
                #ax.set_ylim([0,100])
        else:
            ax=fig.add_subplot(nr,1,c)
            ax.plot(runs[i].step,runs[i].data,label = runs[i].Head[2])
            ax.set_ylabel(label)
            ax.locator_params(axis='x', tight=True, nbins=10)
            ax.ticklabel_format(axis='x', style = 'sci')

        if c==1:
            plt.setp(ax.get_xticklabels(), visible=True)
                
        #plt.figtext(0.5,0.95,'The sim is {}'.format(runs[i].descrip),ha='center')
    fig.tight_layout()
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
