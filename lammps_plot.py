import numpy as np
import scipy
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy import interpolate
import itertools
import re
import os

def log_plots(runs, row, ft):
    # ft is the 'filetype', either rad or eq or...
    # row is the number of rows
    # runs is the data object
    # Determine the filename for figure naming
    fn = runs[0].descrip[0]
    print fn

    kcalToeV = 0.0433641 # convert kcal/mole to eV/atom

    col = len(runs) #num columns

    #print "a\t=\t%d\nb\t=\t%d\na*b\t=\t%d\nn\t=\t%d" % (a,b,a*b,n)
    #fig, axs = plt.subplots(b,a,figsize=(a,a/b), sharex='col', sharey='row')
    #axs.shape = (b,a)

    pemin=0
    pemax=0
    kemin=0
    kemax=0
    #print ft
    if ft==0:
        fig=plt.figure(figsize=(col+1,col/row))
    elif ft==1:
        col = col-1
        fig=plt.figure(figsize=(3*col,3*row))
    ax=[]
    ax2=[]
    if row>1:
        ax3=[]

    for i in range(ft,len(runs)):
        if ft == 1:
            ncol = i-1
        else:
            ncol = i
        time = (runs[i].step-runs[0].step[0])*runs[i].descrip[1]/1000
        num_molec = runs[i].descrip[3]
        
        if i==ft:
            ax.append(plt.subplot(row,col,ncol+1))
        else:
            ax.append(plt.subplot(row,col,ncol+1, sharey=ax[0]))
            
        ax[-1].plot(time,runs[i].peave*kcalToeV*num_molec,color='blue',label='Sys. {} (eV)'.format(runs[i].thermoCol[3]))
        ax[-1].locator_params(axis='x', tight=True) 
        #plt.setp(ax[-1].get_xticklabels(), visible=False)
        plt.setp(ax[-1], xticks=[min(time),max(time)])

        if ft == 0:
            eng=runs[i].tempave
            label='Sys. {} (K)'.format(runs[i].thermoCol[1])
            kemin_t=min(runs[i].tempave)
            kemax_t=max(runs[i].tempave)
            units='Temperature (K)'
        if ft == 1:
            eng=runs[i].ke*kcalToeV
            label='Sys. {} (eV)'.format(runs[i].thermoCol[5])
            kemin_t=min(runs[i].ke)*kcalToeV
            kemax_t=max(runs[i].ke)*kcalToeV
            units='Kin. Energy ($eV$)'
        ax2.append(ax[-1].twinx())            
        ax2[-1].plot(time,eng,color='red',label=label)            
        plt.setp(ax2[-1].get_xticklabels(), visible=False)

        if i == ft:
            plt.setp(ax2[-1].get_yticklabels(), visible=False)
        elif i>ft:
            plt.setp(ax2[-1].get_yticklabels(), visible=False)
            plt.setp(ax[-1].get_yticklabels(), visible=False)

        if ncol==col-1:
            ax[0].set_ylabel('Pot. Energy ($eV$)')
            ax2[-1].set_ylabel(units)
            h1,l1=ax[-1].get_legend_handles_labels()
            h2,l2=ax2[-1].get_legend_handles_labels()
            l = ax[-1].legend(h1+h2, l1+l2, loc=4, fontsize=10)
            #l.set_zorder(20)
            plt.setp(ax2[-1].get_yticklabels(), visible=True)

        if i>ft:
            for ac in ax:
                ac.get_shared_y_axes().join(ac,ax[-1])
            for ac in ax2:
                ac.get_shared_y_axes().join(ac,ax2[-1])

        try:
            ax3
            if i==ft:
                ax3.append(plt.subplot(row,col,ncol+1+col))
                ax3[-1].set_ylabel('{} (eV)'.format(runs[i].thermoCol[7]))
            else:
                ax3.append(plt.subplot(row,col,ncol+1+col, sharey=ax3[-1]))
            ax3[-1].plot(time,runs[i].etot*kcalToeV,label = runs[i].thermoCol[7])
            #ax3[-2].set_xlabel('(ss={})'.format(runs[i].descrip[0]))
            ax3[-1].ticklabel_format(axis='x', style = 'sci')
            ax3[-1].locator_params(axis='x', tight=True, nbins=2)
            plt.setp(ax3[-1], xticks=[min(time), max(time)])
            if i>ft:
                plt.setp(ax3[-1].get_yticklabels(), visible=False)

        except NameError:
            print 'Only one axis defned'
                        
        pemin_t=min(runs[i].peave)*kcalToeV*num_molec
        pemax_t=max(runs[i].peave)*kcalToeV*num_molec
        
        if i == ft:
            pemax = pemax_t
            pemin = pemin_t
            kemax=kemax_t
            kemin=kemin_t
            ax[-1].set_ylim([pemin-5,pemax+5])
            #ax2[-1].set_ylim([kemin-5,kemax+5])
        if ft == 0:
            if kemax_t>kemax:
                kemax=kemax_t
                #kemin=kemin_t
                ax2[-1].set_ylim([kemin-5,kemax+5])
            if pemin_t > pemin:
                pemax = pemin_t
                ax[-1].set_ylim([pemin-5,pemax+5])
        elif ft == 1:
            if i == len(runs)-1:
                pemax=pemin_t
                kemax=kemin_t
                #kemin=kemin_t
                ax2[-1].set_ylim([kemax-5,kemax+5])
                ax[-1].set_ylim([pemax-5,pemax+5])
            
            #plt.figtext(0.5,0.95,'The sim is {}'.format(runs[i].descrip),ha='center')
        plt.subplots_adjust(wspace=0.001)

    fig.text(0.5,0.04,'Time [ps]', ha='center')
    plt.savefig('/Users/spacebob/Work/Simulations/images/{}energy.png'.format(fn), dpi=600)
    return

def msd_plots(runs, lbl, nr, c, ft):
    # Determine the filename for figure naming
    fn = runs[0].descrip[0]
    #print fn
    fig = plt.figure(2)
    for i in range(len(runs)):
        #time = (runs[i].step-runs[i].step[0])*param[i].timesteps/1000
        if 'Shell Avg' in lbl:
            ax=fig.add_subplot(nr,2,c)
            for j in range(len(runs[i].data)):
                print len(runs[i].step), len(runs[i].data[j])
                #print runs[i].descrip
                leglbl = '{} to {} $\AA$'.format(runs[i].descrip[-1]*j,runs[i].descrip[-1]*(j+1))
                ax.plot(runs[i].step,runs[i].data[j],label = leglbl)
            if 'MSD' in lbl:
                #print runs[i].descrip[0]
                ax.set_ylabel('{} ($\AA^2$)'.format(lbl))
                ax.locator_params(axis='x', tight=True, nbins=4)
            elif 'KE' in lbl:
                ax.set_ylabel('{} ($eV$)'.format(lbl))
                ax.locator_params(axis='x', tight=True, nbins=4)
                ax.legend(fontsize=10, loc=1)
                #ax.set_ylim([0,100])
            fig.text(0.5,0.04,'Time [ps]', ha='center')
        else:
            ax=fig.add_subplot(nr,1,c)
            #print runs[i].descrip
            ax.plot(runs[i].step,runs[i].data,label = runs[i].descrip[9])
            ax.set_ylabel(lbl)
            if c==2 and ft==0:
                ax.set_xlabel('Time [ps]')
            #else:
            #    plt.setp(ax.get_xticklabels(), visible=False)
            if c==1 and ft==1:
                #ax.set_xlabel('Time [ps]')
                plt.setp(ax.get_xticklabels(), visible=True)
            #else:
            #    plt.setp(ax.get_xticklabels(), visible=False)
                
            ax.locator_params(axis='x', tight=True, nbins=10)
            ax.ticklabel_format(axis='x', style = 'sci')

    plt.subplots_adjust(wspace=0.25)
    #plt.tight_layout(h_pad=0.05, w_pad=0.05)
    #plt.figtext(0.5,0.95,'The sim is {}'.format(runs[i].descrip),ha='center')
    plt.savefig('/Users/spacebob/Work/Simulations/images/{}dens.png'.format(fn),dpi=600)
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
