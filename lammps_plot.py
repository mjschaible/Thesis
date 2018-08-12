import numpy as np
import scipy
import csv

import matplotlib.pyplot as plt
from scipy import interpolate
import itertools
import re
import os

def log_plots(runs, ft, fignum=None):
    row=2 #number of rows to show in image
    kcalToeV = 0.0433641 # convert kcal/mole to eV/atom
    # ft is the 'filetype', either rad or eq or...
    # row is the number of rows
    # runs is the data object
    # Determine the filename for figure naming
    fn = runs[0].descrip[0]
    fn = fn.replace('tip4p','')
    nMolec = runs[0].descrip[3]
    if ft==1:
        exE= runs[0].descrip[6]/nMolec
    else:
        exE=1

    col = len(runs) #num columns
    if ft==0:
        fig=plt.figure(fignum, figsize=(col+1,col/row))
    elif ft==1:
        fig=plt.figure(fignum, figsize=(3*col,3*row))
        
    pemin=0
    pemax=0
    kemin=0
    kemax=0
    ax=[]
    ax2=[]
    if row>2:
        ax3=[]

    for i,run in enumerate(runs):
        ts=runs[i].descrip[1] #timestep for current run
        time = (runs[i].step-runs[0].step[0])*ts/1000
        label1=runs[i].thermoCol[3]
        peng=runs[i].peave*kcalToeV
        if ft == 1:
            keng=runs[i].keave*kcalToeV
            label2='Sys. {} (eV)'.format(runs[i].thermoCol[5])
            units='Avg. Kin. Energy ($eV$/molec)'
        elif ft == 0:
            keng=runs[i].tempave
            label2='Avg. Sys. {} (K)'.format(runs[i].thermoCol[1])
            units='Temperature (K)'
        ncol=i
        kemin_t=min(keng)
        kemax_t=max(keng)
        pemin_t=min(peng)
        if pemin_t < pemin:
            pemin = pemin_t
        pemax_t=max(peng)
        if pemax_t < pemax:
            pemax = pemax_t
        if kemin_t > kemin:
            kemin = kemin_t
        if kemax_t > kemax:
            kemax = kemax_t

        ax.append(plt.subplot(row,col,ncol+1))    
        ax[-1].plot(time,peng,color='blue',label='Sys. {} (eV)'.format(label1))
        ax[-1].locator_params(axis='x', tight=True) 
        plt.setp(ax[-1], xticks=[min(time),max(time)])
        
        ax2.append(ax[-1].twinx())            
        ax2[-1].plot(time,keng,color='red',label=label2)
        plt.setp(ax2[-1].get_xticklabels(), visible=False)
        plt.setp(ax2[-1].get_yticklabels(), visible=False)

        if i>0:
            plt.setp(ax[-1].get_yticklabels(), visible=False)
            for ac in ax:
                ac.get_shared_y_axes().join(ac,ax[-1])
            for ac in ax2:
                ac.get_shared_y_axes().join(ac,ax2[-1])
            if ft==1:
                ax[-1].set_ylim([pemin-1/nMolec,pemin+exE+3/nMolec])
                ax2[-1].set_ylim([kemax-exE-3/nMolec,kemax+1/nMolec])
            else:
                ax[-1].set_ylim([pemin,pemax])
                ax2[-1].set_ylim([kemin,kemax])
            
        if ncol==col-1:
            #print 'PE {}'.format(pemin+exE+3-pemin+1)
            #print 'KE {}'.format(kemax+1-kemax+exE+3)
            ax[0].set_ylabel('Avg. Pot. Energy ($eV$/molec)')
            ax2[-1].set_ylabel(units)
            h1,l1=ax[-1].get_legend_handles_labels()
            h2,l2=ax2[-1].get_legend_handles_labels()
            l = ax[-1].legend(h1+h2, l1+l2, loc=4, fontsize=10)
            plt.setp(ax2[-1].get_yticklabels(), visible=True)

        try:
            ax3
            if i==ft:
                ax3.append(plt.subplot(row,col,ncol+1+col))
                ax3[-1].set_ylabel('{} [eV]'.format(runs[i].thermoCol[7]))
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
            print
            #print 'Only one axis defned'

    plt.subplots_adjust(wspace=0.001)
    plt.show()
    #fig.text(0.5,0.04,'Time [ps]', ha='center')
    #plt.savefig('/Users/spacebob/Work/Simulations/images/{}energy.png'.format(fn), dpi=600)
    #plt.savefig('/home/mikey/Work/Sinmulations/images/{}energy.png'.format(fn), dpi=600)
    return

def msd_plots(runs, lbl, ft, fignum=None):
    # Determine the filename for figure naming
    fn = runs.descrip[0]
    fn = fn.replace('tip4p','')
    fig = plt.figure(fignum)

    numshells=len(runs.descrip[-1])
    shell_t=runs.descrip[-2]

    #time = (runs[i].step-runs[i].step[0])*param[i].timesteps/1000
    if 'Shell avg.' in lbl:
        ax=fig.add_subplot(111)
        if 'MSD' in lbl:
            ax.set_ylabel('{} [$\AA^2$]'.format(lbl))
            ax.locator_params(axis='x', tight=True, nbins=4)
            tp='disp'
        elif 'KE' in lbl:
            ax.set_ylabel('{} [$eV$]'.format(lbl))
            ax.locator_params(axis='x', tight=True, nbins=4)
            tp='eng'
            #ax.set_ylim([0,100])

        for j in range(len(runs.data)):
            ns = len(runs.descrip[-1][j])
            leglbl = '{} to {} $\AA$ (N={})'.format(shell_t*j,shell_t*(j+1),ns)
            ax.semilogy(runs.step,runs.data[j],label = leglbl)

        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        ax.set_xlabel('Time [ps]')
        #plt.savefig('/Users/spacebob/Work/Simulations/images/{}{}.png'.format(fn,tp),dpi=600)
        #plt.savefig('/home/mikey/Work/Sinmulations/images/{}{}.png'.format(fn,tp),dpi=600)
    else:
        nsp=2
        ax=fig.add_subplot(2,1,nsp)
        ax.plot(runs.step,runs.data,label = runs.descrip[3])
        ax.set_ylabel(r'{} [$\AA^2$]'.format(lbl))
        ax.locator_params(axis='x', tight=True, nbins=6)
        if nsp==2 and ft==0:
            ax.set_xlabel('Time [ps]')
        #else:
        #    plt.setp(ax.get_xticklabels(), visible=False)
        if nsp==1 and ft==1:
            #ax.set_xlabel('Time [ps]')
            plt.setp(ax.get_xticklabels(), visible=True)
        #else:
        #    plt.setp(ax.get_xticklabels(), visible=False)

        ax.ticklabel_format(axis='x', style = 'sci')

    return

def rdf_plots(runs):
    fn = runs[0].descrip[0]
    fn = fn[:-3]
    fn = fn.replace('tip4p','')

    for i in range(len(runs)):
        fig, axs = plt.subplots(2,1, sharex='col')
        axs[0].plot(runs[i].pos,runs[i].RDF)
        axs[0].set_ylabel('Radial dist. fct. $g_{OO}$(r)')
        axs[1].plot(runs[i].pos,runs[i].coordN)
        axs[1].set_ylabel('Coordination number')
        axs[1].set_xlabel(r'Distance [$\AA$]')
        #plt.savefig('/Users/spacebob/Work/Simulations/images/{}rdf.png'.format(fn),dpi=600)
        #plt.savefig('/home/mikey/Work/Sinmulations/images/{}rdf.png'.format(fn),dpi=600)
    return fn

def com_plots(msd):
    plt.figure()
    fn=msd[0].descrip[0]
    fn = fn[:-3]
    fn = fn.replace('tip4p','')
    shell_avg=[]
    
    for i, run in enumerate(msd):
        shell_thickness=run.descrip[7]
        num_shells=len(run.descrip[8])
        dist_arr=np.arange(0,num_shells*shell_thickness,shell_thickness)
        dist_arr=[x+shell_thickness/2 for x in dist_arr]
        shell_avg.append(run.avg)

    shells_avg=np.mean(shell_avg, axis=0)
    shell_stdev=np.std(shell_avg, axis=0)
    
    plt.plot(dist_arr,shells_avg,marker='o',markersize=10,linestyle='')
    plt.errorbar(dist_arr,shells_avg,yerr=shell_stdev,linestyle='')
    plt.yscale('log')
    plt.xlabel(r'Distance from PKA[$\AA$]')
    plt.ylabel(r'Mean MSD [$\AA^2$]')
    #plt.savefig('/Users/spacebob/Work/Simulations/images/{}MSDvDist.png'.format(fn),dpi=600)

    a = [dist_arr,shells_avg,shell_stdev]
    print 'This is being written to text'
    print dist_arr
    print 'then'
    print shells_avg
    print 'and combined'
    print a
    
    np.savetxt('{}_MSD_runAvgs.res'.format(fn), a, delimiter=",")

    with open('{}_MSD_runAvgs.csv'.format(fn), 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        header = ["Shell Distance, Average, StDev"]
        writer.writerow(header)
        for n in range(len(dist_arr)):
            line = [dist_arr[n],shells_avg[n],shell_stdev[n]]
            writer.writerow(line)
    return
