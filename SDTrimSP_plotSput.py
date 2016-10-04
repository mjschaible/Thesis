import numpy as np
#import scipy
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import itertools

def plot_sputExpt(log, tar, lbl=None):
    for i in range(len(log.label[0])):
        expt_tar=log.label[0][i]
        ion=log.label[1][i]
        auth=log.label[3][i]
        mk='s'
        if auth=='Roth_1979':
            mfc='k'
            lauth = 'Roth et al., 1979'
        elif 'Ken' in auth:
            lauth = 'Wehner and KenKnight, 1967'
            mfc='none'
        else:
            mfc='None'
 
        for j,nsim in enumerate(tar):
            nsim=nsim.replace('_', '->', 1)
            if expt_tar in nsim:
                ttar=nsim.replace('->', r'$\rightarrow$')
                ntar=tar[j].replace('->', '_')
                fig = plt.figure(1)
                #fig.suptitle("The simulation run is {0}".format(log[0].label[nf]))
                ax1 = fig.add_subplot(2,2,j+1)
                yerr=np.nanmax(log.totYld[i])*0.2
                ax1.scatter(log.energy,log.totYld[i],
                             label=lauth,marker=mk,facecolor=mfc,color='k',s=75)
                ax1.errorbar(log.energy,log.totYld[i],yerr=yerr,color='k')
                
                #plt.savefig('/Users/spacebob/Work/Simulations/images/{}expt.png'.format(ntar),dpi=600)
    return

def plot_srim(log,tar,lbl=None):
    for i in range(len(log.label[0])):
        expt_tar=log.label[0][i]
        ion=log.label[1][i]
        SBEO=log.label[4][i]
        SBEX=log.label[5][i]

        if 'Al' in expt_tar:
            metal='Al'
        elif 'Si' in expt_tar:
            metal='Si'
        lbl = 'SRIM, $E_{}$({})={} eV, $E_{}$(O)={} eV'.format('{SBE}',metal,SBEX,'{SBE}',SBEO)
        mk='o'
        
        for j, nsim in enumerate(tar):
            nsim=nsim.replace('_', '->', 1)
            if expt_tar in nsim:
                ttar=tar[j].replace('->', '$\rightarrow$')
                ntar=tar[j].replace('->', '_')
                fig = plt.figure(1)
                ax1 = fig.add_subplot(2,2,j+1)
                ax1.plot(log.energy,log.totYld[i], label = lbl, linestyle='--',marker=mk)
                #plt.savefig('/Users/spacebob/Work/Simulations/images/{}srim.png'.format(ntar),dpi=600)
    return

def plot_log(log,nf,ls,c,lbl=None,h=None,mc=None):
    fig = plt.figure(nf)
    ax1 = fig.add_subplot(111)
    w=0
    if lbl=='bar':
        color=iter(plt.cm.rainbow(np.linspace(0,1,len(log))))
        ind=np.arange(len(log[0].label[4]))
        width = float(1)/(len(log)+1)
        for i in log:
            l=i.label[0]
            c=next(color)
            if len(i.Flux) > len(log[0].label[4]):
                ax1.bar(ind+w*width, i.Flux[1:],width,color=c,label=l,log=1, hatch=h)
            else:
                ax1.bar(ind+w*width, i.Flux, width, color=c,label=l,log=1, hatch=h)
                x_labels=i.label[4]
            w+=1
        ax1.set_xticks(ind+width*w/2)
        ax1.set_xticklabels(x_labels)
        ax1.set_ylim(bottom=0.0001)

        plt.legend(bbox_to_anchor=(1,1), loc='upper left', fontsize=10)
        plt.subplots_adjust(right=0.85, wspace=0.25, hspace=0.35)
        fname='../../images/{}_bar.png'.format(mc)
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)
    else:    
        altx=np.linspace(100,1e4,num=50, endpoint=True)
        isbv_pos=[1,3,5]
        mk = ['s','*','^']
        e_sort = [[] for i in range(len(isbv_pos))]
        yld_sort = [[] for i in range(len(isbv_pos))]
        for i in range(len(log)):
            isbv=log[i].label[1]
            for j in range(len(isbv_pos)):
                match="isbv={}".format(isbv_pos[j])
                if match==isbv:
                    e_sort[j].append(log[i].energy)
                    yld_sort[j].append(log[i].totYld)

        for j in range(len(isbv_pos)):
            yld_sort[j] = [x for (y,x) in sorted(zip(e_sort[j],yld_sort[j]), key=lambda pair:pair[0])]
            e_sort[j].sort()

        for j in range(len(isbv_pos)):
            if e_sort[j] != []:
                if lbl != None and w==0:
                    w=1
                    ax1.plot(e_sort[j],yld_sort[j],c=c,label=lbl,ls=ls,marker=mk[j])
                else:
                    ax1.plot(e_sort[j],yld_sort[j],c=c,ls=ls,marker=mk[j])
            leg = ax1.legend()
    return

def plot_out1(sput, nf, c, mk, lbl=None):
    isbv_pos=['isbv1']#,'isbv3','isbv5']
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(2,2,nf)

    iontar=sput[0][0].label[0]
    ttar=iontar.replace('->', r'$\rightarrow$')
    ntar=iontar.replace('->', '')
    ion=sput[0][0].label[1].split('->')[0]
    ax1.text(.9,.25,"{}".format(ttar), fontsize=14, horizontalalignment='right', transform=ax1.transAxes)
    
    e_sort = [[] for i in range(len(isbv_pos))]
    yld_sort = [[] for i in range(len(isbv_pos))]
    yld_yerr = [[] for i in range(len(isbv_pos))]
    for sputy in sput:
        for i in range(len(sputy)):
            if sputy[i].label[1] != iontar:
                print 'Different ions used'
                exit

            isbv=sputy[i].label[3]

            energy = np.mean(sputy[i].energy)
            totYld = np.mean(sputy[i].totYld[len(sputy[i].totYld)-10:])
            totYld_yerr = np.std(sputy[i].totYld[len(sputy[i].totYld)-10:])
            for j in range(len(isbv_pos)):
                if isbv_pos[j]==isbv:
                    e_sort[j].append(energy)
                    yld_sort[j].append(totYld)
                    yld_yerr[j].append(totYld_yerr)
    for j in range(len(isbv_pos)):
        yld_sort[j] = [x for (y,x) in sorted(zip(e_sort[j],yld_sort[j]), key=lambda pair:pair[0])]
        yld_yerr[j] = [x for (y,x) in sorted(zip(e_sort[j],yld_yerr[j]), key=lambda pair:pair[0])]
        e_sort[j].sort()

    lblSBE=''
    for j in range(len(isbv_pos)):
        if yld_sort[j]:
            for k in range(2,len(sput[0][0].label[4])):
                elemSBE=', $E_{}$({})={:.2f} eV'.format('{SBE}', sput[0][0].label[4][k], float(sput[0][0].label[6][k]))
                lblSBE=lblSBE+elemSBE
            lbl=r'SDTrimSP{}'.format(lblSBE)
        else:
            lbl=None
        ax1.plot(e_sort[j],yld_sort[j],label=lbl,marker=mk)
        #ax1.errorbar(e_sort[j],yld_sort[j],yerr=yld_yerr[j])

    if nf==1:
        ax1.set_zorder(100)
        leg = ax1.legend(prop={'size':11}, bbox_to_anchor=(1.15, 1.05)) #handles=[ax1], loc=1, 
    #ax1.set_yscale('log')
    #ax1.set_xscale('symlog')
    ax1.set_xlim([0,10000])
    
    if ion == 'H':
        ax1.set_ylim([0,0.1])
    elif ion == 'He':
        ax1.set_ylim([0,0.65])
    if nf==1 or nf==3:
        ax1.set_ylabel('Yield (atom/ion)', fontsize=18)
    if nf==2 or nf==4:
        ax1.yaxis.set_ticklabels([])
    if nf<3:
        ax1.xaxis.set_ticklabels([])
    if nf>2:
        ax1.set_xlabel('Energy (eV)', fontsize=18)
    if nf==4:
        plt.tight_layout()

    #plt.savefig('/Users/spacebob/Work/Simulations/images/SOxYield_vsE.png'.format(ntar),dpi=600)
    return

#Function to plot the total and elemental yields as a function of flux
def plot_vflu(sput, lbl=None):
    level=0
    for n, i in enumerate(sput):
        ion=i.label[1].split('->')[0]
        energy = i.energy[0]
        Esbe_O = i.label[6][-1]
        tar=i.label[1].split('->')[1]
        ttar= r'{}$\rightarrow${}'.format(ion,tar)
        ntar='{}_{}'.format(ion,tar)
        
        if ion =='H' and energy== 1000 and Esbe_O==2.0:
            if 'Al2O3' in ntar and Esbe_O==2.0:
                fig1 = plt.figure(4)
            elif 'SiO2' in ntar and Esbe_O==2.0:
                fig1 = plt.figure(5)
            level=1
            ax = fig1.add_subplot(1,2,level)
            ax.set_ylim(3e-3,.1)
            ax.set_ylabel('Total Sputter Yield (atoms/ion)')
            ax.set_title(ttar)
            ax.set_xlim(i.fluence[1],i.fluence[-1])
            ax.set_xlabel('Fluence (x$10^{16}$)')
        elif ion == 'He' and energy== 4000 and Esbe_O==2.0:
            if 'Al2O3' in ntar and Esbe_O==2.0:
                fig1 = plt.figure(4)
            elif 'SiO2' in ntar and Esbe_O==2.0:
                fig1 = plt.figure(5)
            level=2
            ax = fig1.add_subplot(1,2,level)
            ax.set_ylim(3e-2,1)
            #ax.get_yaxis().set_visible(False)
            ax.set_title(ttar)
            ax.set_xlim(i.fluence[1],i.fluence[-1])
            ax.set_xlabel('Fluence (x$10^{16}$)')
        elif ion == 'SW' and lbl=='Met':
            level=1
            fig1=plt.figure()
            ax=fig1.add_subplot(111)
            ax.set_ylim(1e-6,2e-3)
            ax.set_ylabel('Secondary Ion Sputter Yield (ions/ion)')

        # -- Plot the yield summed over all species --
        tlbl = 'Tot Yld'
        p=0
        if ion == 'H' and energy== 1000 and Esbe_O==2.0:
            ax.loglog(i.fluence,i.totYld,c='k',lw=1,ls='-',label=tlbl)
            p=1
        elif ion == 'He' and energy== 4000 and Esbe_O==2.0:
            ax.loglog(i.fluence,i.totYld,c='k',lw=1,ls='-',label=tlbl)
            p=1
        elif ion =='SW':
            ax.loglog(i.fluence,i.totYld,c='k',lw=1,ls='-',label=tlbl)
            
        
        # -- Plot the individual species yields -- 
        c2=iter(plt.cm.rainbow(np.linspace(0,1,5)))
        for y,elem in enumerate(i.label[4]):
            if elem=='O' or elem=='Al' or elem=='Si' or elem=='Fe' or elem=='Mn':
                ce=next(c2)
                if not isinstance(i.Flux[y], (int,long)) and y>0:
                    if p==1:
                        elbl=elem
                    else:
                        elbl=None
                    
                    if ion == 'H' and energy== 1000 and Esbe_O==2.0:
                        ax.loglog(i.fluence,i.Flux[y],c=ce,lw=2,ls='--',label=elbl)
                    elif ion == 'He' and energy== 4000 and Esbe_O==2.0:
                        ax.loglog(i.fluence,i.Flux[y],c=ce,lw=2,ls='--',label=elbl)
                    elif ion =='SW':
                        ax.loglog(i.fluence,i.Flux[y],c=ce,lw=2,ls='--',label=elbl)
                        
            handles, labels = plt.gca().get_legend_handles_labels()

    if lbl=='SOx':
        if level==2:
            plt.legend(loc=1,fontsize=14)
            
    if lbl=='Met':
        i =1
        while i<len(labels):
            if labels[i] in labels[:i]:
                del(labels[i])
                del(handles[i])
            else:
                i +=1
        if level==1:
            plt.legend(handles, labels,loc=2,fontsize=12)
        if level==2:
            plt.tight_layout()
    #plt.savefig('/Users/spacebob/Work/Simulations/images/{}_vsF.png'.format(ntar),dpi=600)
    return

def plot_iavg(sput, nf, ct, shift, mk, lbl=None): 
    fig1 = plt.figure(nf)

    if nf==1:
        plotElem=['C_g','O','Na','Mg','Al','Si','S','K','Ca','Mn','Fe']
    elif nf==2:
        plotElem=['Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe']
    # --Calculate elemental averages and variances for each target--
    #elem_iyld=[[0]*len(sput[0].label[4]) for i in range(len(sput))]
    elem_iyld=[[0]*len(plotElem) for i in range(len(sput))]
    elem_iyld_var=[[0]*len(plotElem) for i in range(len(sput))]
    for n,i in enumerate(sput):
        ni=0
        a=i.label[1].split('->')[0]
        if a =='H':
            level=1
            ax = fig1.add_subplot(1,2,level)
            ax.set_ylim(5e-5,1)
            ax.set_ylabel('Sputter Yield (atoms/ion)')
        elif a == 'He':
            level=2
            ax = fig1.add_subplot(1,2,level)
            ax.set_ylim(5e-5,1)
            ax.get_yaxis().set_visible(False)
        elif a == 'SW':
            level=1
            ax=fig1.add_subplot(111)
            ax.set_yscale('log')
            if nf==1:
                ax.set_ylabel('SW Total Sputter Yield (atoms/ion)')
            elif nf==2:
                ax.set_ylabel('SW Secondary Ion Sputter Yield (ions/ion)')
            
        c2=iter(plt.cm.rainbow(np.linspace(0,1,len(i.Flux))))
        for y,elem in enumerate(i.label[4]):
            if not isinstance(i.Flux[y], (int,long)) and y>1 and elem in plotElem:
                ni=plotElem.index(elem)
                dd=len(i.Flux[y])
                navg=10
                elem_iyld[n][ni]=np.mean(i.Flux[y][dd-navg:])
                elem_iyld_var[n][ni]=np.std(i.Flux[y][dd-navg:])
                    
        tot_iyld=np.sum(elem_iyld)

    sw_iavg=np.zeros(len(elem_iyld[0]))
    myarray = np.asarray(elem_iyld)
    metclass_mean=np.mean(myarray, axis=0)
    metclass_std=np.std(myarray, axis=0)
    #for y, elem in enumerate(i.label[4]):
    #    print elem, metclass_mean[y], metclass_std[y]
    tot_iyld=np.sum(metclass_mean)
    
    pos=np.arange(shift,len(metclass_mean)+shift,1.0)
    ax.plot(pos,metclass_mean,c=ct,marker=mk,label=lbl,markersize=10,markeredgewidth=0.0,linestyle='')
    ax.errorbar(pos,metclass_mean,yerr=metclass_std,linestyle='',c=ct)
    #labels = [item.get_text() for item in ax.get_xticklabels()]
    #labels = i.label[4]
    ax.set_xticks(np.arange(0,len(metclass_mean)+1,1.0))
    ax.set_xticklabels(plotElem)
    ax.set_xlim(min(pos)-1, max(pos)+1)
    plt.legend(loc=1,fontsize=12)
    #plt.autoscale(enable=True, axis='x', tight=True)
    if nf==1:
        ax.set_ylim(1e-5,2e-2)
        plt.savefig('/Users/spacebob/Work/Simulations/images/SWtYldComp.png',dpi=600)
        plt.savefig('/Users/spacebob/Box Sync/Thesis/phd/images/SWtYldComp.png',dpi=600)
    if nf==2:
        ax.set_ylim(1e-5,1e-3)
        plt.savefig('/Users/spacebob/Work/Simulations/images/SWiYldComp.png',dpi=600)
        plt.savefig('/Users/spacebob/Box Sync/Thesis/phd/images/SWiYldComp.png',dpi=600)
    
    '''    for k in range(len(isbv_pos)):
            isbv=i.label[3]
            #print isbv
            if isbv_pos[k]==isbv and plotly==1:
                ax.loglog(i.fluence,i.totYld,lw=1,ls='-',label='Tot. yield',c='k')
                if '65901' in descrip:
                    ax2.loglog(i.fluence,i.totYld,lw=1,ls='-',label='Tot. yield',c='k')
                for j in range(1,len(i.Flux)):
                    c=next(c2)
                    #print i.label[4][j]
                    elem='{} yield'.format(i.label[4][j])
                    if tag=='SOx':
                        if j>0:
                            ax.semilogx(i.fluence,i.Flux[j],lw=2,ls='--',label=elem,c=c)
                    elif tag=='Met':
                        #print elem
                        if '65901' in descrip:
                            ax2.loglog(i.fluence,i.Flux[j],lw=2,label=elem,c=c)
                            tttar=descrip.replace('->', r'$\rightarrow$')
                            ax2.set_title(r'{:.0f}eV {}'.format(i.energy[0],tttar))
                        if 'Si' in elem or 'O' in elem:
                            ax.loglog(i.fluence,i.Flux[j],lw=2,label=elem,c=c)

                ax.set_xlabel('Fluence (x$10^{16}$)')
                ax.set_xlim(i.fluence[0],i.fluence[-1])
                ax.set_title(r'{:.0f}eV {}$\rightarrow${}'.format(i.energy[0],ion,ttar))
                ax2.set_xlabel('Fluence (x$10^{16}$)')
                ax2.set_xlim(i.fluence[0],i.fluence[-1])
                ax2.set_ylim(1e-6,1)
                if level ==1:    
                    ax.set_ylabel('Sputter Yield (atoms/ion)')
                    ax2.set_ylabel('Sputter Yield (atoms/ion)')
                elif level==2:
                    ax2.yaxis.set_visible(False)
                    ax2.legend(bbox_to_anchor=(1,1), loc='best', fontsize=12)
                    fig2.subplots_adjust(right=0.8)

        ni+=1

    if tag=='SOx':
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles, labels,loc=1,fontsize=12)#,bbox_to_anchor=(1.2,0.75)

    
    '''
    return
    
def plot_ER(ER, c, nf, lbl, targets):
    #print ER
    #print lbl
    xlabels=['Mg/Si', 'Mg/Si', 'Mg/Si', 'Al/Si', 'Al/Si', 'Ca/Si']
    x_index=[0, 0, 0, 1, 1, 2] 
    ylabels=['Al/Si', 'Ca/Si', 'Fe/Si', 'Ca/Si', 'Fe/Si', 'Fe/Si']
    y_index=[1, 2, 3, 2, 3, 3]
    fig = plt.figure(nf)
    handles=[]
    labels=[]
    figs=[]
    axs=[]
    
    for i in range(len(xlabels)):
        figs.append(plt.figure(nf+i+10))
        axs.append( figs[i].add_subplot(111) )
        ax = fig.add_subplot(2,3,i+1)
                
        for j, ratio in enumerate(ER):
            if j==0:
                ax.scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,label=lbl,s=100)
                axs[i].scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,label=lbl,s=100)
            else:
                ax.scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,s=100)
                axs[i].scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,s=100)
            ax.set_xlabel(xlabels[i])
            ax.set_ylabel(ylabels[i])
            axs[i].set_xlabel(xlabels[i])
            axs[i].set_ylabel(ylabels[i])
            axs[i].legend(bbox_to_anchor=(1.25,1.0), loc='best', fontsize=12)
            figs[i].subplots_adjust(right=0.8)
            
            if lbl=='CCs':
                ax.set_xlim(left=0)
                ax.set_ylim(bottom=0)
                axs[i].set_xlim(left=0)
                axs[i].set_ylim(bottom=0)

                #print targets
                #print ER[y_index[i]]
                #print ER[x_index[i]]
                #if lbl=='MetStds':
                #    ax1.annotate(
                #        targets[j], 
                #        xy = (ratio[x_index[i]],ratio[y_index[i]]), xytext = (-20, 20),
                #        textcoords = 'offset points', ha = 'right', va = 'bottom',
                #        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
                #bbox = dict(boxstyle = 'round,pad=0.5', fc = 'blue', alpha = 0.5),
                
    handles, labels = fig.gca().get_legend_handles_labels()
    #print handles, labels
    fig.legend(handles, labels, bbox_to_anchor=(0.1,.95),loc='center left',ncol=10,fontsize=12)
    plt.subplots_adjust(top=0.9, wspace=0.25, hspace=0.35)
    return
