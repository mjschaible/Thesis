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
            yerr=np.nanmax(log.totYld[i])*0.2
        elif 'Ken' in auth:
            lauth = 'Wehner and KenKnight, 1967'
            mfc='none'
            yerr=np.nanmax(log.totYld[i])*0.5
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
                ax1.scatter(log.energy,log.totYld[i],
                             label=lauth,marker=mk,facecolor=mfc,color='k',s=80)
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
            c='m'
        if SBEO =='2.00':
            mfc=None
            c='c'
        else:
            mfc='None'
        lbl = 'SRIM, $E_{}$({})={} eV, $E_{}$(O)={} eV'.format('{SBE}',metal,SBEX,'{SBE}',SBEO)
        for j, nsim in enumerate(tar):
            nsim=nsim.replace('_', '->', 1)
            if expt_tar in nsim:
                ttar=tar[j].replace('->', '$\rightarrow$')
                ntar=tar[j].replace('->', '_')
                fig = plt.figure(1)
                ax1 = fig.add_subplot(2,2,j+1)
                ax1.plot(log.energy,log.totYld[i],c=c,label = lbl,linestyle='--',marker='o',markerfacecolor=mfc,markeredgewidth=2,markeredgecolor=c)
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

def plot_out1(sput, nf, c, mk,mfc, lbl=None):
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
        ax1.plot(e_sort[j],yld_sort[j],label=lbl,color=c,marker=mk,markerfacecolor=mfc,markeredgewidth=2,markeredgecolor=c)
        #ax1.errorbar(e_sort[j],yld_sort[j],yerr=yld_yerr[j])

    if nf==1:
        ax1.set_zorder(100)
        leg = ax1.legend(prop={'size':11},bbox_to_anchor=(0.35, 1)) #handles=[ax1], loc=1, 
    #ax1.set_yscale('log')
    #ax1.set_xscale('symlog')
    ax1.set_xlim([0,9900])
    
    if ion == 'H':
        ax1.set_ylim([0,0.1])
    elif ion == 'He':
        ax1.set_ylim([0,0.65])
    if nf==1 or nf==3:
        ax1.set_ylabel('Yield (atom/ion)')
    if nf==2 or nf==4:
        ax1.yaxis.set_ticklabels([])
    if nf<3:
        ax1.xaxis.set_ticklabels([])
    if nf>2:
        ax1.set_xlabel('Energy (eV)')
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
            ax.set_ylim(0,.04)
            ax.set_ylabel('Sputter Yield (atoms/ion)')
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
            ax.set_ylim(0,0.4)
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
        tlbl = 'Tot Yield, $Y^{tot}$'
        p=0
        if ion == 'H' and energy== 1000 and Esbe_O==2.0:
            ax.semilogx(i.fluence,i.totYld,c='k',ls='-',label=tlbl)
            p=1
        elif ion == 'He' and energy== 4000 and Esbe_O==2.0:
            ax.semilogx(i.fluence,i.totYld,c='k',ls='-',label=tlbl)
            p=1
        elif ion =='SW':
            ax.loglog(i.fluence,i.totYld,c='k',ls='-',label=tlbl)
            
        
        # -- Plot the individual species yields -- 
        c2=iter(plt.cm.viridis(np.linspace(0.2,0.8,2)))
        for y,elem in enumerate(i.label[4]):
            if elem=='O' or elem=='Al' or elem=='Si' or elem=='Fe' or elem=='Mn':
                ce=next(c2)
                if not isinstance(i.Flux[y], (int,long)) and y>0:
                    if p==1:
                        elbl=r'{} Yield, $Y_{{{}}}$'.format(elem, elem)
                    else:
                        elbl=None
                    
                    if ion == 'H' and energy== 1000 and Esbe_O==2.0:
                        ax.semilogx(i.fluence,i.Flux[y],c=ce,ls='--',label=elbl)
                    elif ion == 'He' and energy== 4000 and Esbe_O==2.0:
                        ax.semilogx(i.fluence,i.Flux[y],c=ce,ls='--',label=elbl)
                    elif ion =='SW':
                        ax.loglog(i.fluence,i.Flux[y],c=ce,ls='--',label=elbl)
                        
            handles, labels = plt.gca().get_legend_handles_labels()

    if lbl=='SOx':
        if level==2:
            handles = [handles[0], handles[2], handles[1]]
            labels = [labels[0], labels[2], labels[1]]
            plt.legend(handles, labels, loc=1)
            
    if lbl=='Met':
        i =1
        while i<len(labels):
            if labels[i] in labels[:i]:
                del(labels[i])
                del(handles[i])
            else:
                i +=1
        if level==1:
            plt.legend(handles, labels,loc=2)
        if level==2:
            plt.tight_layout()
    #plt.savefig('/Users/spacebob/Work/Simulations/images/{}_vsF.png'.format(ntar),dpi=600)
    return

def plot_iavg(sput, nf, ct,shift, mk, lbl=None):
    print lbl
    fig1 = plt.figure(nf)
    allElem=['C_g','O','Na','Mg','Al','Si','S','K','Ca','Mn','Fe']
    if nf==1:
        plotElem=allElem
        t = 'Neutral'
        nSi=5
    elif nf==2:
        plotElem=['Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe']
        t= 'Ion'
        nSi=3
    # --Calculate elemental averages and variances for each target--
    #elem_iyld=[[0]*len(sput[0].label[4]) for i in range(len(sput))]
    elem_iyld=[[0]*len(plotElem) for i in range(len(sput))]
    elem_iyld_var=[[0]*len(plotElem) for i in range(len(sput))]
    for n,i in enumerate(sput):
        ni=0
        eng=i.label[0].split()[0] #variable for the incident energy
        ion=i.label[0].split()[1] # variable for the incident ion
        if ion =='H':
            level=1
            ax = fig1.add_subplot(1,2,level)
            ax.set_ylim(5e-5,1)
            ax.set_ylabel('Sputter Yield (atoms/ion)')
        elif ion == 'He':
            level=2
            ax = fig1.add_subplot(1,2,level)
            ax.set_ylim(5e-5,1)
            ax.get_yaxis().set_visible(False)
        else:
            level=1
            ax=fig1.add_subplot(1,2,1)
            ax2=fig1.add_subplot(1,2,2)
            ax.set_yscale('log')
            ax2.set_yscale('log')

            if nf==1:
                ax.set_ylabel('Total Sputter Yield (atoms/ion)')
                ax2.set_ylabel('Elemental Yield / Si yield')
            elif nf==2:
                ax.set_ylabel('Secondary Ion Sputter Yield (ions/ion)')
                ax2.set_ylabel('Secondary Ion Yield / Si Ion Yield')

        c2=iter(plt.cm.rainbow(np.linspace(0,1,len(i.Flux))))
        for y,elem in enumerate(i.label[3]):
            if not isinstance(i.Flux[y], (int,long)) and elem in plotElem:
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
    elem_vsSi_mean=metclass_mean/metclass_mean[nSi]
    #elem_vsSi_std=elem_vsSi_mean*((metclass_std/metclass_mean)**2+(metclass_std[5]/metclass_mean[5])**2)**0.5
    
    print 'Class Averages, {}'.format(t)
    for y, elem in enumerate(plotElem):
        print '{}, {}$\pm${}: vs. Si = {}'.format(elem, metclass_mean[y], metclass_std[y], elem_vsSi_mean[y])

    tot_iyld=np.sum(metclass_mean)
    print 'Total {} yield = {:.3}'.format(t, tot_iyld)
    
    pos=np.arange(shift,len(metclass_mean)+shift,1.0)
    for y,elem in enumerate(i.label[3]):
        if elem in plotElem:
            if elem != 'Si':
                llbl=None
            else:
                llbl=lbl

            ax.plot(pos+0.4,metclass_mean,c=ct,marker=mk,label=llbl,markersize=10,markeredgewidth=0.0,linestyle='')
            ax.errorbar(pos+0.4,metclass_mean,yerr=metclass_std,linestyle='',c=ct)
            ax2.plot(pos+0.4,elem_vsSi_mean,c=ct,marker=mk,label=llbl,markersize=10,markeredgewidth=0.0,linestyle='')

    ax.set_xlim(min(pos)-0.5, max(pos)+0.2)
    ax.set_xticks(pos+0.5,minor=True)
    ax.xaxis.set(ticks=pos,ticklabels=plotElem)
    ax2.set_xlim(min(pos)-0.5, max(pos)+0.2)
    ax2.set_xticks(pos+0.5,minor=True)
    ax2.xaxis.set(ticks=pos,ticklabels=plotElem)
    # Turn on the grid for the minor ticks
    ax.xaxis.grid(True, which='minor')
    ax2.xaxis.grid(True, which='minor')
    gridlines = ax.get_ygridlines()
    gridlines2 = ax2.get_ygridlines()
    for line in gridlines:
        line.set_linestyle('--')
    for line in gridlines2:
        line.set_linestyle('--')
    
    ax.legend(loc=1,fontsize=10)
    ax2.legend(loc=1,fontsize=10)
    #plt.autoscale(enable=True, axis='x', tight=True)
    if nf==1:
        ax.set_ylim(1e-5,5e-2)
        #plt.savefig('/Users/spacebob/Work/Simulations/images/SWtYldComp.png',dpi=600)
        #plt.savefig('/Users/spacebob/Box Sync/Thesis/phd/images/SWtYldComp.png',dpi=600)
    if nf==2:
        ax.set_ylim(3e-6,5e-3)
        ax2.set_ylim(0.01,15)
        #plt.savefig('/Users/spacebob/Work/Simulations/images/SWiYldComp.png',dpi=600)
        #plt.savefig('/Users/spacebob/Box Sync/Thesis/phd/images/SWiYldComp.png',dpi=600)
        
    return metclass_mean, metclass_std
    
def plot_ER(ER,c,lbl,targets,nf,mk):
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
                ax.scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,label=lbl,s=100,marker=mk)
                axs[i].scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,label=lbl,s=100,marker=mk)
            else:
                ax.scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,s=100,marker=mk)
                axs[i].scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,s=100,marker=mk)
            ax.set_xlabel(xlabels[i])
            ax.set_ylabel(ylabels[i])
            axs[i].set_xlabel(xlabels[i])
            axs[i].set_ylabel(ylabels[i])
            axs[i].legend(loc='best', fontsize=12) #bbox_to_anchor=(1.25,1.0), 
#            figs[i].subplots_adjust(right=0.8)
            
        if lbl=='CCs':
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)
            axs[i].set_xlim(left=0)
            axs[i].set_ylim(bottom=0)
                
    handles, labels = fig.gca().get_legend_handles_labels()
    #print handles, labels
    fig.legend(handles, labels, bbox_to_anchor=(0.1,.95),loc='center left',ncol=10,fontsize=12)
    plt.subplots_adjust(top=0.9, wspace=0.25, hspace=0.35)
    return

def plot_flux(mean_yld,tld_std,met,nf):
    color=iter(plt.cm.viridis(np.linspace(0,0.9,6)))
    marker=iter(['o', 's', 'D', '^', 'v', '<', '>'])
    if met=='Lunar':
        ttl=met
        col = 1
    elif met=='Mars':
        ttl=met
        col=2
    elif met=='CCs':
        ttl = 'Carboneceous Chondrites'
        col=3
    else:
        col=0
    if col>0:
        fig = plt.figure(nf,figsize=(12.0, 7))
        #fig2=plt.figure(nf+1)
        ax = fig.add_subplot(1,3,col)
        ax.set_title(ttl)
        if col==1:
            ax.set_ylabel('Sputtered Ion Flux (ions $cm^{-2}s^{-1}$)')
        if col==2:
            ax.set_xlabel('Distance above ~20 km diameter object (km)')

        elem=['Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe']

        phiSW=1e8 # solar wind flux, ion/cm^2/s
        P=1./3 # porosity reduction factor
        theta=3.78 # accounting for cosine distribution of incident ion angle
        theta_p=0.68404 # accounting for flux ejected within 20deg from surface normal
        r_sb = 10 # radius of small body, km
        dist = np.arange(5.,100.,5)
        dist_fac=(1-(1-(r_sb/(r_sb+dist))**2)**0.5)

        for i, yld in enumerate(mean_yld):
            c=next(color)
            mk=next(marker)
            surf_flux=phiSW*yld*P*theta*theta_p
            dist_flux=surf_flux*dist_fac
            ax.semilogy(dist,dist_flux,label=elem[i],color=c,marker=mk)
            print 'The {} flux at {} km is {} and at {} km is {}'.format(elem[i], dist[0], dist_flux[0], dist[7], dist_flux[7])
        ax.set_ylim(10,1e4)
        if col ==2:
            ax.legend(loc=1,fontsize=12)
        if col>1:
            ax.yaxis.set_ticklabels([])

        if col==3:
            plt.savefig('/home/mikey/Work/Simulations/images/FluxVsDist.png',dpi=600, bbox_inches='tight')
    return
