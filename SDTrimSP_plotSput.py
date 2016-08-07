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
        if ion=='H':
            mk='s'
        elif ion=='He':
            mk='o'
        elif ion=='D':
            mk='^'
        if auth=='Roth_1979':
            mfc='k'
            lauth = 'Roth et al., 1979'
        elif 'Ken' in auth:
            lauth = 'Wehner and KenKnight, 1967'
            mfc='none'
        else:
            mfc='None'
 
        for j in range(len(tar)):
            if expt_tar in tar[j]:
                ttar=tar[j].replace('->', r'$\rightarrow$')
                ntar=tar[j].replace('->', '_')
                fig = plt.figure(j)
                #fig.suptitle("The simulation run is {0}".format(log[0].label[nf]))
                ax1 = fig.add_subplot(111)
                #print expt_tar
                yerr=np.nanmax(log.totYld[i])*0.2
                #print yerr
                ax1.scatter(log.energy,log.totYld[i],
                             label=lauth,marker=mk,facecolor=mfc,color='k',s=75)
                ax1.errorbar(log.energy,log.totYld[i],yerr=yerr,color='k')
                ax1.set_ylabel('Yield (atom/ion)', fontsize=18)
                ax1.set_xlabel('Energy (eV)', fontsize=18)
                fig.suptitle("{}".format(ttar), fontsize=20, y=0.96)
                leg1=ax1.legend(prop={'size':11})
                #plt.savefig('/Users/spacebob/Work/Simulations/images/{}expt.png'.format(ntar),dpi=600)
    return

def plot_srim(log,tar,lbl=None):
    for i in range(len(log.label[0])):
        expt_tar=log.label[0][i]
        ion=log.label[1][i]
        SBEO=log.label[4][i]
        SBEX=log.label[5][i]
        if 'H->'in ion:
            mk='s'
        elif 'He->' in ion:
            mk='o'
        elif 'D->' in ion:
            mk='^'
            
        if 'Al' in expt_tar:
            metal='Al'
        elif 'Si' in expt_tar:
            metal='Si'
            
        lbl = 'SRIM, $E_{}$({})={} eV, $E_{}$(O)={} eV'.format('{SBE}',metal,SBEX,'{SBE}',SBEO)
        if ion=='H':
            mk='s'
        elif ion=='He':
            mk='o'
        elif ion=='D':
            mk='^'
        for j in range(len(tar)):
            if expt_tar in tar[j]:
                ttar=tar[j].replace('->', '$\rightarrow$')
                ntar=tar[j].replace('->', '_')
                fig = plt.figure(j)
                ax1 = fig.add_subplot(111)
                ax1.plot(log.energy,log.totYld[i], label = lbl, linestyle='--',marker=mk)
                leg=ax1.legend(prop={'size':11})
                #plt.savefig('/Users/spacebob/Work/Simulations/images/{}srim.png'.format(ntar),dpi=600)
    return

def plot_log(log,nf,ls,c,lbl=None,h=None,mc=None):
    fig = plt.figure(nf)
    ax1 = fig.add_subplot(111)
    w=0
    if lbl=='bar':
        color=iter(plt.cm.rainbow(np.linspace(0,1,len(log))))
        ind=np.arange(len(log[0].label[4]))
        #print len(ind), len(log[0].Flux)
        width = float(1)/(len(log)+1)
        for i in log:
            l=i.label[0]
            c=next(color)
            #print len(i.Flux), i.label[4]
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
                    #print log[i].label[0], log[i].label[1]
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

def plot_out1(sput, nf, ls, c, lbl=None):
    isbv_pos=['isbv1']#,'isbv3','isbv5']
    #mk = ['s','*','^']

    fig1 = plt.figure(nf)
    ax1 = fig1.add_subplot(111)

    ion=sput[0][0].label[1]
    if 'H->'in ion:
        mk='s'
    elif 'He->' in ion:
        mk='o'
    elif 'D->' in ion:
        mk='^'

    e_sort = [[] for i in range(len(isbv_pos))]
    yld_sort = [[] for i in range(len(isbv_pos))]
    yld_yerr = [[] for i in range(len(isbv_pos))]
    for sputy in sput:
        for i in range(len(sputy)):
            if sputy[i].label[1] != ion:
                print 'Different ions used'
                exit

            isbv=sputy[i].label[3]
            #print isbv
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
        
    tar=sput[0][0].label[0]
    ttar=tar.replace('->', '$\rightarrow$')
    ntar=tar.replace('->', '')

    lblSBE=''
    for j in range(len(isbv_pos)):
        if yld_sort[j]:
            #print sput[0].label[4]
            #print sput[0].label[5]
            for k in range(1,len(sput[0][0].label[4])):
                elemSBE=', $E_{}$({})={:.2f} eV'.format('{SBE}', sput[0][0].label[4][k], float(sput[0][0].label[6][k]))
                lblSBE=lblSBE+elemSBE
            lbl=r'SDTrimSP{}'.format(lblSBE)
        else:
            lbl=None
        ax1.plot(e_sort[j],yld_sort[j],label=lbl,ls=ls,marker=mk)
        #ax1.errorbar(e_sort[j],yld_sort[j],yerr=yld_yerr[j])

    leg = ax1.legend(prop={'size':11}) #handles=[ax1], loc=1, 
    #ax1.set_yscale('log')
    #ax1.set_xscale('symlog')
    ax1.set_xlim([0,10000])
    ax1.set_ylim(bottom=0)
    #ax1.set_yticks([0.001, 0.01, 0.1, 1])
    plt.savefig('/Users/spacebob/Work/Simulations/images/{}_vsE.png'.format(ntar),dpi=600)
    return

def plot_out2(sput, nf, c, lbl=None): 
    isbv_pos=['isbv1']#,'isbv3','isbv5']
    #mk = ['s','*','^']
    fig1 = plt.figure(nf)
    
    for i in sput:
        for j, a in enumerate(i.label[4]):
            #print a
            if a == 'Si':
                nSi=j
            elif a == 'He' and j==0:
                level=1
            elif a =='H' and j==0:
                level=2
        c2=iter(plt.cm.rainbow(np.linspace(0,1,len(i.Flux))))
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        lbl = i.label[0]
        ax = fig1.add_subplot(2,1,level)            
        ax.plot(i.fluence,i.totYld,label=lbl,c=c,lw=1,ls='-')
        #if float(i.Flux[0][0]) == 1.0:
        #    pass
        #elif float(i.Flux[0][0])> 1.0:
        #    i.Flux[0]=[x-1.0 for x in i.Flux[0]]
        #    ax.plot(i.fluence,i.Flux[0],c=c,lw=2,ls=':')
        #else:
        #    ax.plot(i.fluence,i.Flux[0],c=c,lw=2,ls=':')
        ax.plot(i.fluence,i.Flux[nSi],c=c,lw=2,ls='--')
        ax.plot(i.fluence,i.Flux[-1],c=c,lw=2,ls=':')
        ax.set_xlim(i.fluence[1],i.fluence[-1])
        if level==1:
            ax.get_xaxis().set_visible(False)
        if level==2:
            ax.set_xlabel('Fluence (x$10^{16}$)')
        ax.set_ylabel('Sputter Yield (atoms/ion)')
        #ax.set_title(i.label[0])
    handles, labels = plt.gca().get_legend_handles_labels()
    i =1
    while i<len(labels):
        if labels[i] in labels[:i]:
            del(labels[i])
            del(handles[i])
        else:
            i +=1
    if level==1:
        plt.legend(handles, labels,bbox_to_anchor=(1.15,0.45),loc='center',fontsize=12)
    plt.subplots_adjust(right=0.8, wspace=0.25, hspace=0.35)
    return

def plot_out3(sput, mets, tag=None): 
    isbv_pos=['isbv1']#,'isbv3','isbv5']
    #mk = ['s','*','^']
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax = fig1.add_subplot(2,1,1)
    ax2=fig2.add_subplot(111)
    ni=0
    for i in sput:
        c2=iter(plt.cm.rainbow(np.linspace(0,1,len(i.Flux))))
        tar=i.label[0]
        descrip=i.label[1]
        ion=i.label[1].split('->')[0]
        ttar=tar.replace('->', r'$\rightarrow$')
        ntar=tar.replace('->', '')
        #print tar, ion
        for j, a in enumerate(i.label[4]):
            #print a
            if a == 'Si':
                nSi=j
            elif a == 'He' and j==0:
                level=1
            elif a =='H' or 'D' and j==0:
                level=2
            plotly=0
            if i.energy[0] == 1000 and level==2:
                plotly=1
            elif i.energy[0] == 4000 and level==1:
                plotly=1

        for j in range(len(isbv_pos)):
            isbv=i.label[3]
            #print isbv
            if isbv_pos[j]==isbv and plotly==1:
                ax.semilogx(i.fluence,i.totYld,lw=1,ls='-',label='Tot. yield',c='k')
                for j in range(len(i.Flux)):
                    c=next(c2)
                    #print i.label[4][j]
                    elem='{} yield'.format(i.label[4][j])
                    if float(i.Flux[0][0])>i.totYld[0]:
                        pass
                    else:
                        ax.semilogx(i.fluence,i.Flux[0],lw=2,ls=':',label=elem,c=c)
                    if tag=='SOx':
                        if j>0:
                            ax.semilogx(i.fluence,i.Flux[j],lw=2,ls='--',label=elem,c=c)
                    elif tag=='Met':
                        #print elem
                        if '65901' in descrip:
                            ax2.loglog(i.fluence,i.Flux[j],lw=2,label=elem,c=c)
                            ax2.legend(bbox_to_anchor=(1,1), loc='best', fontsize=12)
                            fig2.subplots_adjust(right=0.8)
                        if 'Si' in elem or 'O' in elem:
                            ax.loglog(i.fluence,i.Flux[j],lw=2,label=elem,c=c)
                        if ni==0:
                            ax.legend(loc=3)

                ax.set_xlabel('Fluence (x$10^{16}$)')
                ax.set_ylabel('Sputter Yield (atoms/ion)')
                ax.set_xlim(i.fluence[0],i.fluence[-1])
                ax.set_title(r'{:.0f}eV {}$\rightarrow${}'.format(i.energy[0],ion,ttar))
                ax2.set_xlabel('Fluence (x$10^{16}$)')
                ax2.set_ylabel('Sputter Yield (atoms/ion)')
                ax2.set_xlim(i.fluence[0],i.fluence[-1])
                ax2.set_title(r'{:.0f}eV {}$\rightarrow${}'.format(i.energy[0],ion,ttar))
        ni+=1
    
    if tag=='SOx':
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles, labels,loc=1,fontsize=12)#,bbox_to_anchor=(1.2,0.75)

    #plt.savefig('/Users/spacebob/Work/Simulations/images/{}{}_vsF.png'.format(ntar,ion),dpi=600)
    return

def plot_ER(ER, c, nf, lbl, targets):
    #print ER
    print lbl
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
