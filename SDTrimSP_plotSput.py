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
        else:
            mfc='None'
 
        for j in range(len(tar)):
            if expt_tar in tar[j]:
                fig = plt.figure(j)
                #fig.suptitle("The simulation run is {0}".format(log[0].label[nf]))
                ax1 = fig.add_subplot(111)
                ax1.scatter(log.energy,log.totYld[i],
                            label=auth,marker=mk,facecolor=mfc,color='k')
                ax1.set_ylabel('Yield (atom/ion)')
                ax1.set_xlabel('Energy (eV)')
                fig.suptitle("{}".format(expt_tar))
                leg1=ax1.legend(prop={'size':10})
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
        if ion=='H':
            mk='s'
        elif ion=='He':
            mk='o'
        elif ion=='D':
            mk='^'
        for j in range(len(tar)):
            if expt_tar in tar[j]:
                fig = plt.figure(j)
                ax1 = fig.add_subplot(111)
                ax1.plot(log.energy,log.totYld[i], label = lbl, linestyle='-.')
                leg=ax1.legend(prop={'size':10})
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
    e_sort = [[] for i in range(len(isbv_pos))]
    yld_sort = [[] for i in range(len(isbv_pos))]
    for i in range(len(sput)):
        isbv=sput[i].label[3]
        energy = np.mean(sput[i].energy)
        totYld = np.mean(sput[i].totYld[len(sput[i].totYld)-10:])
        for j in range(len(isbv_pos)):
            if isbv_pos[j]==isbv:
                e_sort[j].append(energy)
                yld_sort[j].append(totYld)
        ion=sput[i].label[1]
        if 'H->'in ion:
            mk='s'
        elif 'He->' in ion:
            mk='o'
        elif 'D->' in ion:
            mk='^'
            
    lblSBE=''
    for j in range(len(isbv_pos)):
        if yld_sort[j]:
            #print sput[0].label[4]
            #print sput[0].label[5]
            for k in range(1,len(sput[0].label[4])):
                elemSBE=', $E_{}$({})={} eV'.format('{SBE}', sput[0].label[4][k], sput[0].label[5][k])
                lblSBE=lblSBE+elemSBE
            lbl=r'SDTrimSP, {} {}'.format(isbv_pos[j],lblSBE)
        else:
            lbl=None
        ax1.plot(e_sort[j],yld_sort[j],label=lbl,ls=ls,c=c,marker=mk)

    leg = ax1.legend(prop={'size':10}) #handles=[ax1], loc=1, 
    #ax1.set_yscale('log')
    #ax1.set_xscale('symlog')
    ax1.set_xlim([0,10000])
    ax1.set_ylim(bottom=0)
    #ax1.set_yticks([0.001, 0.01, 0.1, 1])
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
    plt.legend(handles, labels,bbox_to_anchor=(1.2,0.75),loc='center',fontsize=10)
    plt.subplots_adjust(right=0.8, wspace=0.25, hspace=0.35)
    return

def plot_ER(ER, c, nf, lbl, targets):
    #print ER
    print lbl
    xlabels=['Mg/Si', 'Mg/Si', 'Mg/Si', 'Al/Si', 'Al/Si', 'Ca/Si']
    x_index=[0, 0, 0, 1, 1, 2] 
    ylabels=['Al/Si', 'Ca/Si', 'Fe/Si', 'Ca/Si', 'Fe/Si', 'Fe/Si']
    y_index=[1, 2, 3, 2, 3, 3]
    fig1 = plt.figure(nf)
    handles=[]
    labels=[]
    for i in range(len(xlabels)):
        ax1 = fig1.add_subplot(2,3,i+1)
        for j, ratio in enumerate(ER):
            for k in range(len(ratio)):
                if k==0 and j==0:
                    ax1.scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,label=lbl,s=30)
                else:
                    ax1.scatter(ratio[x_index[i]],ratio[y_index[i]],color=c,s=30)
                ax1.set_xlabel(xlabels[i])
                ax1.set_ylabel(ylabels[i])
            
                if lbl=='CCs':
                    ax1.set_xlim(left=0)
                    ax1.set_ylim(bottom=0)

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
    #print handles
    plt.legend(bbox_to_anchor=(-1,2.5), ncol=10, loc='center', fontsize=10)
    plt.subplots_adjust(top=0.9, wspace=0.25, hspace=0.35)
    return


