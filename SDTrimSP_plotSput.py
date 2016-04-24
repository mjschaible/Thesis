import numpy as np
import scipy
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy import interpolate
import itertools

def plot_sputExpt(log, nf, ls, lbl=None):
    fig = plt.figure(nf)

    fig.suptitle("The simulation run is {0}".format(log[0].label[nf]))
    ax1 = fig.add_subplot(111)
    ax1.scatter(log[0].energy,log[0].totYld[nf], label=lbl, marker=ls)
    ax1.set_ylabel('Yield (atom/ion)')
    ax1.set_xlabel('Energy (eV)')
    leg=ax1.legend()
    return

def plot_srim(log,nf,ls,c,lbl=None):
    fig = plt.figure(nf)

    ax1 = fig.add_subplot(111)
    ax1.plot(log.energy,log.totYld[nf], label = lbl, linestyle=ls)
    leg=ax1.legend()
    return

def plot_log(log,nf,ls,c,lbl=None):
    lll=0
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

    fig = plt.figure(nf)
#    fig.suptitle("The simulation run is {0}".format(descrip))
    ax1 = fig.add_subplot(111)
    for j in range(len(isbv_pos)):
#        spl=interpolate.splrep(e_sort[j],yld_sort[j],s=0)
#        alty=interpolate.splev(altx, spl, der=0)
        if e_sort[j] != []:
            if lbl != None and lll==0:
                lll=1
                ax1.plot(e_sort[j],yld_sort[j],c=c,label=lbl,ls=ls,marker=mk[j])
            else:
                ax1.plot(e_sort[j],yld_sort[j],c=c,ls=ls,marker=mk[j])
        leg = ax1.legend()
    ax1.set_xlim([0,10000])
    ax1.set_ylim(bottom=0)
    return

def plot_flu(sput, nf, ls, labl=None): 
    fig = plt.figure(nf)
    ax1 = fig.add_subplot(111)
    asdf=np.array(range(len(sput)))
    c1=iter(plt.cm.inferno(np.linspace(0,1,5)))
    c2=iter(plt.cm.rainbow(np.linspace(0,1,8)))
    #ls = itertools.cycle(('-','--','-.',':'))
    isbv_pos=[1]
    lbl=[]
    e_sort = []
    yld_sort = []
    for i in range(len(sput)):
        isbv=sput[i].label[1]
        match="isbv={}".format(isbv_pos[0])
        if match==isbv:
            #print log[i].label[0], log[i].label[1]
            lbl.append(sput[i].label[0])
            e_sort.append(sput[i].fluence)
            yld_sort.append(sput[i].totYld)
    for i in [1,3,len(e_sort)-1]:
        c=next(c1)
        ax1.plot(e_sort[i],yld_sort[i],label=lbl[i],c=c,lw=2,ls=ls)
    fig = plt.figure(nf*2)
    ax2 = fig.add_subplot(111)
    for k in range(len(sput[i].Flux)):
        c=next(c2)
        lbl=sput[i].label[4][k]
        ax2.plot(sput[4].fluence,sput[4].Flux[k],label=lbl,c=c,lw=2,ls=ls)

    ax1.set_xlim(0,sput[0].fluence[-1])
    ax1.set_xlabel('Fluence (x10^16)')
    ax1.set_ylabel('Sputter Yield (atoms/ion)')
    IT_pair = sput[0].label[0].split()
    ax1.set_title('{}, {}, {}, {}'.format(IT_pair[1],sput[0].label[1],sput[0].label[2],sput[0].label[3]))
    ax2.set_xlim(0,sput[0].fluence[-1])
    ax2.set_xlabel('Fluence (x10^16)')
    ax2.set_ylabel('Sputter Yield (atoms/ion)')
    ax2.set_title(sput[4].label[0])
    if labl!= None:
        leg=ax1.legend()
        leg=ax2.legend()
    return
