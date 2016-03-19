import numpy as np
import matplotlib.pyplot as plt
import itertools


def plot_sputExpt(log, nf, ls, lbl=None):
    fig = plt.figure(nf)
    fig.suptitle("The simulation run is {0}".format(log[0].label[nf]))
    ax1 = fig.add_subplot(111)
    ax1.scatter(log[0].energy,log[0].totYld[nf])#, label = lbl, linestyle=ls)
    ax1.set_ylabel('Yield (atom/ion)')
    ax1.set_xlabel('Energy (eV)')
    return

def plot_log(log, nf, ls, lbl=None):
    isbv_pos=[1,3,5,7]
    color=iter(plt.cm.rainbow(np.linspace(0,1,len(isbv_pos))))
    e_sort = [[] for i in range(len(isbv_pos))]
    yld_sort = [[] for i in range(len(isbv_pos))]
    for i in range(len(log)):
        isbv=log[i].label[1]
 #       print log[i].label
        for j in range(len(isbv_pos)):
            match="isbv={}".format(isbv_pos[j])
            if match==isbv:
                e_sort[j].append(log[i].energy)
                yld_sort[j].append(log[i].totYld)

    for j in range(len(isbv_pos)):
        yld_sort[j] = [x for (y,x) in sorted(zip(e_sort[j],yld_sort[j]), key=lambda pair:pair[0])]
        e_sort[j].sort()
#    print e_sort[0], yld_sort[0]
    fig = plt.figure(nf)
#    fig.suptitle("The simulation run is {0}".format(descrip))
    ax1 = fig.add_subplot(111)
    for j in range(len(isbv_pos)):
        if e_sort[j] != []:
            if lbl != None:
                lbl = 'isbv={}'.format(isbv_pos[j])
            c=next(color)
            ax1.semilogy(e_sort[j],yld_sort[j], c=c, label = lbl, linestyle=ls)    
    leg = ax1.legend()
    return

def plot_flu(sput, nf): 
    fig = plt.figure(nf+1)
    ax1 = fig.add_subplot(111)
    asdf=np.array(range(len(sput)))
    color=iter(plt.cm.rainbow(np.linspace(0,1,len(sput))))
    ls = itertools.cycle(('--','-.',':'))
    for i in range(len(sput)):
        c=next(color)
        ax1.plot(sput[i].fluence,sput[i].totYld, label=sput[i].label[0], c=c, lw=2, ls='-')
        for k in range(len(sput[i].Flux)):
            ax1.plot(sput[i].fluence,sput[i].Flux[k], c=c, linewidth=1, ls=ls.next())

    ax1.set_xlim(0,sput[0].fluence[-1])
    ax1.set_xlabel('Fluence (x10^16)')
    ax1.set_ylabel('Sputtered Fluence')
    leg=ax1.legend()
    return
