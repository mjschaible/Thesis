''' Program to read lammps log file and extract plot energy data '''

import numpy as np
import matplotlib.pyplot as plt
kcalToeV = 0.0433641 # convert kcal/mole to eV/atom

def generate_plots(n):
    a = n
    b = 3
#    print "a\t=\t%d\nb\t=\t%d\na*b\t=\t%d\nn\t=\t%d" % (a,b,a*b,n)
    fig, axs = plt.subplots(b,a,figsize=(a*b,b*b), sharex='col', sharey='row')
    axs.shape = (b,a)

    for i in range(1,n):
        axs[0,i].plot(step[i],pe[i]*kcalToeV,label = HeadCol[2])
        plt.setp(axs[0,i].get_xticklabels(), visible=False)
        axs[0,i].locator_params(axis='x', tight=True, nbins=2)

        axs[1,i].plot(step[i],ke[i]*kcalToeV,label = HeadCol[3])
        plt.setp(axs[1,i].get_xticklabels(), visible=False)
        axs[1,i].locator_params(axis='x', tight=True, nbins=2)
 
        axs[2,i].plot(step[i],etot[i]*kcalToeV,label = HeadCol[4])
        axs[2,i].set_xlabel('Step #  (%s fs)' % timestep[1])
        axs[2,i].ticklabel_format(axis='x', style = 'sci')
        axs[2,i].locator_params(axis='x', tight=True, nbins=2)

        if i == 0:
            axs[0,i].set_ylabel('%s (eV/atom)' % HeadCol[2])
            axs[1,i].set_ylabel('%s (eV/atom)' % HeadCol[3])
            axs[2,i].set_ylabel('%s (eV/atom)' % HeadCol[4])
            plt.suptitle("The pair potential used is %s" % PairPotent[i])
        if i > 0:
            plt.setp(axs[0,i].get_yticklabels(), visible=False)
            plt.setp(axs[1,i].get_yticklabels(), visible=False)
            plt.setp(axs[2,i].get_yticklabels(), visible=False)
        plt.subplots_adjust(wspace=0.001)
#        axs[2,i].set_ylim(-1700,-1600)

    plt.show()
    return fig

endHeader = "Step Temp"
beginTail = "Loop time"
Ltimestep = "timestep"
pairType = "pair_style"

numLines = 0
nops = 0
nFooter=0
nHeader = []
FootLine = []
timestep = []

# First open and scan the file to determine the total number of lines, where the header lines end, and where the footer lines begin
#filename = raw_input('Enter filename: ')
filename = 'log.lammps'

with open(filename, 'r') as logfile:
    for num, line in enumerate(logfile, 0):
        line = line.strip()
        if endHeader in line:
            nHeader.append(num)
            if nops == 0:
                HeadCol = line.split(None)
            #print num, line
        if beginTail in line:
            FootLine.append(num)
            nops = nops + 1
            #print num, line
        if Ltimestep in line:
            columns = line.split(None)
            timestep.append(columns[1])
        if pairType in line:
            PairCol = line.split(None)
            PairPotent=PairCol[1]
        numLines = numLines + 1
    
    if len(FootLine)<len(nHeader):
        FootLine.append(0)
    else:
        nFooter = 0

step=[]
temp=[]
ke=[]
pe=[]
etot=[]

#step.append([])
#contents = [[],[],[],[],[],[],[]]
for i in range(1,nops):
    nFooter = numLines-FootLine[i]+1
    print range(1,nops), nops, nHeader[i]+1, FootLine[i], nFooter
    contents = np.genfromtxt(filename, unpack=True, skip_header = nHeader[i]+1, skip_footer = nFooter)    
    print contents#[0,:]

#    step.append(contents[0,:])
#    temp.append(contents[1,:])
#    ke.append(contents[2,:])
#    pe.append(contents[3,:])
#    etot.append(contents[4,:])

logfile.close()

#a = generate_plots(nops)
#a.show()

