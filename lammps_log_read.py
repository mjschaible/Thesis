''' Program to read lammps log file and extract plot energy data '''

import numpy as np
import matplotlib.pyplot as plt

endHeader = "Step Temp"
beginTail = "Loop time"
Ltimestep = "timestep"
pairType = "pair_style"

numLines = 0

# First open and scan the file to determine the total number of lines, where the header lines end, and where the footer lines begin
cont = 1
n = 0

while cont != 0:
    n = n+1
    filename = raw_input('Enter filename: ')

    with open(filename, 'r') as logfile:
        for num, line in enumerate(logfile, 0):
            line = line.strip()
            if endHeader in line:
                nHeader = num
                HeadCol = line.split(' ')
                #print num, line
            if beginTail in line:
                FootLine = num
                #print num, line
            if Ltimestep in line:
                columns = line.split(None)
                timestep = columns[1]
            if pairType in line:
                PairCol = line.split(None)
                PairPotent = PairCol[1]
            numLines = numLines + 1
    nFooter = numLines - FootLine

    logfile.close()

    contents = np.genfromtxt(filename, unpack=True, skip_header = nHeader+1, skip_footer = nFooter)
    #print contents[2,:]

    fig = plt.figure(1)
    plt.suptitle("The pair potential used is %s" % PairPotent)
    if n > 1:
        alen = len(fig.axes)
        for i in range(alen):
            fig.axes[i].change_geometry(n+1,1,i+1)
                   
    (ax1, ax2, ax3) = fig.add_subplot(3, n, n+1, sharex=True)
    ax1.plot(contents[0,:],contents[2,:], label = HeadCol[2])
    #ax1.set_ylim(-100000,0)
    ax1.set_title(HeadCol[2])
    ax2.plot(contents[0,:],contents[3,:], label = HeadCol[3])
    ax2.set_title(HeadCol[3])
    ax3.plot(contents[0,:],contents[4,:], label = HeadCol[4])
    #ax3.set_ylim(-100000,0)
    ax3.set_xlabel('Step (step size = %s fs)' % timestep)
    ax3.set_title(HeadCol[4])
    if n == 1:
        ax3.set_ylabel('Energy (Kcal/mole)')
        ax2.set_ylabel('Energy (Kcal/mole)')
        ax1.set_ylabel('Energy (Kcal/mole)')
    
    cont = raw_input('Enter -0- to end')

plt.show()
