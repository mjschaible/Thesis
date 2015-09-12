import numpy as np
import matplotlib.pyplot as plt

## First open and scan the file to determine the total number of lines, where the header lines end, and where the footer lines begin                

# Variable to count how many species to track
SpecHead = "by"
numSpec = 0

# Initialize a list to hold the header steps
fsteps = []
sputdat1 = []
sputdat2 = []
sputdat3 = []
# Initialize variable for number of lines in file, where the header starts, and where the footer (in any) starts
numLines = 0
nHeader = 0
nFooter = 0

# To read in from command line, use:
filename = raw_input('Enter filename: ')                                       

# To input a specific file
#filename = 'pythondata.dat'

# open, read only ('r') file. Pull out sim version, description, and parameters
with open(filename, 'r') as logfile:
    for num, line in enumerate(logfile, 0):
        # First 7 lines should ALWAYS be the same. Read in values as 'line'.
        # Strip removes whitespace/carrage return from EOL
        line = line.strip()
        if num == 0:
            version = line
        elif num == 1:
            descrip = line
        elif num == 2:
            paramLabel = line.split(None)
        elif num == 3:
            params = line.split(None)
            ncp = int(params[2]) # Determine # of particles in system
            print "The number of species is", ncp
        elif num == 4:
            species = line.split(None)
        elif num == 5:
            fluLabel = line
        elif num == 6:
            # should automatically set these to 'float'
            firststep = float(line)
            fsteps.append(firststep)
        elif num == 7:
            what_dis = line
        elif SpecHead in line:
            numSpec = numSpec+1
            nHeader = num+2
        else:
            row = line.split(None)
            num_cols = len(row)
            if num_cols == 1:
                fsteps.append(float(row[0]))
            else:
# I am cheating here. Should be able to account for an arbitrary number of species (not just sputdat1-3)
                sputdat1.append(float(row[0]))
                sputdat2.append(float(row[1]))
                sputdat3.append(float(row[2]))

        numLines = numLines + 1

logfile.close()

#print fsteps, len(fsteps)
#print sputdat

numSteps = len(fsteps)
#print sputdat3

# This needs to be defined as a sequence of ints?
specFlu = np.zeros((numSteps,ncp))

# This part also needs updated to handle an arbitrary number of sample species
for i in range(0,ncp):
    for j in range(0,numSteps):
        cur = j*3
        if i == 0:
            specFlu[j,i] = sputdat1[cur]+sputdat1[cur+1]+sputdat1[cur+2]
        elif i == 1:
            specFlu[j,i] = sputdat2[cur]+sputdat2[cur+1]+sputdat2[cur+2]
        else:
            specFlu[j,i] = sputdat3[cur]+sputdat3[cur+1]+sputdat3[cur+2]  
    print specFlu[:,i]

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

#av0 = movingaverage(specFlu[:,0],10)
#av1 = movingaverage(specFlu[:,1],10)
#av2 = movingaverage(specFlu[:,2],10)

fig = plt.figure()
fig.suptitle("The simulation run is %s" %descrip)
ax1 = fig.add_subplot(111)
ax1.plot(fsteps,specFlu[:,0], label = species[0])
#ax1.plot(fsteps,av0)
#ax1.set_xlim(.001,20)
ax1.set_ylim(.001,0.2)                    
ax1.set_ylabel('Yield (atom/ion)')
ax1.plot(fsteps,specFlu[:,1], label = species[1])
#ax1.plot(fsteps,av1)
ax1.plot(fsteps,specFlu[:,2], label = species[2])
#ax1.plot(fsteps,av2)
ax1.set_xlabel('Fluence (x10^16)')
leg = ax1.legend()

fig2 = plt.figure()
fig2.suptitle("Averaged simulation results for %s" %descrip)
ax2 = fig2.add_subplot(111)
ax2.set_ylim(0.001,0.2)
#ax2.set_xlim(0.001,20)
ax2.set_ylabel('Yield (atom/ion')
ax2.set_xlabel('Fluence (x10^16)')

av0 = movingaverage(specFlu[:,0],20)
av1 = movingaverage(specFlu[:,1],20)
av2 = movingaverage(specFlu[:,2],20)

ax2.plot(fsteps,av0, label = species[0])
ax2.plot(fsteps,av1, label = species[1])
ax2.plot(fsteps,av2, label = species[2])

leg2 = ax2.legend()




plt.show()

