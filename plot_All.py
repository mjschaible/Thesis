#!/usr/bin/Python

import math
import numpy as np
import matplotlib.pyplot as plt
import glob

dict = {}
energies = []
ArFlux = []
TaFlux = []
OFlux = []
path = "./*.dat"

for fn in glob.glob(path):
#    if os.path.isfile(fn):
    print('The file is ' + fn)

    f=open(fn, 'r')

    headers = {}
    for i in range(0,6):
        line = f.readline()
        line = line.strip()
        headers['x{0}'.format(i)] = line
        
    simName = f.readline()
    simName = simName.strip()
        
    print (simName + " simulation run using " + headers['x0'])
        
    elementNum=[]
    incidentF=[]
    reflectedF=[]
    reemittedF=[]
    sputteredF=[]
    depositedF=[]
    readData = 200
        
    for j, line in enumerate(f):
        line = line.strip()
#        print j, line
        if j == 15:
            columns = line.split()
            energy = float(columns[1])

        if line == 'cpt    incident   reflected   reemitted   sputtered   deposited':
            readData = int(j)
            print "The sputter data starts on line ", readData+1

        if j in range(readData+1, readData+4):
            columns = line.split()
            elementNum.append(int(columns[0]))
            incidentF.append(float(columns[1]))
            reflectedF.append(float(columns[2]))
            reemittedF.append(float(columns[3]))
            sputteredF.append(float(columns[4]))
            depositedF.append(float(columns[5]))
            print "Line ", j, ", Sputter Flux = ", columns[4]

    f.close()
    name = str(energy) + 'eV'
    print "The incident ion energy is " + name
    energies.append(energy)
    ArFlux.append(sputteredF[0])
    TaFlux.append(sputteredF[1])
    OFlux.append(sputteredF[2])
    dict[name] = sputteredF


#energies = dict.keys()
sputFlux = dict.values()
#print energies, sputFlux[1]
                        
fig = plt.figure()
                        
ax1 = fig.add_subplot(111)
                     
#ax1.set_title("Sputtered Fluence")    
ax1.set_xlabel('Energy')
ax1.set_ylabel('Sputtered Fluence')
                        
ax1.scatter(energies,OFlux, c='r', label='O flux')
ax1.scatter(energies,TaFlux, c='b', label='Ta flux')    
leg = ax1.legend()

plt.show()


