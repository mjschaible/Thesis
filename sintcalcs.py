#!/usr/bin/python

import math

import numpy as np
import matplotlib.pyplot as plt
pi = math.pi

# Define several meaningful surface parameters
kcem = [0.01, 0.1, 1.0, 2.0, 6.5] # cementation TC
kice = 6.5 # Approximate TC of Ih (hexagonal) crystalline water ice
kMout = 6.9E-4 # Effective TC at Mimas outside the anomaly
kMin = 0.0117 #  Effective TC at Mimas inside the anomaly
rg = 25E-4 # cm, grain radius (25 um)
phi = 0.5 # Assume 50% porosity for now

# Define a log scaled vector of cementaion volume fraction, min is zero
chi = np.logspace(1E-18,0.1,num=100,endpoint=True,base=10.0)
chimin = 0

# Define a list of temperatures at which to evaluate the sintering equations
temp = np.linspace(40,100)
#for i in range (40, 100):
#    temp.append(i)

#Define Wood model structural parameters
Ysc = 0.09 #Obtained from curve fitting experimental data
Zsc1 = 1 #Contact radius dependence (linear vs. quadratic), TBD
Zsc2 = 2 
NC = 0 # Neighbor contacts, depends on the porosity

#Calculate Nc [Yang et al, 2000], from curve fitting?
Nc = 2.02*((1+87.38*(1-phi)**4)/(1+25.81*(1-phi)**4))

print "The number of neighbor contacts is %.2f" % Nc

Rcon = 0 # contact radius, depends on cementation
Rconjkr = 0 # Hertzian (JKR) contact radius due to van der Waals
Rconminlin = 0 # Minimum contact radius calculated with Wood theory
Rconminquad = 0 # 'lin' and 'quad' assume linear and quadratic dependence 

#Determine the Wood TC eqn pre-factor
Apflin = (1/(Ysc*Nc))*(2*math.sqrt(Nc-1)*rg/Nc)**(Zsc1)
Apfquad = (1/(Ysc*Nc))*(2*math.sqrt(Nc-1)*rg/Nc)**(Zsc2)
#Calculate linear and quadratic minimum contact radius
Rconminlin = (Apflin*.5*kMout*(2+phi)/(kice*(1-phi)))
Rconminquad = (Apfquad*.5*kMout*(2+phi)/(kice*(1-phi)))**(.5)
chiminlin = (Rconminlin/(rg*0.4739))**4
chiminquad = (Rconminquad/(rg*.4739))**4
print "Rconminlin = %.3f" % (Rconminlin*1E4), " um"
print "Rconminquad = %.3f" % (Rconminquad*1E4), " um"
print "Chiminlin = %.3e" % (chiminlin/100), " um"
print "Chiminquad = %.3e" % (chiminquad/100), " um"


