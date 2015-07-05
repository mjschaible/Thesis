#!/usr/bin/python

import math

import numpy as np
import matplotlib.pyplot as plt


kcem = [0.01, 0.1, 1.0, 2.0, 6.5] # cementation TC
kice = 6.5 # Approximate TC of Ih (hexagonal) crystalline water ice
kMout = 6.9E-4 # Effective TC at Mimas outside the anomaly
kMin = 0.0117 #  Effective TC at Mimas inside the anomaly

rg = 25E-4 # cm, grain radius (25 um)
Rcon = 0 # contact radius, depends on cementation
Rconjkr = 0 # Hertzian (JKR) contact radius due to van der Waals
Rconminlin = 0 # Minimum contact radius calculated with Wood theory
Rconminquad = 0 # 'lin' and 'quad' assume linear and quadratic dependence 
pi = math.pi

# Define a log scaled vector of cementaion volume fraction, min is zero
chi = np.logspace(1E-18,0.1,num=100,endpoint=True,base=10.0)
chimin = 0 

#Define Wood model structural parameters
Ysc = 0.09 #Obtained from curve fitting experimental data
Zsc1 = 1 #Contact radius dependence (linear vs. quadratic), TBD
Zsc2 = 2 
NC = 0 # Neighbor contacts, depends on the porosity

phi = 0.5 # Assume 50% porosity for now

#Calculate Nc [Yang et al, 2000], from curve fitting?
Nc = 2.02*((1+87.38*(1-phi)**4)/(1+25.81*(1-phi)**4))

print "The number of neighbor contacts is %.2f" % Nc

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

#Determine the Hertzian contact radius
gamma_surf = 65 # erg/cm^2
poissonH2O = 0.33 # Poissons ratio at -5C
youngs = 9E7 # erg/cm^3, Young's modulus

Fjkr = 3*pi*gamma_surf*rg # Determine the JKR contact force (var der Waals)

print "The jkr contact force is %.2f" % Fjkr, " dyne [g cm/s^2]"

# Calculate the JKR (Hertzian, ver der Waals) contact radius
inbetween = (1-poissonH2O*poissonH2O)*rg*Fjkr
inbetween2 = (0.75*inbetween)/youngs
jkrpow = float(1)/3
Rconjkr = inbetween2**(jkrpow)

print "The Hertzian contact radius is %.3e" % Rconjkr, " cm"
Rconjkrum = Rconjkr*1E4
print "The Hertzian contact radius is %.3e" % Rconjkrum, " um"

#Determine the change in contact radius with time due to irradiation
