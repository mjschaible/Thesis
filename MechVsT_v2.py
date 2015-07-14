''' Simple script to convert stiffness coefficients C(T) to compliance S(T). Stiffness values are taken from Proctor (1966) via. Hobbs (1974) for hexagonal water ice. They were measured to be valid from 60-110K. After computing the stiffness, Young's modulus and Poisson's ratio for an isotropic material are created. Finally, the Hertzian contact radius for 25um grains is calculated.'''

import math
import numpy as np
import matplotlib.pyplot as plt
pi = math.pi
scale = 1E6 # Convert bar to erg/cm^3

def fix_coeff(T):
    Aij = np.zeros((5,3))

    if min(T) < 0:
        Aij[0,0] = 12.904E4
        Aij[0,1] = -1.489E-3*Aij[0,0]
        Aij[0,2] = -1.85E-6*Aij[0,0]
    
        Aij[1,0] = 14.075E4
        Aij[1,1] = -1.629E-3*Aij[1,0]
        Aij[1,2] = -2.93E-6*Aij[1,0]
    
        Aij[2,0] = 2.819E4
        Aij[2,1] = -1.601E-3*Aij[2,0]
        Aij[2,2] = -3.62E-6*Aij[2,0]
    
        Aij[3,0] = 6.487E4
        Aij[3,1] = -2.072E-3*Aij[3,0]
        Aij[3,2] = -3.62E-6*Aij[3,0]
    
        Aij[4,0] = 5.622E4
        Aij[4,1] = -1.874E-3*Aij[4,0]
        Aij[4,2] = 0
    else:
        Aij[0,0] = 17.1E4
        Aij[0,1] = -47
        Aij[0,2] = -29.2E-2
    
        Aij[1,0] = 18.21E4
        Aij[1,1] = -42
        Aij[1,2] = -32.2E-2
    
        Aij[2,0] = 3.62E4
        Aij[2,1] = 9
        Aij[2,2] = -15.5E-2
    
        Aij[3,0] = 8.51E4
        Aij[3,1] = 21
        Aij[3,2] = -39E-2
    
        Aij[4,0] = 7.13E4
        Aij[4,1] = -43
        Aij[4,2] = 3E-2
    
    C11 = Aij[0,0] + Aij[0,1]*T + Aij[0,2]*T*T
    C33 = Aij[1,0] + Aij[1,1]*T + Aij[1,2]*T*T
    C44 = Aij[2,0] + Aij[2,1]*T + Aij[2,2]*T*T
    C12 = Aij[3,0] + Aij[3,1]*T + Aij[3,2]*T*T
    C13 = Aij[4,0] + Aij[4,1]*T + Aij[4,2]*T*T
    C66 = 0.5*(C11-C12)
    C0 = np.zeros_like(T)

    Cij = np.array([[C11,C12,C13,C0,C0,C0],[C12,C11,C13,C0,C0,C0],[C13,C13,C33,C0,C0,C0],[C0,C0,C0,C44,C0,C0],[C0,C0,C0,C0,C44,C0],[C0,C0,C0,C0,C0,C66]])

    tempSij = np.zeros_like(Cij)

    for i in range(len(T)):
        tempSij[:,:,i] = np.linalg.inv(Cij[:,:,i])

    return tempSij

T = np.linspace(130,273,num=144) # Temp in Kelvin
TdC = T - 273 # Temp in Celcius

Sij = np.zeros((6,6,len(T)))
Sij = fix_coeff(TdC)

Ymod = np.zeros_like(T)
Prat = np.zeros_like(T)

for i in range(len(T)):
    # Calculate Young's modulus in erg/cm^3 from compliance coefficients
    Ymod[i] = 1/Sij[0,0,i]*scale 
    Prat[i] = -Sij[0,1,i]/Sij[0,0,i]

print
print "Poisson's ratio at ", T[139]-273, "C is", Prat[139]
print "Young's modulus at ", T[139], "K is %.2e" % Ymod[139], "erg/cm^3"

TT = T
T = 0
T = np.linspace(60,110,num=81)
TT = np.concatenate((T,TT))

ltSij = np.zeros((6,6,len(T)))
ltSij = fix_coeff(T)

ltYmod = np.zeros_like(T)
ltPrat = np.zeros_like(T)

for i in range(len(T)):
    # Calculate Young's modulus in erg/cm^3 from compliance coefficients
    ltYmod[i] = 1/ltSij[0,0,i]*scale
    ltPrat[i] = -ltSij[0,1,i]/ltSij[0,0,i]

Ymod = np.concatenate((ltYmod,Ymod))
Prat = np.concatenate((ltPrat,Prat))

print 
print "Poisson's ratio at ", T[40], "K is", Prat[40]
print "Young's modulus at ", T[40], "K is %.2e" % Ymod[40], "erg/cm^3"

#Determine the Hertzian contact radius
rg = 25E-4 # cm, grain radius (25um)
gamma_surf = 65 # erg/cm^2, ice/ice (grain boundary) surface energy
#poissonH2O = 0.33 # Poissons ratio at -5C
#youngs = 9E7 # erg/cm^3, Young's modulus
Fjkr = 3*pi*gamma_surf*rg # Determine the JKR contact force (var der Waals)
print "The jkr contact force is %.2f" % Fjkr, " dyne [g cm/s^2]"

# Calculate the JKR (Hertzian, ver der Waals) contact radius
inbetween = np.zeros_like(TT)
inbetween2 = np.zeros_like(TT)
Rconjkr = np.zeros_like(TT)

for i in range(len(TT)):
    inbetween[i] = (1-Prat[i]*Prat[i])*rg*Fjkr
    inbetween2[i] = (0.75*inbetween[i])/Ymod[i]
    jkrpow = float(1)/3
    Rconjkr[i] = inbetween2[i]**(jkrpow)

print "The Hertzian contact radius at ", T[40], "K is %.3e" % Rconjkr[40], " cm"
Rconjkrum = Rconjkr*1E4
print "The Hertzian contact radius ", T[40], "K is %.3e" % Rconjkrum[40], " um"

plt.figure(1)

plt.plot(TT,Rconjkr*1E4, label = 'Rcon,jkr')
plt.xlabel('Temperature[k]')
plt.ylabel('Contact Radius [um]')
plt.legend()

#plt.show()
