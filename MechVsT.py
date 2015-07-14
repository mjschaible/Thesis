''' Simple script to convert stiffness coefficients C(T) to compliance S(T). Stiffness values are taken from Proctor (1966) via. Hobbs (1974) for hexagonal water ice. They were measured to be valid from 60-110K. After computing the stiffness, Young's modulus and Poisson's ratio for an isotropic material are created. Finally, the Hertzian contact radius for 25um grains is calculated.'''

import math
import numpy as np
import matplotlib.pyplot as plt
pi = math.pi
scale = 1E6 # Convert bar to erg/cm^3

Aij = np.zeros((5,3))

def fix_coeff(Aij, T):
    if 
    A011 = 12.904E4
    A111 = -1.489E-3*A011
    A211 = -1.85E-6*A011
    
    A033 = 14.075E4
    A133 = -1.629E-3*A033
    A233 = -2.93E-6*A033
    
    A044 = 2.819E4
    A144 = -1.601E-3*A044
    A244 = -3.62E-6*A044
    
    A012 = 6.487E4
    A112 = -2.072E-3*A012
    A212 = -3.62E-6*A012
    
    A013 = 5.622E4
    A113 = -1.874E-3*A013
    A213 = 0
    else if
    A011 = 17.1E4
    A111 = -47
    A211 = -29.2E-2
    
    A033 = 18.21E4
    A133 = -42
    A233 = -32.2E-2
    
    A044 = 3.62E4
    A144 = 9
    A244 = -15.5E-2
    
    A012 = 8.51E4
    A112 = 21
    A212 = -39E-2
    
    A013 = 7.13E4
    A113 = -43
    A213 = 3E-2
    
    C11 = A011 + A111*TdC + A211*TdC*TdC
    C33 = A033 + A133*TdC + A233*TdC*TdC
    C44 = A044 + A144*TdC + A244*TdC*TdC
    C12 = A012 + A112*TdC + A212*TdC*TdC
    C13 = A013 + A113*TdC + A213*TdC*TdC
    C66 = 0.5*(C11-C12)
    C0 = np.zeros_like(T)

    Cij = np.array([[C11,C12,C13,C0,C0,C0],[C12,C11,C13,C0,C0,C0],[C13,C13,C33,C0,C0,C0],[C0,C0,C0,C44,C0,C0],[C0,C0,C0,C0,C44,C0],[C0,C0,C0,C0,C0,C66]])

    Sij = np.zeros_like(Cij)
    Ymod = np.zeros_like(T)
    Prat = np.zeros_like(T)
    S11 = np.zeros_like(T)
    S33 = np.zeros_like(T)
    S44 = np.zeros_like(T)
    S12 = np.zeros_like(T)
    S13 = np.zeros_like(T)

    for i in range(Cij.shape[2]):
        Sij[:,:,i] = np.linalg.inv(Cij[:,:,i])

    return Sij

for i in range(Cij.shape[2]):
    Ymod[i] = 1/Sij[0,0,i]*scale # Calculate Young's modulus in erg/cm^3 from compliance coefficients
    Prat[i] = -Sij[0,1,i]/Sij[0,0,i]
    S11[i] = Sij.item(0,0,i)
    S33[i] = Sij.item(2,2,i)
    S44[i] = Sij.item(3,3,i)
    S12[i] = Sij.item(0,1,i)
    S13[i] = Sij.item(0,2,i)
print
print "Poisson's ratio at ", T[139]-273, "K is", Prat[139]
print "Young's modulus at ", T[139], "K is %.2e" % Ymod[139], "erg/cm^3"


T = np.linspace(130,273,num=144) # Temp in Kelvin
TdC = T - 273 # Temp in Celcius




TT = T
T = 0
T = np.linspace(60,110,num=81)
TT = np.concatenate((T,TT))

ltC11 = A011 + A111*T + A211*T*T
ltC33 = A033 + A133*T + A233*T*T
ltC44 = A044 + A144*T + A244*T*T
ltC12 = A012 + A112*T + A212*T*T
ltC13 = A013 + A113*T + A213*T*T
ltC66 = 0.5*(ltC11-ltC12)
ltC0 = np.zeros_like(T)

ltCij = np.array([[ltC11,ltC12,ltC13,ltC0,ltC0,ltC0],[ltC12,ltC11,ltC13,ltC0,ltC0,ltC0],[ltC13,ltC13,ltC33,ltC0,ltC0,ltC0],[ltC0,ltC0,ltC0,ltC44,ltC0,ltC0],[ltC0,ltC0,ltC0,ltC0,ltC44,ltC0],[ltC0,ltC0,ltC0,ltC0,ltC0,ltC66]])

ltSij = np.zeros_like(ltCij)
ltYmod = np.zeros_like(T)
ltPrat = np.zeros_like(T)
ltS11 = np.zeros_like(T)
ltS33 = np.zeros_like(T)
ltS44 = np.zeros_like(T)
ltS12 = np.zeros_like(T)
ltS13 = np.zeros_like(T)

for i in range(ltCij.shape[2]):
    ltSij[:,:,i] = np.linalg.inv(ltCij[:,:,i])
    ltYmod[i] = 1/ltSij[0,0,i]*scale # Calculate Young's modulus in erg/cm^3 from compliance coefficients
    ltPrat[i] = -ltSij[0,1,i]/ltSij[0,0,i]
    ltS11[i] = ltSij.item(0,0,i)
    ltS33[i] = ltSij.item(2,2,i)
    ltS44[i] = ltSij.item(3,3,i)
    ltS12[i] = ltSij.item(0,1,i)
    ltS13[i] = ltSij.item(0,2,i)

Ymod = np.concatenate((ltYmod,Ymod))
Prat = np.concatenate((ltPrat,Prat))
S11 = np.concatenate((ltS11,S11))
S33 = np.concatenate((ltS33,S33))
S44 = np.concatenate((ltS44,S44))
S12 = np.concatenate((ltS12,S12))
S13 = np.concatenate((ltS13,S13))
C11 = np.concatenate((ltC11,C11))
C33 = np.concatenate((ltC33,C33))
C44 = np.concatenate((ltC44,C44))
C12 = np.concatenate((ltC12,C12))
C13 = np.concatenate((ltC13,C13))

print 
#print "The C11 coefficient at ", T[40]-273, "C is", ltCij.item(0,0,40), "bar"
#print "The S11 coefficient at ", T[40], "K is", ltSij.item(0,0,40), "bar"
#print "The S33 coefficient at ", T[40], "K is", ltSij.item(2,2,40), "bar"
#print "The S44 coefficient at ", T[40], "K is", ltSij.item(3,3,40), "bar"
#print "The S12 coefficient at ", T[40], "K is", ltSij.item(0,1,40), "bar"
#print "The S13 coefficient at ", T[40], "K is", ltSij.item(0,2,40), "bar"
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

plt.subplot(211)
plt.plot(TT,S11*scale, label = 'S11')
plt.plot(TT,S33*scale, label = 'S33')
plt.plot(TT,S44*scale, label = 'S44')
plt.plot(TT,S12*scale, label = 'S12')
plt.plot(TT,S13*scale, label = 'S13')
plt.ylabel('Compliance Coefficient \n [cm^3/erg]')
plt.legend()

plt.subplot(212)
plt.plot(TT,Prat, label = '-S12/S11') 
plt.xlabel('Temperature [K]')
plt.ylabel('Poisson\'s ratio')
plt.legend()

plt.figure(2)
plt.plot(TT,Rconjkr*1E4, label = 'Rcon,jkr')
plt.xlabel('Temperature[k]')
plt.ylabel('Contact Radius [um]')
plt.legend()

plt.show()
