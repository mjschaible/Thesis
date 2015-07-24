'''Program to calculate the volumetric/radial sintering rate as a function of time in order to determine the time scale at which the radiation induced sintering can increase the contact radius between grains from some initial (Hertzian, or calculated from TC outside the anomalous region) to the contact radius inside the anomalous region calculated from the Wood theory and the measured thermal inertia of the regolith'''

from decimal import Decimal
from decimal import getcontext
import numpy as np
import matplotlib.pyplot as plt
import csv

pi = np.pi
kb = 1.38e-16 # erg/cm, Boltzman constant
invkb = 1/kb

Escale = 1e6 # Convert bar to erg/cm^3 (for mechanical constants)

# Define several global (moon=Mimas) parameters
kMout = 69 # erg/cm/s/K, Effective TC at Mimas outside the anomaly
kMin = 1170 # erg/cm/s/K, Effective TC at Mimas inside the anomaly
rg = 25E-4 # cm, approximate grain radius (25 um)
phi = 0.5 # Assume 50% porosity for now
Phi_elec = 2e3 # elec/cm^2/s, approximate incident electron flux
dEdz = (1e6)*1.602e-12 # (eV)erg/cm, energy deposition rate

# Define water ice parameters
gamma_ss = 65 # erg/cm^2, ice/ice (grain boundary) surface energy
gamma_sv = 109 # erg/cm^2, ice/vapor (free surface) surface energy
delta_gb = 9.04e-8 # cm, grain boundary thickness
delta_surf = 3.19e-8 # cm, effective surface (atomic layer) thickness (cube root of Omega)
Omega = 3.25e-23 # cm^3, molar volume
rho_ice = 0.934 # g/cm^3, density of water ice @ 80K
WH2O = (27)*1.602e-12 # (eV)erg, average energy deposited in an excitation event

radconst = dEdz*Phi_elec/WH2O

# Calculate the approximate length scale of an electron excitation event 
Ei = 5 # eV, heating energy to lattice
U = 0.6 # eV/atom, characteristic energy of the solid
Dz = 0.25*delta_surf*(2*Ei/U -1) # cm, length-scale of excitation event

# ------ Define Wood model structural parameters -----
Ysc = 0.09 # Obtained from curve fitting experimental data
Zsc1 = 1 # Contact radius dependence (linear vs. quadratic), TBD....
Zsc2 = 2 
# Calculate Nc = Neighbor contacts, depends on the porosity from curve fitting  [Yang et al, 2000]
Nc = 2.02*((1+87.38*(1-phi)**4)/(1+25.81*(1-phi)**4))
Rconminlin = 0 # Minimum contact radius calculated with Wood theory
Rconminquad = 0 # 'lin' and 'quad' assume linear and quadratic dependence 
#Determine the Wood TC eqn pre-factor
Apflin = (1/(Ysc*Nc))*(2*np.sqrt(Nc-1)*rg/Nc)**(Zsc1)
Apfquad = (1/(Ysc*Nc))*(2*np.sqrt(Nc-1)*rg/Nc)**(Zsc2)

# generate a logrithmically scaled array where the initial values (integers) are scaled more rapidly than the later
def gen_log_space(limit, n):
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    return np.array(map(lambda x: round(x)-1, result), dtype=np.uint64)

def fix_coeff(T):
    Aij = np.zeros((5,3))

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

def calc_Rcon(dt, Rcon, n):
    invRcon = 1/Rcon
    K3 = 2/rg
    A = Dgb*delta_gb/(Dsurf*delta_surf)
    dt_prev = 0
    Rconmin = Rcon

    if n == 3:
        Rconmax = Rconmaxlin
    else:
        Rconmax = Rconmaxquad

    for j in range(len(dt)):
        Rnc = Rcon**2/(2*(rg-Rcon))
        invRnc = 1/Rnc
        theta = np.arctan(rg/(Rcon-Rnc))
        Km = invRcon - invRnc

        K1 = np.zeros_like(T)
        K2 = np.zeros_like(T)
        d2 = np.zeros_like(T)
        d1 = np.zeros_like(T)
        for i in range(len(T)):
            getcontext().prec = 20
            DKK1pos = Decimal(3*A[i]*Rnc[i]*invRcon[i])+Decimal(0)
            DKK1neg = Decimal(3*(A[i]*Rnc[i]*invRcon[i])**2 + 4*A[i]*Rnc[i]*invRcon[i]+1).sqrt()
            DKK1 = abs(Decimal(DKK1pos)+Decimal('1.5')-Decimal('1.5')*Decimal(DKK1neg))
            K1[i] = Decimal(Km[i])/np.sqrt(1 + Decimal(DKK1))
            K2[i] = Decimal(K1[i])*(1 + Decimal(DKK1))
            d2[i] = Decimal(Rnc[i])/(1 + np.sqrt(Decimal(4/3)*(Decimal(K2[i])-Decimal(K1[i]))/Decimal(K2[i])))
            d1[i] = Decimal(Rnc[i]) - Decimal(d2[i])
        
    
        Vol_to_rad = 2*pi*theta*Rcon*Rnc
        # ----- Calculate the thermal volumetric sintering rates -----
        dVdtsurf = 3*pi*Rcon*Dsurf*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invT/d2
        dVdtlat = 3*pi*Rcon*DlatB*gamma_sv*Omega*(K3-Km)*invkb*invT
        vapDiff = Pev*Rnc*theta*np.sqrt(Omega*invkb*invT/(2*pi*rho_ice))
        dVdtvap = 2*pi*Rcon*vapDiff*gamma_sv*Omega*invkb*invT*(K3-Km)

        dVdttherm_tot = dVdtsurf + dVdtlat + dVdtvap
        dRdt_therm = dVdttherm_tot/Vol_to_rad

        # ----- Calculate radiation induced diffusion coefficients ----
        Veff_surf = 4*pi*delta_surf*Rcon*Dz
        Veff_lat = 4*pi*Rcon*Dz*Dz
        Veff_sput = Dz*rg*rg*(0.25-0.125*(np.sqrt(rg*rg-Rcon*Rcon)+rg)/rg)
            
        inv_tau_surf = radconst*Veff_surf
        inv_tau_lat = radconst*Veff_lat
        inv_tau_sput = radconst*Veff_sput
            
        Dradsurf = Dz*Dz*inv_tau_surf
        Dradlat = Dz*Dz*inv_tau_lat
        Dradsput = Dz*Dz*inv_tau_sput
        dVdtradsurf = 3*pi*Rcon*Dradsurf*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invT/d2
        dVdtradlat = 3*pi*Rcon*Dradlat*gamma_sv*Omega*(K3-Km)*invkb*invT
        dVdtradsput = 2*pi*Rcon*Dradsput*gamma_sv*Omega*invkb*invT*(K3-Km)
    
        dVdtrad_tot = dVdtradsurf + dVdtradlat + dVdtradsput 
        dRdt_rad = dVdtrad_tot/Vol_to_rad

        # ----- Calculate radius increase using ALL sintering mechanisms -----
        dVdt_tot = dVdtrad_tot + dVdttherm_tot
        dRdt_tot = dVdt_tot/Vol_to_rad
        
        for k in range(len(T)):
            Rcon_new[j,k] = Rcon[k] + dRdt_tot[k] * (dt[j]-dt_prev)/3.171

        Rcon = Rcon_new[j,:]

        dt_prev = dt[j]
        invRcon = 1/Rcon
        
        if round(abs(Rcon[20]-Rconmax[20]), 7) == 0:
            print dt[j]
        
    fig = plt.figure(1)
    ax1 = fig.add_subplot(1, 3, n)
    ax1.semilogy(T,Rconmin*1E4, linewidth=2, color='red', label = 'Rcon, out')
    ax1.semilogy(T,Rcon_new[len(dt)/10000,:]*1e4, label= 't=1 kyr')
    ax1.semilogy(T,Rcon_new[len(dt)/100,:]*1e4, label= 't=100 kyr')
    ax1.semilogy(T,Rcon_new[len(dt)/10,:]*1e4, label= 't=1 Myr')
    ax1.semilogy(T,Rcon_new[len(dt)-1,:]*1e4, label= 't=10 Myr')
    ax1.semilogy(T,Rconmax*1e4, linewidth=2, color='blue', label= 'Rcon, in')
    ax1.set_xlim(60,105)
    major_ticks=np.arange(60,101,10)
    ax1.set_xticks(major_ticks)
    box = ax1.get_position()
    ax1.set_position([box.x0,box.y0+0.05,box.width,box.height*0.95])
    if n == 1:        
        ax1.set_title('Hertzian')
        ax1.set_ylabel('Contact Radius [um]')
        ax1.set_ylim(1e-1,10)
    if n == 2:
        ax1.set_title('Wood Quadratic')
        plt.setp(ax1.get_yticklabels(), visible=False)
        ax1.set_xlabel('Temperature[K]')
        ax1.set_ylim(1e-1,10)
    if n == 3:
        ax1.set_title('Wood Linear')
        ax1.yaxis.tick_right()
        ax1.set_ylim(1e-2,1)
 #       plt.setp(ax1.get_yticklabels(), visible=False)
        plt.legend(loc = 'upper center', bbox_to_anchor=(-.8,-0.1), ncol=6, prop={'size':10})
    plt.savefig('./Rcon_vs_t.png',dpi=600) 
    return Rcon_new, dVdt_rad, dVdt_therm

# ----------- Begin main program -------------------------------------------

T = np.linspace(60,110,num=51) # Define temperature range for calculations
invT = 1/T # Calculate the inverse of temperature
Tnorm = T/273 # Calculate the normalized temprature (normalize to melting temp)

# ----- Determine the thermal conductivity of ice -----
kcem = [0.01, 0.1, 1.0, 2.0, 6.5] # [J/m/s/K], cementation TC
kcem = map(lambda x: x*1e5, kcem)  # convert to [erg/cm/s/K]
kxice = (488.19/T+0.4685)*1e5 # [erg/cm/s/K], TC of Ih (hexagonal) crystalline water ice
kaice = 7.1e2*T # [erg/cm/s/K], TC of amorphous water ice

# ----- Calculate thermal diffusion coefficients for ice -----
invR = 1/8.314 # [J/mol/K], gas constant
#DsurfIc = 1.74e5 * np.exp(-38.2e3*invR*invT) # cm^2/s, from Kouchi 1994
Dsurf = 3e2 * np.exp(-44.13e3*invR*invT) # cm^2/s, from Maeno and Ebinume, 1983
DlatLDA = 7e-6*np.exp(-15e3*invR*invT) # from Ghesquiere 2015
DlatIh = 1e-10*np.exp(-9e3*invR*invT) # from Ghesquiere 2015
DlatB = 4.2e8*np.exp(-71.1e3*invR*invT) # from Livingston 1998
Dgb = 8.4*np.exp(-49e3*invR*invT) # from ....

# ----- Calculate the equilibrium vapor pressure of ice -----
Pev = 3.65e13 * np.exp(-6141.7 * invT)

# ----- Determine the compliance coefficients -----
ltSij = np.zeros((6,6,len(T)))
ltSij = fix_coeff(T)

# ----- Calculate Young's modulus, Poisson's ratio, and Hertzian radius
Ymod = np.zeros_like(T)
Prat = np.zeros_like(T)
Rconjkr = np.zeros_like(T)
jkrpow = float(1)/3

Fjkr = 3*pi*gamma_ss*rg # Determine the JKR contact force (var der Waals)

for i in range(len(T)):
    # Calculate Young's modulus in erg/cm^3 from compliance coefficients
    Ymod[i] = 1/ltSij[0,0,i]*Escale 
    Prat[i] = -ltSij[0,1,i]/ltSij[0,0,i]
    Rconjkr[i] = ((0.75*(1-Prat[i]*Prat[i])*rg*Fjkr)/Ymod[i])**(jkrpow)

# ----- Calculate linear and quadratic minimum contact radius from the Wood model -----
Rconminlin = (Apflin*.5*kMout*(2+phi)/(kxice*(1-phi)))
Rconmaxlin = (Apflin*.5*kMin*(2+phi)/(kxice*(1-phi)))
Rconminquad = np.sqrt(Apfquad*.5*kMout*(2+phi)/(kxice*(1-phi)))
Rconmaxquad = np.sqrt(Apfquad*.5*kMin*(2+phi)/(kxice*(1-phi)))
chiminlin = (Rconminlin/(rg*0.4739))**4
chiminquad = (Rconminquad/(rg*.4739))**4

print "Rconminlin = %.2e" % (Rconminlin[20]*1E4), " um"
print "Rconminquad = %.2e" % (Rconminquad[20]*1E4), " um"
print "Rconjkr = %.2e" % (Rconjkr[20]*1e4), " um"
'''print "Chiminlin = %.3e" % (chiminlin[1]/100), " um"
print "Chiminquad = %.3e" % (chiminquad[1]/100), " um" 

plt.figure(1)
plt.subplot(211)
plt.semilogy(T,Rconjkr*1E4, label = 'Rcon,jkr')
plt.semilogy(T,Rconminlin*1e4, label = 'Rcon,Zsc=1')
plt.semilogy(T,Rconminquad*1e4, label = 'Rcon,Zsc=2')
plt.ylabel('Contact Radius [um]')
plt.axis([60,110,8e-3,5e-1])
plt.legend() '''
# dt determines the time step for the simulation. This can be log or linear scaled.
# for log scaling (num=1e4), Tlin = 4.4E10s, Tquad=1.58E14s, Thertz=1.59E14s
#dt = np.logspace(3,15,num=1e2)*3.171
# for lin scaling (num=1e5), Tlin = s, Tquad = 1.58E14s, Thertz=1.59E14s
dt = np.linspace(1e6,1e15,num=1e5)*3.171
#print dt/3.171e8
print dt[len(dt)/10000]/3.171e8, dt[len(dt)/100]/3.171e8, dt[len(dt)/10]/3.171e8, dt[len(dt)-1]/3.171e8
Rcon_new = np.zeros((len(dt), len(T)))
dVdt_rad = np.zeros((len(dt), len(T)))
dVdt_therm = np.zeros((len(dt), len(T)))

# ----- Set initial contact radius as Hertzian (JKR) -----
n = 1
Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(dt, Rconjkr, n)
print "One time"
with open('Rconjkr.csv', 'rb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow('Rconjkr', 'dVdt_rad', 'dVdt_therm')
    for i in range(0, len(dt)):
        spamwriter.writerow(Rcon_new[i,20]

# ----- Set initial contact radius as Wood/quadratic -----
n = 2
Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(dt, Rconminquad, n)
print "Two time"

# ----- Set initial contact radius as Wood/linear -----
n = 3
dtl = np.linspace(1e4,1e8,num=7e4)*3.171
dt = np.linspace(1e8,1e15,num=3e4)*3.171
dt=np.append(dtl,dt)
print dt[len(dt)/20000]/3.171e8, dt[len(dt)/200]/3.171e8, dt[len(dt)/20]/3.171e8, dt[len(dt)-1]/3.171e8
Rcon_new = np.zeros((len(dt), len(T)))
dVdt_rad = np.zeros((len(dt), len(T)))
dVdt_therm = np.zeros((len(dt), len(T)))

Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(dt, Rconminlin, n)

print "Three time"

#    print dVdtsurf[20], dVdtlat[20], dVdtvap[20]
#    print
# print Rcon_therm
#    print Rcon[20]

'''plt.subplot(212)
plt.semilogy(T,Dsurf, label = 'Dsurf')
plt.semilogy(T,DlatB, label = 'Dlat')
plt.semilogy(T,Dgb, label = 'Dgb')
#plt.semilogy(T,Dradsurf, label = 'Dradsurf')
#plt.semilogy(T,Dradlat, label = 'Dradlat')
#plt.semilogy(T,Dradsput, label = 'Dradsput')
plt.xlabel('Temperature [K]')
plt.ylabel('Diffusion Coeff. [cm^2/s]')
plt.legend() '''


# ----- Determine the variation in the sintering rates based on the contact radius -----
Rcon_var = np.logspace(-6, -3, num=100, base=10.0)
invRcon_var = 1/Rcon_var

Rnc = Rcon_var**2/(2*(rg-Rcon_var))
invRnc = 1/Rnc
theta = np.arctan(rg/(Rcon_var-Rnc))
Km = invRcon_var - invRnc
K3 = 2/rg
A = Dgb[20]*delta_gb/(Dsurf[20]*delta_surf)

K1 = np.zeros_like(Rcon_var)
K2 = np.zeros_like(Rcon_var)
d2 = np.zeros_like(Rcon_var)
d1 = np.zeros_like(Rcon_var)
for i in range(len(Rcon_var)):
    getcontext().prec = 25
    DKK1pos = Decimal(3*A*Rnc[i]*invRcon_var[i])+Decimal(0)
    DKK1neg = Decimal(3*(A*Rnc[i]*invRcon_var[i])**2 + 4*A*Rnc[i]*invRcon_var[i]+1).sqrt()
    DKK1 = Decimal(DKK1pos)+Decimal('1.5')-Decimal('1.5')*Decimal(DKK1neg)
    K1[i] = Decimal(Km[i])/np.sqrt(1 + Decimal(DKK1))
    K2[i] = Decimal(K1[i])*(1 + Decimal(DKK1))
    d2[i] = Decimal(Rnc[i])/(1 + np.sqrt(Decimal(4/3)*(Decimal((K2[i]-K1[i])/K2[i]))))
    d1[i] = Decimal(Rnc[i]) - Decimal(d2[i])

Vol_to_rad = 2*pi*theta*Rcon_var*Rnc
# ----- Calculate the thermal volumetric sintering rates -----
dVdtsurf = 3*pi*Rcon_var*Dsurf[20]*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invT[20]/d2
dVdtlat = 3*pi*Rcon_var*DlatB[20]*gamma_sv*Omega*(K3-Km)*invkb*invT[20]
dVdtvap = 2*pi*Rcon_var*Rnc*theta*Pev[20]*gamma_sv*Omega*invkb*invT[20]*(K3-Km)*np.sqrt(Omega*invkb*invT[20]/(2*pi*rho_ice))

dVdttherm_tot = dVdtsurf + dVdtlat + dVdtvap
dRdt_therm = dVdttherm_tot/Vol_to_rad

# ----- Calculate radiation induced diffusion coefficients ----

Veff_surf = 4*pi*delta_surf*Rcon_var*Dz
Veff_lat = 4*pi*Rcon_var*Dz*Dz
Veff_sput = Dz*rg*rg*(0.25-0.125*(np.sqrt(rg*rg-Rcon_var*Rcon_var)+rg)/rg)

inv_tau_surf = radconst*Veff_surf
inv_tau_lat = radconst*Veff_lat
inv_tau_sput = radconst*Veff_sput

Dradsurf = Dz*Dz*inv_tau_surf
Dradlat = Dz*Dz*inv_tau_lat
Dradsput = Dz*Dz*inv_tau_sput
dVdtradsurf = 3*pi*Rcon_var*Dradsurf*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invT[20]/d2
dVdtradlat = 3*pi*Rcon_var*Dradlat*gamma_sv*Omega*(K3-Km)*invkb*invT[20]
dVdtradsput = 2*pi*Rcon_var*Dradsput*gamma_sv*Omega*invkb*invT[20]*(K3-Km)

dVdtrad_tot = dVdtradsurf + dVdtradlat + dVdtradsput 
dRdt_rad = dVdtrad_tot/Vol_to_rad

fig = plt.figure(2)
ax2 = fig.add_subplot(2,1,2)
ax2.loglog(Rcon_var*1e4,dVdttherm_tot, linewidth=2, color='r', label = 'Thermal total')
ax2.loglog(Rcon_var*1e4,dVdtrad_tot, linewidth=2, color='b', label = 'Radiation total')
ax2.loglog(Rcon_var*1e4,dVdtradsurf, color='m', label = 'Rad. surface')
ax2.loglog(Rcon_var*1e4,dVdtradlat, color = 'g', label = 'Rad. lattice')
ax2.loglog(Rcon_var*1e4,dVdtradsput, color='c', label = 'Rad. sputter')
box = ax2.get_position()
ax2.set_position([box.x0,box.y0,box.width,box.height*1.9])
ax2.set_xlabel('Contact Radius [um]')
ax2.set_ylabel('Volume Flow Rate [cm^3/s]')
ax2.legend(loc = 3, prop={'size':10})

ax3 = fig.add_subplot(2,1,1)
ax3.semilogx([Rconminlin[20]*1e4,Rconminlin[20]*1e4],[1e-34,1e-24], linestyle='--', color='blue', label= 'Rcon, out (Zsc=1)', linewidth=2)
ax3.semilogx([Rconmaxlin[20]*1e4,Rconmaxlin[20]*1e4],[1e-34,1e-24], linestyle=':', color='blue', label= 'Rcon, in (Zsc=1)', linewidth=2)
ax3.semilogx([Rconminquad[20]*1e4,Rconminquad[20]*1e4],[1e-34,1e-24], linestyle='--', color='green', label= 'Rcon, out (Zsc=2)', linewidth=2)
ax3.semilogx([Rconmaxquad[20]*1e4,Rconmaxquad[20]*1e4],[1e-34,1e-24], linestyle=':', color='green', label= 'Rcon, in (Zsc=2)', linewidth=2)
ax3.semilogx([Rconjkr[20]*1e4,Rconjkr[20]*1e4],[1e-34,1e-24], linestyle='--', color='red', label= 'Rcon,Hertzian', linewidth=2)
box = ax3.get_position()
ax3.set_position([box.x0,box.y0+0.275,box.width,box.height*0.2])
ax3.legend(prop={'size':8},ncol=5, bbox_to_anchor=(-0.1,2),loc=6)
plt.setp(ax3.get_yticklabels(), visible=False)
ax3.xaxis.tick_top()

plt.savefig('./Vol_Sint_Rates.png',dpi=600) 

#plt.show()

