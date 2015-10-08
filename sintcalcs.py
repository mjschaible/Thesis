'''Program to calculate the volumetric/radial sintering rate as a function of time in order to determine the time scale at which the radiation induced sintering can increase the contact radius between grains from some initial (Hertzian, or calculated from TC outside the anomalous region) to the contact radius inside the anomalous region calculated from the Wood theory and the measured thermal inertia of the regolith'''

from decimal import Decimal
from decimal import getcontext
import numpy as np
import matplotlib.pyplot as plt
import csv

pi = np.pi
kb = 1.38e-16 # erg/cm, Boltzman constant
invkb = 1/kb
stoy = 3.16888e-8 # convert seconds to years

Escale = 1e6 # Convert bar to erg/cm^3 (for mechanical constants)

# Define several global (moon=Mimas) parameters
kMout = 69 # erg/cm/s/K, Effective TC at Mimas outside the anomaly
kMin = 1170 # erg/cm/s/K, Effective TC at Mimas inside the anomaly
kTout = 6.53 # erg/cm/s/K, Effective TC at Tethys outside the anomaly
kTin = 163 # erg/cm/s/K, Effective TC at Tethys inside the anomaly
kDtrail = 17.1 # erg/cm/s/K, Effective TC at Dione on the trailing hemisphere
kDlead = 32.4 # erg/cm/s/K, Effective TC at Dione on the leading hemisphere

rg = 25E-4 # cm, approximate grain radius (25 um) (all bodies)
phi = 0.5 # Assume 50% porosity for now (all bodies)
Phi_elec_Mimas = 8.2e3 # elec/cm^2/s, incident electron flux at Mimas
Phi_elec_Tethys = 2.1e2 # elec/cm^2/s, incident electron flux at Tethys
Phi_elec_Dione = 2.8e2 # elec/cm^2/s, incident electron flux at Dione

# For now assume an average energy loss for all depths
# Should eventually use Penelope results to give vs. depth
dEdz = (1e6)*1.602e-12 # (eV)erg/cm, energy deposition rate at Mimas

# Define water ice parameters
gamma_ss = 65 # erg/cm^2, ice/ice (grain boundary) surface energy
gamma_sv = 109 # erg/cm^2, ice/vapor (free surface) surface energy
delta_gb = 9.04e-8 # cm, grain boundary thickness
delta_surf = 3.19e-8 # cm, effective surface (atomic layer) thickness (cube root of Omega)
Omega = 3.25e-23 # cm^3, molar volume
rho_ice = 0.934 # g/cm^3, density of water ice @ 80K
WH2O = (27)*1.602e-12 # (eV)erg, average energy deposited in an excitation event

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

# ------ Define Sirono and Yamamoto model parameters ------
g = 4
pc = 0.333
p = 6*(1-phi)/pi
SYpf = ((1-pc)/(p-pc))*g*rg*rg/pi

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

def calc_Rcon(Rcon, n, dzswitch):

    K3 = 2/rg
    A = Dgb*delta_gb/(Dsurf*delta_surf)
    dt_prev = 0

    if n  > 2:
        Rconmax = Rconmaxlin
        Rconmin = Rconminlin
        prec = 7
    if n == 1:
        Rconmax = Rconmaxquad
        Rconmin = Rconjkr
        prec = 6
#    if n == 2:
    Rconmax = Rconmaxquad
    Rconmin = Rconminquad
    prec = 7

    kyr1 = 10
    kyr10 = len(dt)*1e-2
    kyr100 = len(dt)*1e-1
    Myr1 = len(dt)-1
    print dt[kyr1]*stoy, dt[kyr10]*stoy, dt[kyr100]*stoy, dt[Myr1]*stoy

    for j in range(len(dt)):
        invRcon = 1/Rcon
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
            ARR = A[i]*Rnc[i]*invRcon[i]
            DKK1mid = Decimal(DKK1pos*DKK1pos)-Decimal(6.75*ARR*ARR)-Decimal(9*ARR)
            DKK1 = np.sqrt(DKK1mid+Decimal('4.5'))
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
            Rcon_new[j,k] = Rcon[k] + dRdt_tot[k] * (dt[j]-dt_prev)

        Rcon = Rcon_new[j,:]
        dt_prev = dt[j]

        if round(abs(Rcon[20]-Rconmax[20]), prec) == 0:
            print 'The max contact radius was reached in %.3e seconds' % dt[j]

    lineT = '-'
    
    if n < 4:
        col = dzswitch
    else:
        col = 3

    if n == 0:
        print "timestep = %d" % (dt[j]*stoy)
    elif n == 3:
        print "timestep = %d" % (dt[j]*stoy)
    elif n == 4:
        print "timestep = %d" % (dt[j]*stoy)
    else:
        fig = plt.figure(1)
        ax1 = fig.add_subplot(1, 3, col)
        ax1.semilogy(T,Rconmin*1E4, linewidth=2, color='k', label = 'Rcon, out', linestyle = '-')
        ax1.semilogy(T,Rcon_new[kyr1,:]*1e4, label= 't=1 kyr', linestyle = lineT, color='b')
        ax1.semilogy(T,Rcon_new[kyr10,:]*1e4, label= 't=10 kyr', linestyle = lineT, color='c')
        ax1.semilogy(T,Rcon_new[kyr100,:]*1e4, label= 't=100 kyr', linestyle = lineT, color='g')
        ax1.semilogy(T,Rcon_new[Myr1,:]*1e4, label= 't=1 Myr', linestyle = lineT,color='m')
        ax1.semilogy(T,Rconmax*1e4, linewidth=2, color='r', label= 'Rcon, in', linestyle = '-')
        ax1.set_xlim(60,105)
        major_ticks=np.arange(60,101,10)
        ax1.set_xticks(major_ticks)
        box = ax1.get_position()
        ax1.set_position([box.x0,box.y0+0.05,box.width,box.height*0.95])
        if dzswitch == 1:
            ax1.set_title('Mimas')
            ax1.set_ylabel('Contact Radius [um]')
#            ax1.set_ylim(1e-1,10)
        if dzswitch == 2:
            ax1.set_title('Tethys')
            plt.setp(ax1.get_yticklabels(), visible=False)
            ax1.set_xlabel('Temperature[K]')
#            ax1.set_ylim(1e-1,10)
        if dzswitch == 3:
            ax1.set_title('Dione')
            ax1.yaxis.tick_right()
#            ax1.set_ylim(1e-2,10)
         #       plt.setp(ax1.get_yticklabels(), visible=False)
            ax1.legend(loc = 'upper center', bbox_to_anchor=(-.8,-0.1), ncol=6, prop={'size':10})

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

# dt determines the time step for the simulation. This can be log or linear scaled.

#dt = np.logspace(3,15,num=1e2)
# Mimas (lin num=1e5): Tlin = 2.7E9s, Tquad = 1.216E13s, Thertz = 1.224E13s
# Tethys: Tlin = 1.9E7s, Tquad = 5.33E13s, Thertz = 5.33E13s
# Dione: Tlin = 188s, Tquad = 1.49E12s , Thertz = 1.49E12s
# for Dz = 20nm, Tlin = s, Tquad = s, Thertz = s
dt = np.linspace(0,3.154e10,num=1e4)
Rcon_new = np.zeros((len(dt), len(T)))
dVdt_rad = np.zeros((len(dt), len(T)))
dVdt_therm = np.zeros((len(dt), len(T)))

switches = [1,2,3]
for dzswitch in switches:
    dtstart = 0
    if dzswitch == 1:
    # Perform calculations for which body?
        kin = kMin
        kout = kMout
        Phi_elec = Phi_elec_Mimas
    if dzswitch == 2:
        kin = kTin
        kout = kTout
        Phi_elec = Phi_elec_Tethys
    if dzswitch == 3:
        kin = kDlead
        kout = kDtrail
        Phi_elec = Phi_elec_Dione

    radconst = dEdz*Phi_elec/WH2O

# ----- Calculate linear and quadratic minimum contact radius from the Wood model -----
    Rconminlin = (Apflin*kout*(1+0.5*phi)/(kxice*(1-phi)))
    Rconmaxlin = (Apflin*kin*(1+0.5*phi)/(kxice*(1-phi)))
    Rconminquad = np.sqrt(Apfquad*kout*(1+0.5*phi)/(kxice*(1-phi)))
    Rconmaxquad = np.sqrt(Apfquad*kin*(1+0.5*phi)/(kxice*(1-phi)))
    RconminSY = np.sqrt((kout/kxice)*SYpf)
    RconmaxSY = np.sqrt((kin/kxice)*SYpf)
    chiminlin = (Rconminlin/(rg*0.4739))**4
    chiminquad = (Rconminquad/(rg*.4739))**4
    
    print "Rconminlin = %.2e" % (Rconminlin[20]*1E4), " um"
    print "Rconmaxlin = %.2e" % (Rconmaxlin[20]*1e4), " um"
    print "Rconminquad = %.2e" % (Rconminquad[20]*1E4), " um"
    print "Rconmaxquad = %.2e" % (Rconmaxquad[20]*1e4), " um"
    print "RconSYmax = %.2e" % (RconmaxSY[20]*1e4), " um"
    print "RconSYmin = %.2e" % (RconminSY[20]*1e4), " um"
    print "Rconjkr = %.2e" % (Rconjkr[20]*1e4), " um"

    print "The interaction length scale being used is ", Dz

# ----- Set initial contact radius as Hertzian (JKR) -----
    '''
    n = 0
    dt = np.linspace(0,3.154e9,num=1e4)
    Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(Rconjkr, n, dzswitch)
    dtstart=len(dt)-1
    Rconfinal=Rcon_new[dtstart,:]

    n = 1
    dt = np.linspace(3.154e9,3.154e13,num=1e4)

    Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(Rconfinal, n, dzswitch)
    print "One time"
    with open('Rconjkr.csv', 'rb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow('Rconjkr', 'dVdt_rad', 'dVdt_therm')
        for i in range(0, len(dt)):
            spamwriter.writerow(Rcon_new[i,20]
    '''

    # ----- Set initial contact radius as Wood/quadratic -----
    n = 0
    dt = np.linspace(0,3.154e9,num=1e4)
    Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(Rconminquad, n, dzswitch)
    dtstart=len(dt)-1
    Rconfinal=Rcon_new[dtstart,:]

    n = 2
    dt = np.linspace(3.154e9,3.154e13,num=1e4)
    Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(Rconfinal, n, dzswitch)
    print "Two time"

    # ----- Set initial contact radius as Wood/linear -----
    '''
    n = 3
    dt = np.linspace(0,3.154e5,num=1e4)
    Rcon_new = np.zeros((len(dt), len(T)))
    dVdt_rad = np.zeros((len(dt), len(T)))
    dVdt_therm = np.zeros((len(dt), len(T)))
    Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(Rconminlin, n, dzswitch)
    dtstart=len(dt)-1
    Rconfinal=Rcon_new[dtstart,:]
    print "Three time"

    n = 4
    dt = np.linspace(3.154e5,3.154e9,num=1e4)
    Rcon_new = np.zeros((len(dt), len(T)))
    dVdt_rad = np.zeros((len(dt), len(T)))
    dVdt_therm = np.zeros((len(dt), len(T)))
    Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(Rconfinal, n, dzswitch)
    dtstart=len(dt)-1
    Rconfinal=Rcon_new[dtstart,:]
    print "Three time"

    n = 5
    dt = np.linspace(3.154e9,3.154e13,num=1e4)
    Rcon_new = np.zeros((len(dt), len(T)))
    dVdt_rad = np.zeros((len(dt), len(T)))
    dVdt_therm = np.zeros((len(dt), len(T)))
    Rcon_new, dVdt_rad, dVdt_therm = calc_Rcon(Rconfinal, n, dzswitch)
    print "Three time" '''
