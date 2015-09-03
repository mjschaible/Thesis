'''Program to calculate the volumetric/radial sintering rate as a function of time in order to determine the time scale at which the radiation induced sintering can increase the contact radius between grains from some initial (Hertzian, or calculated from TC outside the anomalous region) to the contact radius inside the anomalous region calculated from the Wood theory and the measured thermal inertia of the regolith'''

# dt determines the time step for the simulation.
# for log scaling (num=1e4), Tlin = 4.4E10s, Tquad=1.58E14s, Thertz=1.59E14s
#dt = np.logspace(3,15,num=1e2)
# Mimas (lin num=1e5): Tlin = 8E9s, Tquad = 3.85E13s, Thertz = 3.885E13s
# Tethys: Tlin = 1.9E7s, Tquad = 5.33E13s, Thertz = 5.33E13s
# Dione: Tlin = 188s, Tquad = 1.49E12s , Thertz = 1.49E12s


from decimal import Decimal
from decimal import getcontext
import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
kb = 1.38e-16 # erg/cm, Boltzman constant
stoyear = 3.16888e-8 # multiply seconds by this to get years
invkb = 1/kb

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

# Perform calculations for which body?
kin = kMin	
kout = kMout
Phi_elec = Phi_elec_Mimas

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

# ------ Define Sirono and Yamamoto model parameters ------
g = 4
pc = 0.333
p = 6*(1-phi)/pi
SYpf = ((1-pc)/(p-pc))*g*rg*rg/pi

def find_nearest(Rcon,Rcur):
    index = (np.abs(Rcon-Rcur)).argmin()
    return index

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
ltSij = fix_coeff(T) # Make call to calculate mechanical coefficients

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
'''print "Chiminlin = %.3e" % (chiminlin[1]/100), " um"
print "Chiminquad = %.3e" % (chiminquad[1]/100), " um"

fig = plt.figure(1)
ax1 = fig.add_subplot(2,1,1)
#ax1.semilogy(T,Rconjkr*1E4, label = 'Rcon,jkr')
ax1.semilogy(T,Rconminlin*1e4, label = 'Rcon,linear')
ax1.semilogy(T,Rconminquad*1e4, label = 'Rcon,quad.')
ax1.set_ylabel('Contact Radius [um]')
ax1.axis([60,110,5e-3,1])
ax1.legend(loc = 3, prop={'size':14})

ax2= fig.add_subplot(2,1,2)
ax2.semilogy(T,Dsurf, label = 'Dsurf')
ax2.semilogy(T,DlatB, label = 'Dlat')
ax2.semilogy(T,Dgb, label = 'Dgb')
#ax2.semilogy(T,Dradsurf, label = 'Dradsurf')
#ax2.semilogy(T,Dradlat, label = 'Dradlat')
#ax2.semilogy(T,Dradsput, label = 'Dradsput')
ax2.set_xlabel('Temperature [K]')
ax2.set_ylabel('Diffusion Coeff. [cm^2/s]')
ax2.legend(loc = 4, prop={'size':14})

plt.savefig('./DiffusionRates.png',dpi=600) '''

# ----- Determine the variation in the sintering rates based on the contact radius -----
Rcon_var = np.logspace(-6, -3, num=1000, base=10.0) # in cm
invRcon_var = 1/Rcon_var

Rnc = Rcon_var**2/(2*(rg-Rcon_var))
invRnc = 1/Rnc
theta = np.arctan(rg/(Rcon_var-Rnc))
Km = invRcon_var - invRnc
K3 = 2/rg

# ----- Calculate radiation induced diffusion coefficients ----
switches = [1,2,3]
for dzswitch in switches:
    if dzswitch == 1:
        lineT = '-'
        Rmin1=Rconminlin[0]
        Rmin2=Rconminquad[0]
        Rmax1=Rconmaxlin[0]
        Rmax2=Rconmaxquad[0]
        Dsurfc=Dsurf[0]
        Dgbc=Dgb[0]
        DlatBc=DlatB[0]
        invTc=invT[0]
        Tc=T[0]
        Pevc=Pev[0]
#        A = Dgb*delta_gb/(Dsurf*delta_surf)
        
    if dzswitch == 2:
        lineT = '--'
        Rmin1=Rconminlin[20]
        Rmin2=Rconminquad[20]
        Rmax1=Rconmaxlin[20]
        Rmax2=Rconmaxquad[20]
        Dsurfc=Dsurf[20]
        Dgbc=Dgb[20]
        DlatBc=DlatB[20]
        invTc=invT[20]
        Tc=T[20]
        Pevc=Pev[20]
  #       A = Dgb*delta_gb/(Dsurf*delta_surf)
        
    if dzswitch == 3:
        lineT = '-.'
        Rmin1=Rconminlin[40]
        Rmin2=Rconminquad[40]
        Rmax1=Rconmaxlin[40]
        Rmax2=Rconmaxquad[40]
        Dsurfc=Dsurf[40]
        Dgbc=Dgb[40]
        DlatBc=DlatB[40]
        invTc=invT[40]
        Tc=T[40]
        Pevc=Pev[40]
     
    A = Dgbc*delta_gb/(Dsurfc*delta_surf)
    
    K1 = np.zeros_like(Rcon_var)
    K2 = np.zeros_like(Rcon_var)
    d2 = np.zeros_like(Rcon_var)
    d1 = np.zeros_like(Rcon_var)

    for i in range(len(Rcon_var)):
        getcontext().prec = 20
        DKK1pos = Decimal(3*A*Rnc[i]*invRcon_var[i])+Decimal(0)
        ARR = A*Rnc[i]*invRcon_var[i]
        DKK1mid = Decimal(DKK1pos*DKK1pos)-Decimal(6.75*ARR*ARR)-Decimal(9*ARR)
        DKK1 = np.sqrt(DKK1mid+Decimal('4.5'))
#        DKK1neg = Decimal(3*(A*Rnc[i]*invRcon_var[i])**2 + 4*A*Rnc[i]*invRcon_var[i]+1).sqrt()
#        DKK1 = Decimal(DKK1pos)+Decimal('1.5')-Decimal('1.5')*Decimal(DKK1neg)
        K1[i] = Decimal(Km[i])/np.sqrt(1 + Decimal(DKK1))
        K2[i] = Decimal(K1[i])*(1 + Decimal(DKK1))
        d2[i] = Decimal(Rnc[i])/(1 + np.sqrt(Decimal(4/3)*(Decimal((K2[i]-K1[i])/K2[i]))))
        d1[i] = Decimal(Rnc[i]) - Decimal(d2[i])
        
    Vol_to_rad = 2*pi*theta*Rcon_var*Rnc
# ----- Calculate the thermal volumetric sintering rates -----
    dVdtsurf = 3*pi*Rcon_var*Dsurfc*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invTc/d2
    dVdtlat = 3*pi*Rcon_var*DlatBc*gamma_sv*Omega*(K3-Km)*invkb*invTc
    vapDiff = Pevc*Rnc*theta*np.sqrt(Omega*invkb*invTc/(2*pi*rho_ice))
    dVdtvap = 2*pi*Rcon_var*vapDiff*gamma_sv*Omega*invkb*invTc*(K3-Km)

    dVdttherm_tot = dVdtsurf + dVdtlat + dVdtvap
    dRdttherm_tot = dVdttherm_tot/Vol_to_rad

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

    dVdtradsurf = 3*pi*Rcon_var*Dradsurf*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invTc/d2
    dVdtradlat = 3*pi*Rcon_var*Dradlat*gamma_sv*Omega*(K3-Km)*invkb*invTc
    dVdtradsput = 2*pi*Rcon_var*Dradsput*gamma_sv*Omega*invkb*invTc*(K3-Km)

    dRdtradsurf=dVdtradsurf/Vol_to_rad
    dRdtradlat=dVdtradlat/Vol_to_rad
    dRdtradsput=dVdtradsput/Vol_to_rad
    
    dVdtrad_tot = dVdtradsurf + dVdtradlat + dVdtradsput 
    dRdtrad_tot = dVdtrad_tot/Vol_to_rad

    Rcur = Rmin2
 #   print 'The min rad is %.2e' % Rmin2, ' and the max is %.2e' % Rmax2
    ttot=0
    R_t = []
    timesteps = []
    dt = 3.157e7
#    print 'Rcon, dRdtrad_tot, dRdttherm_tot'
#    for i in range(10):
    while Rcur<Rmax2:
        index = find_nearest(Rcon_var,Rcur)
        DR = dRdtrad_tot[index]+dRdttherm_tot[index]
        Rcur+=DR*dt
        R_t.append(Rcur)
        ttot+=dt
        timesteps.append(ttot)
#        print '%.3e' % Rcur, '%.3e' % dRdtrad_tot[index], '%.3e' % dRdttherm_tot[index]
        
    R_tarr=np.array(R_t,dtype='float')
    ttot_years = ttot*stoyear

    print 'the total time for %d K' % Tc, ' is %.2e' % ttot_years, ' years' 

    fig = plt.figure(2)
    ax2 = fig.add_subplot(2,1,2)
    ax2.loglog(Rcon_var*1e4,dRdttherm_tot, linewidth=2, color='r', label = 'Therm. total', linestyle = lineT)
    ax2.loglog(Rcon_var*1e4,dRdtrad_tot, linewidth=2, color='b', label = 'Rad. total', linestyle = lineT)

    ax3 = fig.add_subplot(2,1,1)
    ax3.semilogx([Rmin1*1e4,Rmin1*1e4],[0,1], color='k', label= 'T = %d' % Tc, linewidth=2, linestyle = lineT)
    ax3.semilogx([Rmax1*1e4,Rmax1*1e4],[0,1], color='k', linewidth=2, linestyle = lineT)
    ax3.semilogx([Rmin2*1e4,Rmin2*1e4],[0,1], color='b', linewidth=2, linestyle = lineT)
    ax3.semilogx([Rmax2*1e4,Rmax2*1e4],[0,1], color='b', linewidth=2, linestyle = lineT)
    ax3.legend(prop={'size':18},loc=4, bbox_to_anchor=(1,-3.5))
 
    if dzswitch == 1:
        ax2.loglog(Rcon_var*1e4,dRdtradsurf, color='m', label = 'Rad. surface', linestyle = lineT)
        ax2.loglog(Rcon_var*1e4,dRdtradlat, color = 'g', label = 'Rad. lattice', linestyle = lineT)
        ax2.loglog(Rcon_var*1e4,dRdtradsput, color='c', label = 'Rad. sputter', linestyle = lineT)

        box = ax2.get_position()
        ax2.set_position([box.x0,box.y0,box.width,box.height*1.8])
        ax2.set_xlabel('Contact Radius [um]')
        ax2.set_ylabel('Radial Sintering Rate [cm/s]')
        ax2.legend(loc = 3, prop={'size':18})
        
        box = ax3.get_position()
        ax3.set_position([box.x0,box.y0+0.275,box.width,box.height*0.2])
        plt.setp(ax3.get_yticklabels(), visible=False)
        ax3.xaxis.tick_top()
        ax3.xaxis.set_label_position('top')
        ax3.set_xlabel('Contact Radius [um]')

    
    if dzswitch == 3:
        plt.savefig('./Rad_Sint_Rates.png',dpi=600)

    fig = plt.figure(3)
    ax = fig.add_subplot(1,1,1)
    line, = ax.plot(timesteps,R_tarr*1e4, linewidth=2, color='r', label = 'Rad. sint @ %d' %Tc, linestyle = lineT)
    ax.annotate('Rmax @ %d K' % Tc,xy=(ttot,Rcur*1e4),xytext=(ttot,Rcur+2.4),arrowprops=dict(facecolor='black',shrink=0.05), horizontalalignment='center')
    if dzswitch == 1:
        ax.set_ylabel('Contact Radius [um]')
        ax.set_xlabel('Timestep [s]')
    if dzswitch == 3:
        ax.legend(loc = 3, prop={'size':18})
        plt.savefig('./sint_timescales.png',dpi=600)

#plt.show()

