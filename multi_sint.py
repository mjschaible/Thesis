'''Program to calculate the volumetric/radial sintering rate as a function of time in order to determine the time scale at which the radiation induced sintering can increase the contact radius between grains from some initial (Hertzian, or calculated from TC outside the anomalous region) to the contact radius inside the anomalous region calculated from the Wood theory and the measured thermal inertia of the regolith'''

from decimal import Decimal
from decimal import getcontext
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 16})

pi = np.pi
kb = 1.38e-16 # erg/K, Boltzman constant
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

rg = 25e-4 # cm, approximate grain radius (25 um) (all bodies)
roundness = 0.2
sphericity = 0.8
rloc = rg*sphericity*roundness

#phi = np.linspace(40,99,num=100)
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
lat_const = 3e-8 # cm, approximate lattice spacing
Omega = 3.25e-23 # cm^3, molar volume
rho_ice = 0.934 # g/cm^3, density of water ice @ 80K
WH2O = (27)*1.602e-12 # (eV)erg, average energy deposited in an excitation event
nH2O = 3E22 # molecules/cm^3

Nvac = 6 # number of vacancies produced by a single excited water molecule
radconst = dEdz*Phi_elec/(WH2O*nH2O)

c_geo = 4 # geometric factor to determine sputtering/deposition relation

# Calculate the approximate length scale of an electron excitation event 
Ei = 5 # eV, heating energy to lattice
U = 0.6 # eV/atom, characteristic energy of the solid
Dz = 0.25*delta_surf*(2*Ei/U -1) # cm, length-scale of excitation event assuming binary collision

# ---- Functions ----
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

# ----- Determine the variation in the sintering rates based on the contact radius -----
Rcon_tot = np.logspace(-6, -3.61, num=1000, base=10.0) # in cm
N = 6 # number of contacts
Rcon_var = Rcon_tot/N
invRcon_var = 1/Rcon_var

Rnc = Rcon_var**2/(2*(rloc-Rcon_var))
invRnc = 1/Rnc
theta = np.arctan(rg/(Rcon_var-Rnc))
Km = invRcon_var - invRnc
K3 = 2/rloc

# Calculate Nc = Neighbor contacts, depends on the porosity from curve fitting  [Yang et al, 2000]
Nc = 2.02*((1+87.38*(1-phi)**4)/(1+25.81*(1-phi)**4))

# Calculate the cementation volume fraction
chi_var_single = (3/16)*(Rcon_var/rg)**4 # cementation volume fraction of a single contact
chi_var = chi_var_single * N * Nc # cementation volume fraction of N contacts per Nc neighbors  
V_var = 2*pi*N*Rcon_var*Rcon_var*Rnc

# ----- Determine the thermal conductivity of ice -----
#kcem = [0.01, 0.1, 1.0, 2.0, 6.5] # [J/m/s/K], cementation TC
#kcem = map(lambda x: x*1e5, kcem)  # convert to [erg/cm/s/K]
kxice = (488.19/T+0.4685)*1e5 # [erg/cm/s/K], TC of Ih (hexagonal) crystalline water ice
kaice = 7.1e2*T # [erg/cm/s/K], TC of amorphous water ice
kg=kxice
kcem=kaice
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

# ------ Define Sirono and Yamamoto model parameters ------
g = 4
pc = 0.333
p = 6*(1-phi)/pi
SYpf = ((1-pc)/(p-pc))*g*rg*rg/pi
# ------ Define Wood model structural parameters -----
Ysc = 0.09 # Obtained from curve fitting experimental data
Zsc1 = 1 # Contact radius dependence (linear vs. quadratic), TBD....
Zsc2 = 2 
#Determine the Wood TC eqn pre-factor
Apflin = (1/(Ysc*Nc))*(2*np.sqrt(Nc-1)*rg/Nc)**(Zsc1)
Apfquad = (1/(Ysc*Nc))*(2*np.sqrt(Nc-1)*rg/Nc)**(Zsc2)
# ----- Calculate linear and quadratic minimum contact radius from the Wood model -----
Rconminlin = (Apflin*kout*(1+0.5*phi)/(kxice*(1-phi)))
Rconmaxlin = (Apflin*kin*(1+0.5*phi)/(kxice*(1-phi)))
Rconminquad = np.sqrt(Apfquad*kout*(1+0.5*phi)/(kxice*(1-phi)))
Rconmaxquad = np.sqrt(Apfquad*kin*(1+0.5*phi)/(kxice*(1-phi)))
RconminSY = np.sqrt((kout/kxice)*SYpf)
RconmaxSY = np.sqrt((kin/kxice)*SYpf)
chiminlin = (Rconminlin/(rg*0.4739))**4
chiminquad = Nc*3/16*(Rconminquad/rg)**4
chimaxquad = Nc*3/16*(Rconmaxquad/rg)**4
r_conminquad = Rconminquad/N
r_conmaxquad = Rconmaxquad/N 
Vminquad = pi*N*r_conminquad*r_conminquad*r_conminquad*r_conminquad/(rloc-r_conminquad)
Vmaxquad = pi*N*r_conmaxquad*r_conmaxquad*r_conmaxquad*r_conmaxquad/(rloc-r_conmaxquad)

print "Rconminlin = %.2e" % (Rconminlin[20]*1E4), " um"
print "Rconmaxlin = %.2e" % (Rconmaxlin[20]*1e4), " um"
print "Rconminquad = %.2e" % (Rconminquad[20]*1E4), " um"
print "Rconmaxquad = %.2e" % (Rconmaxquad[20]*1e4), " um"
print "RconSYmax = %.2e" % (RconmaxSY[20]*1e4), " um"
print "RconSYmin = %.2e" % (RconminSY[20]*1e4), " um"
print "Rconjkr = %.2e" % (Rconjkr[20]*1e4), " um"
print "Vminquad = %.3e" % (Vminquad[20]*1e4**3), "um^3"
print "Vmaxquad = %.3e" % (Vmaxquad[20]*1e4**3), "um^3"

#Veff_surf = 4*pi*delta_surf*Rconjkr*Dz
#Veff_lat = 2*pi*Rconjkr*Rconjkr*Dz
#inv_tau_surf = radconst*Veff_surf
#inv_tau_lat = radconst*Veff_lat
#Dradsurf = Dz*Dz*inv_tau_surf
#Dradlat = Dz*Dz*inv_tau_lat
inv_tau_vac = Nvac*radconst
Drad_eff = lat_const*lat_const*inv_tau_vac/6
phi_sput = c_geo*Dz*radconst*nH2O # sputtered flux divided by the molecular number density

'''
fig1 = plt.figure(1)
#ax1 = fig.add_subplot(2,1,1)
#ax1.semilogy(T,Rconjkr*1E4, label = 'Rcon,jkr')
#ax1.semilogy(T,Rconminlin*4e3, label = 'Rcon,linear')
#ax1.semilogy(T,Rconminquad*4e3, label = 'Rcon,quad.')
#ax1.set_ylabel('Contact Radius [um]')
#ax1.axis([60,110,5e-3,1])
#ax1.legend(loc = 3, prop={'size':14})
ax1= fig1.add_subplot(1,1,1)
ax1.semilogy(T,Dsurf, label = 'Dsurf')
ax1.semilogy(T,DlatB, label = 'Dlat')
ax1.semilogy(T,Dgb, label = 'Dgb')
ax1.semilogy([60,110],[Drad_eff,Drad_eff], label = 'Dradsurf')
#ax1.semilogy(T,Dradlat, label = 'Dradlat')
#ax2.semilogy(T,Dradsput, label = 'Dradsput')
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('Diffusion Coeff. [cm^2/s]')
ax1.legend(loc = 4, prop={'size':14})
plt.savefig('./DiffusionRates.png',dpi=600)
'''

# ----- Calculate radiation induced diffusion coefficients ----
switches = [2]
for dzswitch in switches:
    if dzswitch == 1:
        lineT = '-'
        Rmin1=Rconminlin[0]
        Rmin2=Rconminquad[0]
        Rmax1=Rconmaxlin[0]
        Rmax2=Rconmaxquad[0]
        chi_min2 = chiminquad[0]
        chi_max2 = chimaxquad[0]
        Vmin2=Vminquad[0]
        Vmax2=Vmaxquad[0]
        Dsurfc=Dsurf[0]
        Dgbc=Dgb[0]
        DlatBc=DlatB[0]
        invTc=invT[0]
        Tc=T[0]
        Pevc=Pev[0]
        
    if dzswitch == 2:
        lineT = '--'
        Rmin1=Rconminlin[20]
        Rmin2=Rconminquad[20]
        Rmax1=Rconmaxlin[20]
        Rmax2=Rconmaxquad[20]
        chi_min2 = chiminquad[20]
        chi_max2 = chimaxquad[20]
        Vmin2=Vminquad[20]
        Vmax2=Vmaxquad[20]
        Dsurfc=Dsurf[20]
        Dgbc=Dgb[20]
        DlatBc=DlatB[20]
        invTc=invT[20]
        Tc=T[20]
        Pevc=Pev[20]
        
    if dzswitch == 3:
        lineT = '-.'
        Rmin1=Rconminlin[40]
        Rmin2=Rconminquad[40]
        Rmax1=Rconmaxlin[40]
        Rmax2=Rconmaxquad[40]
        chi_min2 = chiminquad[40]
        chi_max2 = chimaxquad[40]
        Vmin2=Vminquad[40]
        Vmax2=Vmaxquad[40]
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
    dVdtradsurf = 3*pi*Rcon_var*Drad_eff*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invTc/d2
    dVdtradlat = 3*pi*Rcon_var*Drad_eff*gamma_sv*Omega*(K3-Km)*invkb*invTc
    dVdtradsput = 2*pi*Rcon_var*theta*Rnc*phi_sput*gamma_sv*Omega*invkb*invTc*(K3-Km)
#    dVdtradsput = 2*pi*Rcon_var*theta*Dradsput*gamma_sv*Omega*invkb*invTc*(K3-Km)

    dRdtradsurf=dVdtradsurf/Vol_to_rad
    dRdtradlat=dVdtradlat/Vol_to_rad
    dRdtradsput=dVdtradsput/Vol_to_rad
    
    dVdtrad_tot = dVdtradsurf + dVdtradlat + dVdtradsput 
    dRdtrad_tot = dVdtrad_tot/Vol_to_rad

    index = find_nearest(Rcon_var,0.0001)
#    print radconst, Veff_lat[index],inv_tau_lat[index]
#    print Drad_eff, dRdtradlat[index], dRdtrad_tot[index], dRdttherm_tot[index]

    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(1,1,1)
    ax2.semilogy(Rcon_var/rg,dRdttherm_tot, linewidth=2, color='r', label = 'Therm. total', linestyle = lineT)
    ax2.semilogy(Rcon_var/rg,dRdtrad_tot, linewidth=2, color='b', label = 'Rad. total', linestyle = lineT)

#    ax3 = fig.add_subplot(2,1,1)
#    ax2.semilogy([Rmin1/rg,Rmin1/rg],[1e-28,1e-10], color='k', linewidth=2, linestyle = lineT) #label= 'T = %d' % Tc,
#    ax2.semilogy([Rmax1/rg,Rmax1/rg],[1e-28,1e-10], color='k', linewidth=2, linestyle = lineT)
    ax2.semilogy([Rmin2/rg,Rmin2/rg],[1e-28,1e-10], color='b', linewidth=2, linestyle = lineT)
    ax2.semilogy([Rmax2/rg,Rmax2/rg],[1e-28,1e-10], color='b', linewidth=2, linestyle = lineT) 
#    ax3.legend(prop={'size':18},loc=3, bbox_to_anchor=(1,-3.5))

    if dzswitch == 2:
#        ax2.set_xlim([0,0.2])
        ax2.semilogy(Rcon_var/rg,dRdtradsurf, color='m', label = 'Rad. surface', linestyle = lineT)
        ax2.semilogy(Rcon_var/rg,dRdtradlat, color = 'g', label = 'Rad. lattice', linestyle = lineT)
        ax2.semilogy(Rcon_var/rg,dRdtradsput, color='c', label = 'Rad. sputter', linestyle = lineT)
        box = ax2.get_position()
#        ax2.set_position([box.x0,box.y0,box.width,box.height*1.8])
        ax2.set_xlabel('Contact Ratio R_{con}/r_g')
        ax2.set_ylabel('Radial Sintering Rate [cm/s]')
        ax2.legend(loc = 1, prop={'size':18})
        
#        box = ax3.get_position()
#       ax3.set_position([box.x0,box.y0+0.275,box.width,box.height*0.2])
#        plt.setp(ax3.get_yticklabels(), visible=False)
#        ax3.xaxis.tick_top()
#        ax3.xaxis.set_label_position('top')
#        ax3.set_xlabel('Contact Radius [um]')

    
#    if dzswitch == 2:
#        plt.savefig('./Rad_Sint_Rates.png',dpi=600)
    
#Rmins = [1e-6, 1e-5, 5e-5] # beginning radius in cm
#Rmax = 10e-4 # final radius in cm (10um) 
#for Rmin in Rmins:     
    Rcur = Rmin2
    chi_cur = chi_min2
    Vcur = Vmin2
    ttot=0
    check = 0
    V_t = []
    timesteps = []
    
#    print 'Rcon, dRdtrad_tot, dRdttherm_tot'
#    for i in range(1000):
    print '%.3e' % Vcur, '%.3e' % dVdtrad_tot[index], '%.3e' % dVdttherm_tot[index]
    while Vcur<Vmax2:
#        if Rcur<1e-5:
#            dt = 3.157e5
#        else:
        dt = 3.157e11

        index = find_nearest(V_var,Vcur)
        DV = dVdtrad_tot[index]+dVdttherm_tot[index]
        Vcur+= N*DV*dt
        V_t.append(Vcur)
        ttot+=dt
        if Vcur > Vmax2 and check == 0:
            ttot_years = ttot*stoyear
            check = 1
        timesteps.append(ttot*stoyear)
    print '%.3e' % Vcur, '%.3e' % dVdtrad_tot[index], '%.3e' % dVdttherm_tot[index]
    
    V_tarr=np.array(V_t,dtype='float')

    print 'the total time for %d K' % Tc, ' is %.4e' % ttot_years, ' years' 

#    print len(Rcon_var), len(dRdtradlat), len(dRdtrad_tot)

    fig = plt.figure(3)
    ax = fig.add_subplot(1,1,1)
    ax.plot(timesteps,V_tarr,linewidth=2, label = 'chi_min = {:.1e}'.format(chi_min2), linestyle = lineT)
#    line, = ax.plot(R_tarr*1e4,timesteps, linewidth=2, color='r', label = 'Rad. sint @ %d' %Tc, linestyle = lineT)
#    ax.annotate('Rmax @ %d K' % Tc,xy=(ttot,Rcur*1e4),xytext=(ttot,Rcur+2.4),arrowprops=dict(facecolor='black',shrink=0.05), horizontalalignment='center')
    if dzswitch == 2:
        ax.set_xlabel('Contact Ratio R_{con}/r_g')
        ax.set_ylabel('Time [years]')
    if dzswitch == 2:
        ax.legend(loc = 2, prop={'size':18})
#        plt.savefig('./sint_timescales.png',dpi=600)

plt.show()

