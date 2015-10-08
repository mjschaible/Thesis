'''Program to calculate the volumetric/radial sintering rate as a function of time in order to determine the time scale at which the radiation induced sintering can increase the contact radius between grains from some initial (Hertzian, or calculated from TC outside the anomalous region) to the contact radius inside the anomalous region calculated from the Wood theory and the measured thermal inertia of the regolith'''

# dt determines the time step for the simulation.
# for log scaling (num=1e4), Tlin = 4.4E10s, Tquad=1.58E14s, Thertz=1.59E14s
#dt = np.logspace(3,15,num=1e2)
# Mimas (lin num=1e5): Tlin = 2.7e9s, Tquad = 1.2e13s, Thertz = 1.2e13s
# Tethys: Tlin = 1.9E7s, Tquad = 4.32e11s, Thertz = s
# Dione: Tlin = 188s, Tquad = 1.49E12s , Thertz = 1.49E12s


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

rg = 25E-4 # cm, approximate grain radius (25 um) (all bodies)
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
inv_tau_vac = Nvac*radconst
Drad_eff = lat_const*lat_const*inv_tau_vac/6

# Calculate the approximate length scale of an electron excitation event 
Ei = 5 # eV, heating energy to lattice
U = 0.6 # eV/atom, characteristic energy of the solid
Dz = 0.25*delta_surf*(2*Ei/U -1) # cm, length-scale of excitation event assuming binary collision
phi_sput = c_geo*Dz*radconst # sputtered flux divided by the molecular number density


# ---- Functions ----
def find_nearest(Rcon,Rcur):
    index = (np.abs(Rcon-Rcur)).argmin()
    return index

# ----------- Begin main program -------------------------------------------

T = np.linspace(60,110,num=51) # Define temperature range for calculations
invT = 1/T # Calculate the inverse of temperature
Tnorm = T/273 # Calculate the normalized temprature (normalize to melting temp)

# ----- Determine the variation in the sintering rates based on the contact radius -----
Rcon_var = np.logspace(-6, -3, num=1000, base=10.0) # in cm
invRcon_var = 1/Rcon_var

Rnc = Rcon_var**2/(2*(rg-Rcon_var))
invRnc = 1/Rnc
theta = np.arctan(rg/(Rcon_var-Rnc))
Km = invRcon_var - invRnc
K3 = 2/rg

# Calculate Nc = Neighbor contacts, depends on the porosity from curve fitting  [Yang et al, 2000]
Nc = 2.02*((1+87.38*(1-phi)**4)/(1+25.81*(1-phi)**4))

# Calculate the cementation volume fraction
chi_var = (3/16)*Nc*(Rcon_var/rg)**4

# ----- Determine the thermal conductivity of ice -----
#kcem = [0.01, 0.1, 1.0, 2.0, 6.5] # [J/m/s/K], cementation TC
#kcem = map(lambda x: x*1e5, kcem)  # convert to [erg/cm/s/K]
kxice = (488.19/T+0.4685)*1e5 # [erg/cm/s/K], TC of Ih (hexagonal) crystalline water ice
kaice = 7.1e2*T # [erg/cm/s/K], TC of amorphous water ice

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
fsc_lin = 1/Apflin
Apfquad = (1/(Ysc*Nc))*(2*np.sqrt(Nc-1)*rg/Nc)**(Zsc2)
fsc_quad = 1/Apfquad
# ----- Calculate linear and quadratic minimum contact radius from the Wood model -----
Rconminlin = (Apflin*kout*(1+0.5*phi)/(kxice*(1-phi)))
Rconmaxlin = (Apflin*kin*(1+0.5*phi)/(kxice*(1-phi)))
Rconminquad = np.sqrt(Apfquad*kout*(1+0.5*phi)/(kxice*(1-phi)))
Rconmaxquad = np.sqrt(Apfquad*kin*(1+0.5*phi)/(kxice*(1-phi)))
RconminSY = np.sqrt((kout/kxice)*SYpf)
RconmaxSY = np.sqrt((kin/kxice)*SYpf)
chiminlin = (Rconminlin/(rg*0.4739))**4
chiminquad = (Rconminquad/(rg*.4739))**4


switches = [1,2,3]
for dzswitch in switches:
    if dzswitch == 1:
        lineT = '-'
        kg=kxice[20]
        kcem=kxice[20]
        
    if dzswitch == 2:
        lineT = '--'
        kg=kxice[20]
        kcem=kaice[20]
        
    if dzswitch == 3:
        lineT = '-.'
        kg=kaice[20]
        kcem=kaice[20]

    # ----- Calculate the effective TC from the Wood model -----
    structure_factor = (kcem*chi_var + kg*(1-phi)*3*kcem/(2*kcem+kg))/(chi_var+(1-phi)*3*kcem/(2*kcem+kg)+1.5*phi)
    keff_lin = fsc_lin*structure_factor*Rcon_var**Zsc1
    keff_quad = fsc_quad*structure_factor*Rcon_var**Zsc2
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(1,1,1)
    ax1.loglog(Rcon_var*1e4,keff_lin/1e5, label = 'kg={:.2f},kcem={:.2f}'.format(kg/1e5, kcem/1e5), linestyle = lineT, color = 'g')
    ax1.loglog(Rcon_var*1e4,keff_quad/1e5,linestyle = lineT, color = 'b')
    if dzswitch == 3:
        ax1.loglog([1e-2,10],[kout/1e5,kout/1e5], linewidth=2, color='r')
        ax1.loglog([1e-2,10],[kin/1e5,kin/1e5], linewidth=2, color='k')
        ax1.set_xlabel('Contact Radius [um]')
        ax1.set_ylabel('Effective Thermal Conductivity [W/mK]')
        ax1.set_ylim([1e-5,1])
    ax1.legend(loc = 2, prop={'size':14})
plt.savefig('./TCeff.png',dpi=600)
plt.show()

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
#    Veff_surf = 4*pi*delta_surf*Rcon_var*Dz
#    Veff_lat = 2*pi*Rcon_var*Rcon_var*Dz
#    Veff_sput = Dz*rg*rg*(0.25-0.125*(np.sqrt(rg*rg-Rcon_var*Rcon_var)+rg)/rg)

#    inv_tau_surf = radconst*Veff_surf
#    inv_tau_lat = radconst*Veff_lat
#    inv_tau_sput = radconst*Veff_sput
    
#    Dradsurf = Dz*Dz*inv_tau_surf
#    Dradlat = Dz*Dz*inv_tau_lat
#    Dradsput = Dz*Dz*inv_tau_sput

    dVdtradsurf = 3*pi*Rcon_var*Drad_eff*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invTc/d2
    dVdtradlat = 3*pi*Rcon_var*Drad_eff*gamma_sv*Omega*(K3-Km)*invkb*invTc
    dVdtradsput = 2*pi*Rcon_var*theta*Rnc*phi_sput*gamma_sv*Omega*invkb*invTc*(K3-Km)
#    dVdtradsput = 2*pi*Rcon_var*Dradsput*gamma_sv*Omega*invkb*invTc*(K3-Km)

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
    ax2.semilogy(Rcon_var*4e3,dRdttherm_tot, linewidth=2, color='r', label = 'Therm. total', linestyle = lineT)
    ax2.semilogy(Rcon_var*4e3,dRdtrad_tot, linewidth=2, color='b', label = 'Rad. total', linestyle = lineT)

    #   ax3 = fig.add_subplot(2,1,1)
    ax2.semilogy([Rmin1*4e3,Rmin1*4e3],[1e-28,1e-10], color='k', linewidth=2, linestyle = lineT) #label= 'T = %d' % Tc,
    ax2.semilogy([Rmax1*4e3,Rmax1*4e3],[1e-28,1e-10], color='k', linewidth=2, linestyle = lineT)
 #   ax2.semilogy([Rmin2*4e3,Rmin2*4e3],[1e-28,1e-10], color='b', linewidth=2, linestyle = lineT)
 #   ax2.semilogy([Rmax2*4e3,Rmax2*4e3],[1e-28,1e-10], color='b', linewidth=2, linestyle = lineT) 
 #   ax3.legend(prop={'size':18},loc=3, bbox_to_anchor=(1,-3.5))

    if dzswitch == 2:
        ax2.set_xlim([0,0.2])
        ax2.semilogy(Rcon_var*4e3,dRdtradsurf, color='m', label = 'Rad. surface', linestyle = lineT)
        ax2.semilogy(Rcon_var*4e3,dRdtradlat, color = 'g', label = 'Rad. lattice', linestyle = lineT)
        ax2.semilogy(Rcon_var*4e3,dRdtradsput, color='c', label = 'Rad. sputter', linestyle = lineT)

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

    
    if dzswitch == 2:
        plt.savefig('./Rad_Sint_Rates.png',dpi=600)
    
#Rmins = [1e-6, 1e-5, 5e-5] # beginning radius in cm
#Rmax = 10e-4 # final radius in cm (10um) 
#for Rmin in Rmins:     
    Rcur = Rmin1
#    print 'The min rad is %.2e' % Rmin, ' and the max is %.2e' % Rmax
    ttot=0
    check = 0
    R_t = []
    timesteps = []

#    print 'Rcon, dRdtrad_tot, dRdttherm_tot'
#    for i in range(10):
    while Rcur<Rmax1+1e-5:
        if Rcur<1e-5:
            dt = 3.157e3
        else:
            dt = 3.157e7

        index = find_nearest(Rcon_var,Rcur)
        DR = dRdtrad_tot[index]+dRdttherm_tot[index]
        Rcur+=DR*dt
        R_t.append(Rcur)
        ttot+=dt
        if Rcur > Rmax1 and check == 0:
            ttot_years = ttot*stoyear
            check = 1
        timesteps.append(ttot*stoyear)
#        print '%.3e' % Rcur, '%.3e' % dRdtrad_tot[index], '%.3e' % dRdttherm_tot[index]
    
    R_tarr=np.array(R_t,dtype='float')

    print 'the total time for %d K' % Tc, ' is %.4e' % ttot_years, ' years' 

#    print len(Rcon_var), len(dRdtradlat), len(dRdtrad_tot)

    fig = plt.figure(3)
    ax = fig.add_subplot(1,1,1)
    ax.plot(timesteps,R_tarr*4e3,linewidth=2, label = 'Rmin = %d' %Rmin1, linestyle = lineT)
#    line, = ax.plot(R_tarr*1e4,timesteps, linewidth=2, color='r', label = 'Rad. sint @ %d' %Tc, linestyle = lineT)
#    ax.annotate('Rmax @ %d K' % Tc,xy=(ttot,Rcur*1e4),xytext=(ttot,Rcur+2.4),arrowprops=dict(facecolor='black',shrink=0.05), horizontalalignment='center')
    if dzswitch == 2:
        ax.set_xlabel('Contact Ratio R_{con}/r_g')
        ax.set_ylabel('Time [years]')
    if dzswitch == 2:
        ax.legend(loc = 2, prop={'size':18})
        plt.savefig('./sint_timescales.png',dpi=600)
'''

