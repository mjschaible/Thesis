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


def find_nearest(Rcon,Rcur):
    index = (np.abs(Rcon-Rcur)).argmin()
    return index

# --------------- Define necessary model parameters -------------------------
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
nH2O = 3e22 # molec/cm^3, molecular number density
alat = 3e-8 # cm, approximate lattice spacing for H2O
rho_ice = 0.934 # g/cm^3, density of water ice @ 80K
WH2O = (27)*1.602e-12 # (eV)erg, average energy deposited in an excitation event
UH2O = (0.7)*1.602e-12 # (eV)erg, bulk cohesive energy of water
DeltaE = 0.2*WH2O

c_geo = 1 # geometric factor to determine sputtering/deposition relation
beta = 6/(pi*pi) # factor for calculating radiation stuff

# Calculate the radiation time constant
invtau = dEdz*Phi_elec/(WH2O*nH2O)

# ------ Define Wood model structural parameters -----
Ysc = 0.22 # Obtained from curve fitting experimental data
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

# ----------- Begin main program -------------------------------------------

T = np.linspace(60,110,num=51) # Define temperature range for calculations
invT = 1/T # Calculate the inverse of temperature
Tnorm = T/273 # Calculate the normalized temprature (normalize to melting temp)

# ----- Calculate the equilibrium vapor pressure of ice -----
Pev = 3.65e13 * np.exp(-6141.7 * invT)

# ----- Determine the thermal conductivity of ice -----
#kcem = [0.01, 0.1, 1.0, 2.0, 6.5] # [J/m/s/K], cementation TC
#kcem = map(lambda x: x*1e5, kcem)  # convert to [erg/cm/s/K]
kxice = (488.19/T+0.4685)*1e5 # [erg/cm/s/K], TC of Ih (hexagonal) crystalline water ice
kaice = 7.1e2*T # [erg/cm/s/K], TC of amorphous water ice

# ----- Calculate thermal diffusion coefficients for ice -----
invR = 1/8.314 # [R] = [J/mol K], gas constant
DsurfIc = 1.74e5 * np.exp(-38.2e3*invR*invT) # cm^2/s, from Kouchi 1994
Dsurf = 3e2 * np.exp(-44.13e3*invR*invT) # cm^2/s, from Maeno and Ebinume, 1983
#Dsurf = 1.4e-4*np.exp(-23.2e3*invR*invT)  # cm^2/s, from Nasello et al 1994
DlatLDA = 7e-6*np.exp(-15e3*invR*invT) # from Ghesquiere 2015
DlatIh = 1e-10*np.exp(-9e3*invR*invT) # from Ghesquiere 2015
DlatB = 4.2e8*np.exp(-71.1e3*invR*invT) # from Livingston 1998
Dlat = 1.5e-11*np.exp(-59.8e3*invR*invT)  # cm^2/s, from Nasello et al 1994
Dvac = 1e-3*np.exp(-0.62*invT/8.62e-5) # from Hobbs, 1974
Dgb = 8.4*np.exp(-49e3*invR*invT) # from ....

# ----- Calculate radiation induced diffusion coefficients ----
# Calculate the approximate length scale of an electron excitation event 
l_c = 0.25*alat*(beta*DeltaE/UH2O) # cm, length-scale of excitation event
# or
Dz_min = 2*alat
Dz_max = 5*alat
R_iv = alat/2 # GUESS, separation distance for spontaneous vacancy/interstitial recombination
alpha_v = 0.1 # GUESS, average number of point defects which escape cascade
d_vac = 5e-6 # cm, approximate average sink diffusion length scale

Ea_bulk_min = (0.6)*1.602e-12 # (eV) erg, min bulk activation energy
Ea_bulk_max = (0.75)*1.602e-12 # (eV) erg, max bulk activation energy 

Ea_surf_min = (0.2)*1.602e-12 # (eV) erg, min surface diffusion activation energy
Ea_surf_max = (0.5)*1.602e-12 # (eV) erg, max surface diffusion activation energy

Nmob_surf_max = beta * DeltaE / (Ea_surf_max) # number of surface molecules mobilized per ionization
Nmob_surf_min = beta * DeltaE / (Ea_surf_min) # number of surface molecules mobilized per ionization

Nmob_bulk_max = beta * DeltaE / (Ea_bulk_max) # number of surface molecules mobilized per ionization
Nmob_bulk_min = beta * DeltaE / (Ea_bulk_min) # number of surface molecules mobilized per ionization

N_FP = 0.8*DeltaE/(2*UH2O)

Deff_mc_surf_min = Nmob_surf_min*Dz_min*Dz_min *invtau/4
Deff_mc_surf_max = Nmob_surf_max*Dz_max*Dz_max *invtau/4
Deff_mc_bulk_min = Nmob_bulk_min*Dz_min*Dz_min *invtau/6
Deff_mc_bulk_max = Nmob_bulk_max*Dz_max*Dz_max *invtau/6

Deff_esd_mid = 2*alpha_v*d_vac*d_vac*N_FP*invtau 
Deff_esd_low = (alpha_v*N_FP*Dvac*invtau/(pi*nH2O*R_iv))**0.5

Deff_surf_min = Deff_mc_surf_min + Deff_esd_mid + Deff_esd_low
Deff_bulk_min = Deff_mc_bulk_min + Deff_esd_mid + Deff_esd_low
Deff_surf_max = Deff_mc_surf_max + Deff_esd_mid + Deff_esd_low
Deff_bulk_max = Deff_mc_bulk_max + Deff_esd_mid + Deff_esd_low

phi_sput = c_geo*l_c*invtau*nH2O # sputtered flux

# ----- Calculate linear and quadratic minimum contact radius from the Wood model -----
Rconminlin = (Apflin*kout*(1+0.5*phi)/(kxice*(1-phi)))
Rconmaxlin = (Apflin*kin*(1+0.5*phi)/(kxice*(1-phi)))
Rconminquad = np.sqrt(Apfquad*kout*(1+0.5*phi)/(kxice*(1-phi)))
Rconmaxquad = np.sqrt(Apfquad*kin*(1+0.5*phi)/(kxice*(1-phi)))
RconminSY = np.sqrt((kout/kxice)*SYpf)
RconmaxSY = np.sqrt((kin/kxice)*SYpf)
'''
print "Rconminlin = %.2e" % (Rconminlin[20]*1E4), " um"
print "Rconmaxlin = %.2e" % (Rconmaxlin[20]*1e4), " um"
print "Rconminquad = %.2e" % (Rconminquad[20]*1E4), " um"
print "Rconmaxquad = %.2e" % (Rconmaxquad[20]*1e4), " um"
print "RconSYmax = %.2e" % (RconmaxSY[20]*1e4), " um"
print "RconSYmin = %.2e" % (RconminSY[20]*1e4), " um"

fig = plt.figure(1)
#ax1 = fig.add_subplot(2,1,1)
#ax1.semilogy(T,Rconjkr*1E4, label = 'Rcon,jkr')
#ax1.semilogy(T,Rconminlin*1e4, label = 'Rcon,linear')
#ax1.semilogy(T,Rconminquad*1e4, label = 'Rcon,quad.')
#ax1.set_ylabel('Contact Radius [um]')
#ax1.axis([60,110,5e-3,1])
#ax1.legend(loc = 3, prop={'size':14})

ax2= fig.add_subplot(1,1,1)
#ax2.semilogy(T,DsurfIc, label = 'DsurfIc')
ax2.semilogy(T,Dsurf, label = 'Dsurf')
ax2.semilogy(T,DlatB, label = 'Dlat')
ax2.semilogy(T,Dgb, label = 'Dgb')
ax2.semilogy(T,Dvac, label = 'Dvac')
ax2.semilogy([60,110],[Deff_mc_bulk,Deff_mc_bulk], label = 'Drad_mc')
ax2.semilogy([60,110],[Deff_esd_mid,Deff_esd_mid], label = 'Drad_esd, mid')
ax2.semilogy(T,Deff_esd_low, label = 'Drad_esd, low')
ax2.set_xlabel('Temperature [K]')
ax2.set_ylabel('Diffusion Coeff. [cm^2/s]')
ax2.legend(loc = 4, prop={'size':14})

#plt.savefig('./DiffusionRates.png',dpi=600)
'''

# ----- Determine the variation in the sintering rates based on the contact radius -----
Rcon_var = np.logspace(-6, -3, num=100, base=10.0) # in cm
invRcon_var = 1/Rcon_var

Rnc = Rcon_var**2/(2*(rg-Rcon_var))
invRnc = 1/Rnc
theta = np.arctan(rg/(Rcon_var-Rnc))
Km = invRcon_var - invRnc
K3 = 2/rg

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
        Dlatc=Dlat[0]
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
        Dlatc=Dlat[20]
        DlatIhc=DlatIh[20]
        Dvacc=Dvac[20]
        Deff_bulk_minc = Deff_bulk_min[20]
        Deff_surf_minc = Deff_surf_min[20]
        Deff_bulk_maxc = Deff_bulk_max[20]
        Deff_surf_maxc = Deff_surf_max[20]
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
        Dlatc=Dlat[40]
        DlatIhc=DlatIh[40]
        Dvacc=Dvac[40]
        Deff_bulkc = Deff_bulk[40]
        Deff_surfc = Deff_surf[40]
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
    dVdtlat = 3*pi*Rcon_var*Dlatc*gamma_sv*Omega*(K3-Km)*invkb*invTc
    vapDiff = Pevc*Rnc*theta*np.sqrt(Omega*invkb*invTc/(2*pi*rho_ice))
    dVdtvap = 2*pi*Rcon_var*vapDiff*gamma_sv*Omega*invkb*invTc*(K3-Km)

    dVdttherm_tot = dVdtsurf + dVdtlat + dVdtvap
    dRdttherm_tot = dVdttherm_tot/Vol_to_rad
    
    dVdtrad_tot = []
    Surf = np.linspace(Deff_surf_minc, Deff_surf_maxc, num=20)
    Bulk = np.linspace(Deff_bulk_minc, Deff_bulk_maxc, num=20)
    for Deff_surf in Surf:
        for Deff_bulk in Bulk:
            dVdtradsurf = 3*pi*Rcon_var*Deff_surf*delta_surf*gamma_sv*Omega*(K3-K2)*invkb*invTc/d2
            dVdtradlat = 3*pi*Rcon_var*Deff_bulk*gamma_sv*Omega*(K3-Km)*invkb*invTc
            dVdtradsput = 2*pi*Rcon_var*theta*Rnc*(phi_sput/nH2O)*gamma_sv*Omega*invkb*invTc*(K3-Km)
            dVdtrad_tot.append(dVdtradsurf + dVdtradlat + dVdtradsput)

            dRdtradsurf=dVdtradsurf/Vol_to_rad
            dRdtradlat=dVdtradlat/Vol_to_rad
            dRdtradsput=dVdtradsput/Vol_to_rad
            dRdtrad_tot = dVdtrad_tot/Vol_to_rad

'''
            fig = plt.figure(2)
            ax2 = fig.add_subplot(1,1,1)
            ax2.semilogy(Rcon_var/rg,dVdttherm_tot, linewidth=2, color='r', label = 'Therm. total')
            ax2.semilogy(Rcon_var/rg,dVdtrad_tot, linewidth=2, color='b', label = 'Rad. total', linestyle = lineT)

            ax2.semilogy([Rmin2/rg,Rmin2/rg],[1e-34,1e-24], color='k', linewidth=2, linestyle = '-.')
            ax2.semilogy([Rmax2/rg,Rmax2/rg],[1e-34,1e-24], color='k', linewidth=2, linestyle = '-.')

            if dzswitch == 2:
                ax2.semilogy(Rcon_var/rg,dVdtradsurf, color='m', label = 'Rad. surface', linestyle = lineT)
                ax2.semilogy(Rcon_var/rg,dVdtradlat, color = 'g', label = 'Rad. lattice', linestyle = lineT)
                ax2.semilogy(Rcon_var/rg,dVdtradsput, color='c', label = 'Rad. sputter', linestyle = lineT)

        #        ax2.loglog(Rcon_var/rg,dVdtsurf, color='m', label = 'Therm. surface')
        #        ax2.loglog(Rcon_var/rg,dVdtlat, color = 'g', label = 'Therm. lattice')
        #        ax2.loglog(Rcon_var/rg,dVdtvap, color='c', label = 'Therm. vapor')

                box = ax2.get_position()
                ax2.set_xlabel('Contact/Grain radius (Rcon/rg)')
                ax2.set_ylabel('Volumetric Sintering Rate [cm^3/s]')
                ax2.set_xlim([0,0.1])
                ax2.legend(loc = 3, prop={'size':12})

        #    plt.savefig('./Rad_Sint_Rates.png',dpi=600)


            Rcur = Rmin2
         #   print 'The min rad is %.2e' % Rmin2, ' and the max is %.2e' % Rmax2
            ttot=0
            R_t = []
            timesteps = []
            dt = 3.157e9
            print 'Rcon, dRdtrad_tot, dRdttherm_tot'
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

            fig = plt.figure(3)
            ax = fig.add_subplot(1,1,1)
            line, = ax.plot(timesteps,R_tarr*1e4, linewidth=2, color='r', label = 'Rad. sint @ %d' %Tc, linestyle = lineT)
            ax.annotate('Rmax @ %d K' % Tc,xy=(ttot,Rcur*1e4),xytext=(ttot,Rcur+2.4),arrowprops=dict(facecolor='black',shrink=0.05), horizontalalignment='center')
            if dzswitch == 1:
                ax.set_ylabel('Contact Radius [um]')
                ax.set_xlabel('Timestep [s]')
            if dzswitch == 3:
                ax.legend(loc = 3, prop={'size':18})
        #        plt.savefig('./sint_timescales.png',dpi=600)

        #plt.show()
'''
