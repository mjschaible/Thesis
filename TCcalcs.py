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
import matplotlib.cm as cm
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

# Calculate the cementation volume fraction
def chi_calc(X, Nc):
    chi = 0.1875*Nc*(X)**4
    return chi

def tick_function(X):
    chi = 0.1875*5*(X)**4
    return ["%.1e" % z for z in chi]


# ----------- Begin main program -------------------------------------------
T = np.linspace(60,110,num=51) # Define temperature range for calculations
invT = 1/T # Calculate the inverse of temperature
Tnorm = T/273 # Calculate the normalized temprature (normalize to melting temp)

# ----- Determine the variation in the sintering rates based on the contact radius -----
Rcon_var = np.logspace(-7, -3, num=1000, base=10.0) # in cm
invRcon_var = 1/Rcon_var

Rrat = Rcon_var/rg

Rnc = Rcon_var**2/(2*(rg-Rcon_var))
invRnc = 1/Rnc
theta = np.arctan(rg/(Rcon_var-Rnc))
Km = invRcon_var - invRnc
K3 = 2/rg

# ----- Determine the thermal conductivity of ice -----
#kcem = [0.01, 0.1, 1.0, 2.0, 6.5] # [J/m/s/K], cementation TC
#kcem = map(lambda x: x*1e5, kcem)  # convert to [erg/cm/s/K]
kxice = (488.19/T+0.4685)*1e5 # [erg/cm/s/K], TC of Ih (hexagonal) crystalline water ice
kaice = 7.1e2*T # [erg/cm/s/K], TC of amorphous water ice

phi_arr = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95] # Assume 50% porosity for now (all bodies)

colors = iter(cm.rainbow(np.linspace(0, 1, len(phi_arr))))

for phi in phi_arr:
    # ------ Define Sirono and Yamamoto model parameters ------
    g_SY = 4
    pc_SY = 0.333
    p_SY = 6*(1-phi)/pi
    SYpf = ((1-pc_SY)/(p_SY-pc_SY))*g_SY*rg*rg/pi

    # ------ Define Wood model structural parameters -----
    # Calculate Nc = Neighbor contacts, depends on the porosity from curve fitting  [Yang et al, 2000]
    Nc = 2.02*((1+87.38*(1-phi)**4)/(1+25.81*(1-phi)**4))
    chi_var = chi_calc(Rrat, Nc)

    Ysc = 0.22 # Obtained from curve fitting experimental data
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

    switches = [1]
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

        Rmin1=Rconminlin[20]
        Rmin2=Rconminquad[20]

        # ----- Calculate the effective TC from the Wood model -----
        structure_factor = (kcem*chi_var + kg*(1-phi)*3*kcem/(2*kcem+kg))/(chi_var+(1-phi)*3*kcem/(2*kcem+kg)+1.5*phi)
        keff_lin = fsc_lin*structure_factor*Rcon_var**Zsc1
        keff_quad = fsc_quad*structure_factor*Rcon_var**Zsc2

        keff_lin_nz = fsc_lin*structure_factor*(Rmin1 + Rcon_var)**Zsc1
        keff_quad_nz = fsc_quad*structure_factor*(Rmin2 + Rcon_var)**Zsc2

        keff_SY = kg*pi/g_SY * (p_SY - pc_SY)/(1 - pc_SY) * (Rcon_var/rg)**Zsc2

        fig = plt.figure(1)
        ax1 = fig.add_subplot(1,1,1)
#        ax2 = ax1.twiny()

        curcol = next(colors)

    #    ax1.loglog(Rcon_var/rg,keff_lin/kg, label = 'kg={:.2f},kcem={:.2f}'.format(kg/1e5, kcem/1e5), color = 'g')
        ax1.loglog(Rcon_var/rg,keff_lin/kg,linestyle = '--', color = curcol)
        ax1.loglog(Rcon_var/rg,keff_quad/kg,linestyle = '-', color = curcol, label = 'phi={:.2f}'.format(phi))
#        ax1.loglog(Rcon_var/rg,keff_lin_nz/kg,linestyle = '--', color = curcol)
#        ax1.loglog(Rcon_var/rg,keff_quad_nz/kg,linestyle = '-', color = curcol, label = 'phi={:.2f}'.format(phi) )
#        ax1.loglog(Rcon_var/rg,keff_SY/kg,label = 'S&Y, keff~Rcon^2', color = 'c')

        xmin = 5e-5
        xmax = 0.5
        ax1.loglog([1e-5,.50],[kout/kg,kout/kg], linewidth=2, color='r')
        ax1.loglog([1e-5,.50],[kin/kg,kin/kg], linewidth=2, color='k')
        ax1.set_xlabel('Contact to grain radius ratio, Rcon/rg')
        ax1.set_xlim([xmin,xmax])
        ax1.text(6e-5,2e-3,'Mimas, In')
        ax1.text(1.75e-2,8e-5,'Mimas, Out')
        ax1.set_ylabel('keff/kice')
        ax1.set_ylim([5e-5,.01])
        ax1.legend(loc = 2, prop={'size':14})

  #      ax1Ticks = ax1.get_xticks()
  #      ax2Ticks = ax1Ticks
  #      new_tick_locations=np.logspace(-2,5,num=7)

  #      ax2.set_xticks(new_tick_locations)
  #      ax2.set_xbound(ax1.get_xbound())
  #      ax2.set_xticklabels(tick_function(new_tick_locations))
        
  #           ax2.loglog(chi_var, keff_quad/kg, alpha=0)
  #      ax2.set_xlabel('Cementation Volume Fraction, Vcem/Vg')

plt.savefig('./TCeff.png',dpi=600)
#plt.show()
