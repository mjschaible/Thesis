'''Program to calculate the volumetric/radial sintering rate as a function of time in order to determine the time scale at which the radiation induced sintering can increase the contact radius between grains from some initial (Hertzian, or calculated from TC outside the anomalous region) to the contact radius inside the anomalous region calculated from the Wood theory and the measured thermal inertia of the regolith'''

from decimal import Decimal
from decimal import getcontext
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from sintfuncs import fix_coeff
from sintfuncs import plot_JKR
from sintfuncs import plot_Diff
from sintfuncs import plot_VSR

def find_nearest(Rcon,Rcur):
    index = (np.abs(Rcon-Rcur)).argmin()
    return index

def tick_function(X):
    Tax = 1/X
    return Tax

def calc_TC(TI):
    TC = (TI*TI)/(c_H2O*rho_H2O_I*(1-phi))
    return TC

def calc_Drad(Ea, Dz, alpha_v, J):
    Nmob = beta * DeltaE / (Ea) # number of surface molecules mobilized per ionization

    Deff_mc = Nmob*Dz*Dz*invtau/J

    d = 5*Dz # approximate average sink diffusion length scale
    Deff_esd_mid = 2*alpha_v*d*d*N_FP*invtau 
    Deff_esd_low = (alpha_v*N_FP*Dvac*invtau/(pi*nH2O*R_iv))**0.5

    if J==4:
        Deff = Deff_mc# + Deff_esd_mid_min + Deff_esd_low_min
    elif J==6:
        Deff = Deff_mc + Deff_esd_mid + Deff_esd_low
    else:
        'Please select either Bulk or Surface diffusion'

    return Deff

# ----------- Begin main program -------------------------------------------
fignum = 0
T = np.linspace(60,110,num=51) # Define temperature range for calculations
invT = 1/T # Calculate the inverse of temperature
Tnorm = T/273 # Calculate the normalized temprature (normalize to melting temp)

# --------------- Define necessary model parameters -------------------------
pi = np.pi
kb = 1.38e-16 # erg/K, Boltzman constant
stoyear = 3.16888e-8 # multiply seconds by this to get years
invkb = 1/kb

# For now assume an average energy loss for all depths
# Should eventually use Penelope results to give vs. depth
dEdz = (1e6)*1.602e-12 # (eV)erg/cm, energy deposition rate at Mimas

# Define water ice parameters
gamma_ss = 65 # erg/cm^2, ice/ice (grain boundary) surface energy
gamma_sv = 109 # erg/cm^2, ice/vapor (free surface) surface energy
delta_gb = 9.04e-8 # cm, grain boundary thickness
delta_surf = 3.19e-8 # cm, effective surface (atomic layer) thickness (cube root of Omega)
Omega = 3.25e-23 # cm^3, molar volume
nH2O = 3e22 # molec/cm^3, molecular number density
alat = 3e-8 # cm, approximate lattice spacing for H2O
rho_ice = 0.934 # g/cm^3, ~density of water ice @ 80K
WH2O = (27)*1.602e-12 # (eV)erg, average energy deposited in an excitation event
UH2O = (0.7)*1.602e-12 # (eV)erg, bulk cohesive energy of water
DeltaE = 0.2*WH2O

N_vacH2O = nH2O*np.exp(-UH2O*invkb*invT) # Equilibrium vacancy concentration

# ----- Calculate the temp dependent specific heat capacity of ice -----
c_H2O = (0.09+0.00749*T)*1E7  # [erg/g/K], (Gutierez et al, 2001)

# ----- Calculate thermal diffusion coefficients for ice -----
invR = 1/8.314 # [R] = [J/mol K], gas constant
DsurfIc = 1.74e5 * np.exp(-38.2e3*invR*invT) # cm^2/s, from Kouchi 1994
#Dsurf = 3e2 * np.exp(-44.13e3*invR*invT) # cm^2/s, from Maeno and Ebinume, 1983
Dsurf = 1.4e-4*np.exp(-48.2e3*invR*invT) # cm^2/s, max estimate from Nie et al 2009
Dsurf_max = 1.4e-4*np.exp(-23.2e3*invR*invT)  # cm^2/s, from Nasello et al 2007
DlatLDA = 7e-6*np.exp(-15e3*invR*invT) # from Ghesquiere 2015
DlatIh = 1e-10*np.exp(-9e3*invR*invT) # from Ghesquiere 2015
Dbulk = 4.2e8*np.exp(-71.1e3*invR*invT) # from Livingston 1998
Dbulk_max = 1.49e-11*np.exp(-59.8e3*invR*(1/263-invT))  # cm^2/s, from Nasello et al 2007
Dvac = 1e-3*np.exp(-59.8e3*invR*invT) # from Hobbs, 1974
Dgb = 8.4*np.exp(-49e3*invR*invT) # from ....

# ----- Calculate the equilibrium vapor pressure of ice -----
Pev = 3.65e13 * np.exp(-6141.7 * invT) # [dyne/cm^2], (Fanale and Salvail, 1984)

# ----- Calculate the temp dependent density of ice -----
rho_H2O_I = (940.3 + 0.1143*(10*Pev*1E-7 -1) - 0.08585*T)*1E-3 # [g/cm^3], density (Ellsworth and Schubert, 1983)
rho_H2O_II = (1171.1 + 8.96E-2*Pev*1E-7 - 7.9E-2*T)*1E-3 # [g/cm^3], density (Ellsworth and Schubert, 1983)

# ----- Determine the thermal conductivity of ice -----
#kcem = [0.01, 0.1, 1.0, 2.0, 6.5] # [J/m/s/K], cementation TC
#kcem = map(lambda x: x*1e5, kcem)  # convert to [erg/cm/s/K]
kxice = (488.19/T+0.4685)*1e5 # [erg/cm/s/K], TC of Ih (hexagonal) crystalline water ice
kaice = 7.1e2*T # [erg/cm/s/K], TC of amorphous water ice
print 'The thermal conducivity [erg/cm/s/K] at {}K is {} for crystalline ice and {} for amorphous ice'.format(T[20], kxice[20], kaice[20])
 
# ----- Calculate Young's modulus, Poisson's ratio, Shear modulus, and Hertzian radius
rg_jkr = 25e-4
Fjkr = 3*pi*gamma_ss*rg_jkr
Ymod, Prat, Gmod, Rconjkr = fix_coeff(T, Fjkr, rg_jkr)
fignum +=1
#afig = plot_JKR(rg_jkr, T, Ymod, Prat, Rconjkr, fignum)

# ----- Define measured icy moon quantities -----

phi = 0.5 # Assume 50% porosity for now (all bodies)

# Define measured TI and determine TC  (slight underestimates?)
TIMout = 16*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Mimas outside the anomaly (< 16, J/m^2/s^0.5/K)
TIMin = 66*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Mimas inside the anomaly (66+/- 23, J/m^2/s^0.5/K)
TITout = 5*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Tethys outside the anomaly (5+/- 1, J/m^2/s^0.5/K)
TITin = 25*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Tethys inside the anomaly (25+/- 3, J/m^2/s^0.5/K)
TIDout = 8*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Dione on the trailing hemisphere (8 +/-1 J/m^2/s^0.5/K)
TIDin = 11*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Dione on the leading hemisphere (11+/- 4, J/m^2/s^0.5/K)

kMout = calc_TC(TIMout) # erg/cm/s/K, Effective TC at Mimas outside the anomaly
print 'For {}, TC={} [erg/cm/s/K]'.format('Mimas, in', kMout[20])
kMin = calc_TC(TIMin) # erg/cm/s/K, Effective TC at Mimas inside the anomaly
print 'For {}, TC={} [erg/cm/s/K]'.format('Mimas, out', kMin[20])
kTout = calc_TC(TITout) # erg/cm/s/K, Effective TC at Tethys outside the anomaly
print 'For {}, TC={} [erg/cm/s/K]'.format('Tethys, in', kTout[20])
kTin = calc_TC(TITin) # erg/cm/s/K, Effective TC at Tethys inside the anomaly
print 'For {}, TC={} [erg/cm/s/K]'.format('Tethys, out', kTin[20])
kDtrail = calc_TC(TIDout) # erg/cm/s/K, Effective TC at Dione on the trailing hemisphere
print 'For {}, TC={} [erg/cm/s/K]'.format('Dione, in', kDtrail[20])
kDlead = calc_TC(TIDin) # erg/cm/s/K, Effective TC at Dione on the leading hemisphere
print 'For {}, TC={} [erg/cm/s/K]'.format('Dione, out', kDlead[20])

r_Mimas = 25E-4 # cm, approximate grain radius (25 um); 
r_Tethys = 25E-4 # cm
r_Dione = 25E-4 # cm
# Mimas: 20-80um diameter, leading; 50-100um diameter, Hershel (Hendrix 2012)
# Mimas: D=58um Filacchione (2012) 
# Tethys: D=69um 
# Dione: D=59um

# Best Fitting Contour
Eflux_Mimas = 5.6E4 # MeV/cm^2/s, to Mimas thermal inertia anomaly boundary
Eflux_Tethys = 1.8E4 # MeV/cm^2/s, to Tethys thermal inertia anomaly boundary
Eflux_Dione = 5.6E4 # MeV/cm^2/s, to ..... (Mimas?) thermal inertia anomaly boundary

# Approximate Electron Fluxes
Phi_elec_Mimas = 8.2e3 # elec/cm^2/s, incident electron flux at Mimas
Phi_elec_Tethys = 2.1e2 # elec/cm^2/s, incident electron flux at Tethys
Phi_elec_Dione = 2.8e2 # elec/cm^2/s, incident electron flux at Dione

# Perform calculations for which body?
Bodies = ['Mimas']#, 'Tethys', 'Dione']
for body in Bodies:
    if body == 'Mimas':
        kin = kMin
        kout = kMout
        Phi_elec = Phi_elec_Mimas
        rg = r_Mimas
    elif body =='Tethys':
        kin = kTin
        kout = kTout
        Phi_elec = Phi_elec_Tethys
        rg = r_Tethys
    elif body =='Dione':
        kin = kDlead
        kout = kDtrail
        Phi_elec = Phi_elec_Dione
        rg = r_Dione
    else:
        print "Please enter a valid body"

    print 'For {}'.format(body)

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
    g = 1
    pc = 1/3
    p = 6*(1-phi)/pi
    SYpf = ((1-pc)/(p-pc))*g*rg*rg/pi
#    SYpf_lin = ((1-pc)/(p-pc))*g*rg/pi

    # ----- Calculate linear and quadratic minimum contact radius from the Wood model -----
    Rconminlin = (Apflin*kout*(1+0.5*phi)/(kxice*(1-phi)))
    Rconmaxlin = (Apflin*kin*(1+0.5*phi)/(kxice*(1-phi)))
    Rconminquad = np.sqrt(Apfquad*kout*(1+0.5*phi)/(kxice*(1-phi)))
    Rconmaxquad = np.sqrt(Apfquad*kin*(1+0.5*phi)/(kxice*(1-phi)))
    RconminSY = np.sqrt((kout/kxice)*SYpf)
    RconmaxSY = np.sqrt((kin/kxice)*SYpf)
#    RconminSY = (kout/kxice)*SYpf_lin # for linear
#    RconmaxSY = (kin/kxice)*SYpf_lin # for linear
    kg_out_lin = Apflin*kout*(1+0.5*phi)/(Rconjkr*2*(1-phi))
    kg_in_lin = Apflin*kin*(1+0.5*phi)/(Rconjkr*2*(1-phi))
    kg_out_quad = Apfquad*kout*(1+0.5*phi)/(Rconjkr*Rconjkr*(1-phi))
    kg_in_quad = Apfquad*kin*(1+0.5*phi)/(Rconjkr*Rconjkr*(1-phi))
    kg_out_SY = kout*SYpf/(Rconjkr*Rconjkr*4)
    kg_in_SY = kin*SYpf/(Rconjkr*Rconjkr*4)

    print "Rconminlin = %.2e" % (Rconminlin[20]*1E4), " um"
    print "Rconmaxlin = %.2e" % (Rconmaxlin[20]*1e4), " um"
    print "Rconminquad = %.2e" % (Rconminquad[20]*1E4), " um"
    print "Rconmaxquad = %.2e" % (Rconmaxquad[20]*1e4), " um"
    print "RconSYmax = %.2e" % (RconmaxSY[20]*1e4), " um"
    print "RconSYmin = %.2e" % (RconminSY[20]*1e4), " um"
    print "Rconjkr = %.2e" % (Rconjkr[20]*1e4), " um"
    print "kg_out_lin = {} [erg/cm^2/s^0.5/K]".format(kg_out_lin[20])
    print "kg_in_lin = {} [erg/cm^2/s^0.5/K]".format(kg_in_lin[20])
    print "kg_out_SY = {} [erg/cm^2/s^0.5/K]".format(kg_out_SY[20])
    print "kg_in_SY = {} [erg/cm^2/s^0.5/K]".format(kg_in_SY[20])

    # ----- Calculate radiation induced diffusion coefficients ----

    c_geo = 1 # geometric factor to determine sputtering/deposition relation
    beta = 6/(pi*pi) # factor for calculating radiation stuff
    # Calculate the radiation time constant
    invtau = dEdz*Phi_elec/(WH2O*nH2O)

    # Calculate the approximate length scale of an electron excitation event 
    l_c = 0.25*alat*(beta*DeltaE/UH2O) # cm, length-scale of excitation event
    EF_min = 1.0
    EF_max = 10.0
    # or
    Dz_min = 1*alat
    Dz_max = 10*alat
    alpha_v_min = 0.01
    alpha_v_max = 0.1 # GUESS, average number of point defects which escape cascade
    R_iv = alat/2 # GUESS, separation distance for spontaneous vacancy/interstitial recombination
    N_FP = 0.8*DeltaE/(2*UH2O)
    phi_sput = c_geo*l_c*invtau*nH2O # sputtered flux

    Jsurf = 4
    Jbulk = 6
    Ea_surf_min = (0.2)*1.602e-12 # (eV) erg, min surface diffusion activation energy
    Ea_surf_max = (0.5)*1.602e-12 # (eV) erg, max surface diffusion activation energy
    Ea_bulk_min = (0.6)*1.602e-12 # (eV) erg, min bulk activation energy
    Ea_bulk_max = (0.75)*1.602e-12 # (eV) erg, max bulk activation energy 
    Ea_gb_min = (0.6)*1.602e-12 # (eV) erg, min gb activation energy
    Ea_gb_max = (0.75)*1.602e-12 # (eV) erg, max gb activation energy 

    Dsurf_min = calc_Drad(Ea_surf_min, Dz_min, alpha_v_min, Jsurf)
    Dsurf_max = calc_Drad(Ea_surf_max, Dz_max, alpha_v_max, Jsurf)
    Dbulk_min = calc_Drad(Ea_bulk_min, Dz_min, alpha_v_min, Jbulk)
    Dbulk_max = calc_Drad(Ea_bulk_max, Dz_max, alpha_v_max, Jbulk)
    Dgb_min = calc_Drad(Ea_gb_min, Dz_min, alpha_v_min, Jbulk)
    Dgb_max = calc_Drad(Ea_gb_max, Dz_max, alpha_v_max, Jbulk)

    fignum +=1
#    a = plot_Diff(T, Dsurf, Dbulk, Dgb, Dvac, Dsurf_min, Dbulk_min, Dgb_min, fignum)
    fignum +=1
#    b = plot_Diff(T, Dsurf, Dbulk, Dgb, Dvac, Dsurf_max, Dbulk_max, Dgb_max, fignum)

    # ----- Determine the variation in the sintering rates based on the contact radius -----
    Rcon_var = np.logspace(-6, -3, num=100, base=10.0) # in cm
    invRcon_var = 1/Rcon_var
    NNc = 1 # Number of contacts
    Acon_var = pi*Rcon_var*Rcon_var/NNc

    Rnc = Rcon_var**2/(2*(rg-Rcon_var))
    invRnc = 1/Rnc
    theta = np.arctan(rg/(Rcon_var-Rnc))
    Km = invRcon_var - invRnc
    K3 = 2/rg

    # ----- Calculate radiation induced diffusion coefficients ----
    switches = ['80Kmin', '80Kmax']
    for dzswitch in switches:
        if dzswitch == '80Kmin':
            lineT = '--'
            Rmin1=Rconminlin[20]
            Rmin2=RconminSY[20]
            Rmax1=Rconmaxlin[20]
            Rmax2=RconmaxSY[20]
            EFc=EF_min
            invTc=invT[20]
            Tc=T[20]
            Pevc=Pev[20]
            Gmodc=Gmod[20]
            dt = 3.157e11

            Dsurfc=Dsurf[20] # Maeno and Ebinume
            Dgbc=Dgb[20]
            Dbulkc = Dbulk[20]
            Dvacc=Dvac[20]
            Deff_surfc = Dsurf_min
            Deff_bulkc = Dbulk_min[20]
            Deff_gbc = Dgb_min[20]
            Nvacc = N_vacH2O[20]

        if dzswitch == '80Kmax':
            lineT = '-'
            Rmin1=Rconminlin[20]
            Rmin2=RconminSY[20]
            Rmax1=Rconmaxlin[20]
            Rmax2=RconmaxSY[20]
            EFc=EF_max
            invTc=invT[20]
            Tc=T[20]
            Pevc=Pev[20]
            dt = 3.157e8

            Gmodc=Gmod[20]

            Dsurfc=Dsurf[20]
            Dgbc=Dgb[20]
            Dbulkc = Dbulk[20]
            Deff_surfc = Dsurf_max
            Deff_bulkc = Dbulk_max[20]
            Deff_gbc = Dgb_max[20]
            Nvacc = N_vacH2O[20]

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
        dVdtlat = 3*pi*Rcon_var*Dbulkc*gamma_sv*Omega*(K3-Km)*invkb*invTc

        vapDiff = Pevc*Rnc*theta*np.sqrt(Omega*invkb*invTc/(2*pi*rho_ice))
        dVdtvap = 2*pi*Rcon_var*vapDiff*gamma_sv*Omega*invkb*invTc*(K3-Km)

        dVdtgb = 16*pi*invRcon_var*Dgbc*gamma_sv*delta_gb*Omega*(1-0.5*K1*Rcon_var)*invkb*invTc
        dVdtlat_bound = 32*pi*invRcon_var*Rnc*theta*Dbulkc*gamma_sv*Omega*(1-0.5*Km*Rcon_var)*invkb*invTc
        dVdtlat_dloc = 8/9*pi*Rcon_var*Rcon_var*Rnc*theta*Nvacc*Dbulkc*gamma_sv*Omega*(-Km-3*Gmodc*Rcon_var/(2*rg*gamma_sv))*invkb*invTc

        dVdttherm_tot = dVdtsurf + dVdtlat + dVdtvap + dVdtgb + dVdtlat_bound + dVdtlat_dloc
        dRdttherm_tot = dVdttherm_tot/Vol_to_rad

    # ----- Calculate Radiation Enhanced volumetric sintering rates -----
        SurfEF = (delta_surf/d2)*gamma_sv*Omega*(K3-K2)*invkb*invTc*EFc
        dVdtradsurf = 3*pi*Rcon_var*Deff_surfc*SurfEF

        BulkEF = Omega*(K3-Km)*invkb*invTc*EFc
        dVdtradlat = 3*pi*Rcon_var*Deff_bulkc*gamma_sv*BulkEF

        SputEF = gamma_sv*Omega*invkb*invTc*(K3-Km)*EFc
        dVdtradsput = 2*pi*Rcon_var*theta*Rnc*(phi_sput/nH2O)*SputEF

        GBEF = delta_gb*gamma_sv*Omega*(1-0.5*K1*Rcon_var)*invkb*invTc*EFc
        dVdtradgb = 16*pi*invRcon_var*Dgbc*GBEF

        lat_boundEF = gamma_sv*Omega*(1-0.5*Km*Rcon_var)*invkb*invTc*EFc
        dVdtradlat_bound = 32*pi*invRcon_var*Rnc*theta*Dbulkc*lat_boundEF

        lat_dlocEF = gamma_sv*Omega*(-Km-3*Gmodc*Rcon_var/(2*rg*gamma_sv))*invkb*invTc*EFc
        dVdtradlat_dloc = 8/9*pi*Rcon_var*Rcon_var*Rnc*theta*Nvacc*Dbulkc*lat_dlocEF

        dVdtrad_tot = dVdtradsurf + dVdtradlat + dVdtradsput + dVdtradgb + dVdtradlat_bound + dVdtradlat_dloc
        dRdtrad_tot = dVdtrad_tot/Vol_to_rad

        fignum+=1
        a = plot_VSR(rg,Rcon_var,dVdtrad_tot,dVdtradsurf,dVdtradlat,dVdtradsput,dVdtradgb,dVdtradlat_bound,dVdtradlat_dloc,fignum,lineT)

        Amax = pi*Rmax2*Rmax2/NNc
        Rmax = np.sqrt(Amax/pi)
        Amin = pi*Rmin2*Rmin2/NNc
        Rmin = np.sqrt(Amin/pi)
        Rcur = Rmin
        Acur = pi*Rcur*Rcur/NNc
     #   print 'The min rad is %.2e' % Rmin2, ' and the max is %.2e' % Rmax2
        ttot=0
        R_t = []
        timesteps = []

        print 'Rcon, dRdtrad_tot, dRdttherm_tot'
    #    for i in range(10):
        while Rcur<Rmax:
            index = find_nearest(Rcon_var,Rcur)
            DR = dRdtrad_tot[index]+dRdttherm_tot[index]
            Rcur+=DR*dt
            Acur = pi*Rcur*Rcur/NNc
            R_t.append(Rcur)
            ttot+=dt
            timesteps.append(ttot)
    #        print SurfEF[index], BulkEF[index], SputEF[index]
    #        print '%.3e' % Rcur, '%.3e' % dRdtrad_tot[index], '%.3e' % dRdttherm_tot[index]

        R_tarr=np.array(R_t,dtype='float')
        ttot_years = ttot*stoyear

        print 'the total time for %d K' % Tc, ' is %.2e' % ttot_years, ' years' 

        fig = plt.figure(4)
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
