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

def calc_TC(TI, phi):
    TC = (TI*TI)/(c_H2O*rho_H2O_I*(1-phi))
    return TC

def calc_Drad(Ea, UH2O, J):
    bleta = 0.8
    Nmob = bleta * DeltaE / (Ea) # number of surface molecules mobilized per ionization
    Dz = alat*Nmob*0.25
    d = 10*Dz # approximate average sink diffusion length scale    
    N_FP = bleta*DeltaE/(UH2O)

    Deff_mc = Nmob*Dz*Dz*invtau/J
    Deff_esd_mid = 2*alpha_v*d*d*N_FP*invtau
    Deff_esd_low = np.sqrt(alpha_v*N_FP*Dvac*invtau/(pi*nH2O*R_iv))

    if J==4:
        Deff = Deff_mc
    elif J==6:
        Deff = Deff_mc + Deff_esd_low# + Deff_esd_mid
    else:
        'Please select either Bulk or Surface diffusion'
    return Deff

def calc_timescale(Rmin, Rmax):
    Amax = pi*Rmax*Rmax/NNc
    Rmax = np.sqrt(Amax/pi)
    Amin = pi*Rmin*Rmin/NNc
    Rmin = np.sqrt(Amin/pi)
    Rcur = Rmin
    Acur = pi*Rcur*Rcur/NNc
#   print 'The min rad is %.2e' % Rmin2, ' and the max is %.2e' % Rmax2
    ttot=0
    R_t = []
    timesteps = []
    index = find_nearest(Rcon_var,Rmax)
    Vcem_max = V_cem[index]
    index = find_nearest(Rcon_var,Rmin)
    Vcem_min = V_cem[index]
    Vcem=Vcem_max-Vcem_min
    DV = dVdtrad_tot[index]+dVdttherm_tot[index]
    dt = Vcem/DV/1E5
#    for i in range(10):
    while Acur<Amax:
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

    return ttot_years

# ----------- Begin main program -------------------------------------------
T = np.linspace(60,110,num=51) # Define temperature range for calculations
invT = 1/T # Calculate the inverse of temperature
Tnorm = T/273 # Calculate the normalized temprature (normalize to melting temp)

# --------------- Define necessary model parameters -------------------------
pi = np.pi
kb = 1.38e-16 # erg/K, Boltzman constant
stoyear = 3.16888e-8 # year/sec
ev2erg = 1.602e-12 # erg/eV
invkb = 1/kb

# For now assume an average energy loss for all depths
# Should eventually use Penelope results to give vs. depth
dEdz = (1e6)*ev2erg # (eV)erg/cm, energy deposition rate at Mimas

# Define water ice parameters
alat = 3e-8 # cm, approximate lattice spacing for H2O
gamma_ss = 65 # erg/cm^2, ice/ice (grain boundary) surface energy (Ketcham and Hobbs, 1969)
gamma_sv = 300 # erg/cm^2, ice/vapor (free surface) surface energy (Henry, 2003)
# = 190 (Gundlach et al., 2011)
# = 109 (Ketcham and Hobbs, 1969)
gamma_sv_max = 370 # erg/cm^2 (Aumatell and Wurm, 2014; Henry, 2003)
delta_gb = 9.04e-8 # cm, grain boundary thickness
#delta_surf = 3.19e-8 # cm, effective surface (atomic layer) thickness
delta_surf=2*alat 
R_iv = alat/2 # GUESS, separation distance for spontaneous vacancy/interstitial recombination
rho_ice = 0.934 # g/cm^3, ~density of water ice @ 80K
nH2O = rho_ice*6.022e23/18 # cm^3, molar volume
Omega = 1/nH2O # Approximate molecular volume

Jsurf = 4 # Degrees of freedom for surface diffusion mechanisms
Jbulk = 6 # Degrees of freedom for bulk diffusion processes

WH2O = (27)*ev2erg # (eV)erg, average energy deposited in an excitation event
DeltaE = 0.2*WH2O
c_geo = 1 # geometric factor to determine sputtering/deposition relation
beta = 6/(pi*pi) # factor for calculating radiation stuff
alpha_v = 1 # GUESS, average number of point defects which escape cascade

# cohesive energy of water [Watkins etal, 2011]
UH2O_surf_min = 0.1*ev2erg
UH2O_surf_max = 0.5*ev2erg
UH2O_bulk_min = 0.6*ev2erg
UH2O_bulk_max = 0.75*ev2erg

Ea_surf_min = (0.1)*ev2erg # (eV) erg, min surface diffusion activation energy [Nasello et al., 2007]
Ea_surf_max = (0.4)*ev2erg # (eV) erg, max surface diffusion activation energy [Nie et al., 2009]
Ea_bulk_min = (0.6)*ev2erg # (eV) erg, min bulk activation energy [Nasello et al., 2007]
Ea_bulk_max = (0.75)*ev2erg # (eV) erg, max bulk activation energy [Livingston et al., 1998]
Ea_gb_min = (0.5)*ev2erg # (eV) erg, min gb activation energy [Goldsby and Kohlstedt, 2001]
Ea_gb_max = (0.5)*ev2erg # (eV) erg, max gb activation energy 

N_vacH2O = nH2O*np.exp(-UH2O_bulk_min*invkb*invT) # Equilibrium vacancy concentration

# ----- Calculate the temp dependent specific heat capacity of ice -----
c_H2O = (0.09+0.00749*T)*1E7  # [erg/g/K], (Gutierez et al, 2001)

# ----- Calculate thermal diffusion coefficients for ice -----
invR = 1/8.314 # [R] = [J/mol K], gas constant
Dsurf = 6.1e-3* np.exp(-44.3e3*invR*invT) # cm^2/s, from Kouchi 1994
#Dsurf = 3e2 * np.exp(-44.13e3*invR*invT) # cm^2/s, from Maeno and Ebinume, 1983
#Dsurf = 1.4e-4*np.exp(-48.2e3*invR*invT) # cm^2/s, max estimate from Nie et al 2009
Dsurf_HT = 1.4e-4*np.exp(-23.2e3*invR*invT)  # cm^2/s, high temp [Nasello et al 2007]
DlatLDA = 7e-6*np.exp(-15e3*invR*invT) # from Ghesquiere 2015
DlatIh = 1e-10*np.exp(-9e3*invR*invT) # from Ghesquiere 2015
Dbulk = 4.2e8*np.exp(-71.1e3*invR*invT) # from Livingston 1998
Dbulk_HT = 1.49e-11*np.exp(-59.8e3*invR*(1/263-invT))  # cm^2/s, high temp  [Nasello et al 2007]
Dvac = 1e-3*np.exp(-59.8e3*invR*invT) # from Hobbs, 1974
Dgb = 8.4*np.exp(-49e3*invR*invT) # from Goldsby and Kohlstedt 2001

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
#print 'The thermal conducivity [erg/cm/s/K] at {}K is {} for crystalline ice and {} for amorphous ice'.format(T[20], kxice[20], kaice[20])

# ----- Define measured icy moon quantities -----
# Define measured TI and determine TC  (slight underestimates?)
# Howett et al. 2011, 2012, 2014
TIMout = 16*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Mimas outside the anomaly (< 16, J/m^2/s^0.5/K)
Mouterr = 4*1E3 # [erg/cm^2/s^0.5/K], Approx error at Mimas outside the anomaly
TIMin = 66*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Mimas inside the anomaly (66+/-23, J/m^2/s^0.5/K)
Minerr = 23*1E3 # [erg/cm^2/s^0.5/K], Approx error at Mimas inside the anomaly

TITout = 5*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Tethys outside the anomaly (5+/-1, J/m^2/s^0.5/K)
Touterr = 1*1E3 # [erg/cm^2/s^0.5/K], Approx error at Tethys outside the anomaly
TITin = 25*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Tethys inside the anomaly (25+/-3, J/m^2/s^0.5/K)
Tinerr = 3*1E3 # [erg/cm^2/s^0.5/K], Approx error at Tethys inside the anomaly

TIDout = 8*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Dione on the trailing hemisphere (8+/-1 J/m^2/s^0.5/K)
Douterr = 1*1E3 # [erg/cm^2/s^0.5/K], Approx error at Tethys outside the anomaly
TIDin = 11*1E3 # [erg/cm^2/s^0.5/K], Measured TI at Dione on the leading hemisphere (11+/-4 J/m^2/s^0.5/K)
Dinerr = 4*1E3 # [erg/cm^2/s^0.5/K], Approx error at Dione inside the anomaly

# Best Fitting Contour [Howett et al., 2011, 2012, 2014]
Eflux_Mimas = 5.6E4 # MeV/cm^2/s, to Mimas thermal inertia anomaly boundary
Eflux_Tethys = 1.8E4 # MeV/cm^2/s, to Tethys thermal inertia anomaly boundary
Eflux_Dione = 1.8E4 # MeV/cm^2/s, to ..... (Mimas?) thermal inertia anomaly boundary

# Approximate Volumetric energy deposition rate (taken from Penelope results)
# ED = (dE/dz) * Electron Fluxes = (dE/dz)*Eflux/<Eelec>
# <Eelec> = average electron energy [Paranicas et al., 2014]
ED_Mimas = 5e10*ev2erg # eV/cm^3/s, incident electron flux at Mimas
ED_Tethys = 4e10*ev2erg # eV/cm^3/s, incident electron flux at Tethys
ED_Dione = 4e10*ev2erg # eV/cm^3/s, incident electron flux at Dione


# Mimas: 20-80um diameter, leading; 50-100um diameter, Hershel (Hendrix 2012)
# Mimas: D=58um Filacchione (2012) 
# Tethys: D=69um
# Dione: D=59um

# Select the grain sizes to be considered
r_Mimas = [5E-4]#, 25e-4]

# Select the regolith porosity
phi_Mimas = [0.5]#, 0.65, 0.85]

# Select the activation energies to be used
switches = ['80KminEa']#,'80KmaxEa','80KBoost']

# Decide bulk ice thermal conductivity
kbulk = [kxice]#,kaice]

# Perform calculations for which body?
Bodies = ['Mimas']#,'Tethys','Dione'] #

for body in Bodies:
    if body == 'Mimas':
        TIin = TIMin
        TIout = TIMout
        TIinerr = Minerr
        TIouterr = Mouterr
        ED = ED_Mimas
    elif body =='Tethys':
        TIin = TITin
        TIout = TITout    
        TIinerr = Tinerr
        TIouterr = Touterr
        ED = ED_Tethys
    elif body =='Dione':
        TIin = TIDin
        TIout = TIDout
        TIinerr = Dinerr
        TIouterr = Douterr
        ED = ED_Dione
    else:
        print "Please enter a valid body"

    for rg in r_Mimas:
        for phi in phi_Mimas:
            for kice in kbulk:
                print
                print '--------------------------------------------------------------------------'
                print 'For {}, rg = {}um, phi = {}, and kice={}[erg/cm/s/K] @ 80K'.format(body, rg*1e4, phi, kice[20])

                # Calculate effective TC using measured TI
                kout = calc_TC(TIout, phi) # erg/cm/s/K, outside the anomaly
                kerr_min = calc_TC(TIout-TIouterr, phi)
                kerr_max = calc_TC(TIout+TIouterr, phi)
                kouterr_m = kout-kerr_min
                kouterr_p = kerr_max-kout
                print 'Outside anomaly, keff={:.1f}+{:.1f}\-{:.1f} [erg/cm/s/K]'.format(kout[20], kouterr_p[20], kouterr_m[20])
                kin = calc_TC(TIin, phi) # erg/cm/s/K, inside the anomaly
                kerr_min = calc_TC(TIin-TIinerr, phi)
                kerr_max = calc_TC(TIin+TIinerr, phi)
                kinerr_m = kin-kerr_min
                kinerr_p = kerr_max-kin
                print 'Inside anomaly, keff={:.1f}+{:.1f}\-{:.1f} [erg/cm/s/K]'.format(kin[20], kinerr_p[20], kinerr_m[20])

                # ----- Calculate the approximate "Hertzian" contact radius for the given body grain size -----
                # Young's modulus, Poisson's ratio, Shear modulus, and Hertzian radius
                Fjkr = 3*pi*gamma_ss*rg
                Ymod, Prat, Gmod, Rconjkr = fix_coeff(T, Fjkr, rg)
                hfig=1
                #afig = plot_JKR(rg_jkr, T, Ymod, Prat, Rconjkr, hfig)
                print "Rconjkr = {:.3f} um".format(Rconjkr[20]*1e4)

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

                # ----- Calculate linear and quadratic minimum contact radius from the Wood model -----
                Rconminlin = (Apflin*kout*(1+0.5*phi)/(kice*(1-phi)))
                Rconmaxlin = (Apflin*kin*(1+0.5*phi)/(kice*(1-phi)))
                Rconminquad = np.sqrt(Apfquad*kout*(1+0.5*phi)/(kice*(1-phi)))
                Rconmaxquad = np.sqrt(Apfquad*kin*(1+0.5*phi)/(kice*(1-phi)))

                # ------ Define Sirono and Yamamoto model parameters ------
                g = 1
                pc = 1/3
                p = 6*(1-phi)/pi
                SYpf = ((1-pc)/(p-pc))*g*rg*rg/pi
            #    SYpf_lin = ((1-pc)/(p-pc))*g*rg/pi
                # ----- Calculate linear and quadratic minimum contact radius from the Sirono and Yamamoto model -----
                RconminSY = np.sqrt((kout/kice)*SYpf)
                RminSYerr_m = RconminSY-np.sqrt(((kout-kouterr_m)/kice)*SYpf)
                RminSYerr_p = np.sqrt(((kout+kouterr_p)/kice)*SYpf)-RconminSY
                RconmaxSY = np.sqrt((kin/kice)*SYpf)
                RmaxSYerr_m = RconmaxSY-np.sqrt(((kin-kinerr_m)/kice)*SYpf)
                RmaxSYerr_p = np.sqrt(((kin+kinerr_p)/kice)*SYpf)-RconmaxSY

                chi_min = 3*Nc*(RconminSY/rg)**4 /16
                chi_min_m = chi_min-3*Nc*((RconminSY-RminSYerr_m)/rg)**4 /16
                chi_min_p = 3*Nc*((RconminSY+RminSYerr_p)/rg)**4 /16-chi_min 
                chi_max = 3*Nc*(RconmaxSY/rg)**4 /16
                chi_max_m = chi_max-3*Nc*((RconmaxSY-RmaxSYerr_m)/rg)**4 /16
                chi_max_p = 3*Nc*((RconmaxSY+RmaxSYerr_p)/rg)**4 /16-chi_max

            #    print "Rconminlin = %.2e" % (Rconminlin[20]*1E4), " um"
            #    print "Rconmaxlin = %.2e" % (Rconmaxlin[20]*1e4), " um"
            #    print "Rconminquad = %.2e" % (Rconminquad[20]*1E4), " um"
            #    print "Rconmaxquad = %.2e" % (Rconmaxquad[20]*1e4), " um"
                print "RconSY, out = {:.0f}+/-{:.0f} nm (chi={:.0e} +/-{:.0e})".format(RconminSY[20]*1e7, RminSYerr_p[20]*1E7, chi_min[20], chi_min_m[20])
                print "RconSY, in = {:.0f}+/-{:.0f} nm (chi={:.0e} +/-{:.0e})".format(RconmaxSY[20]*1e7, RmaxSYerr_p[20]*1E7, chi_max[20], chi_max_m[20])

                # ----- Determine the variation in the sintering rates based on the contact radius -----
                Rcon_ = np.logspace(-5, -1, num=1000, base=10.0) # in cm
                Rcon_var=Rcon_*rg
                invRcon_var = 1/Rcon_var
                NNc = 1 # Number of contact
                Acon_var = pi*Rcon_var*Rcon_var/NNc

                V_g = 4*pi*rg*rg*rg/3
                Rcon_rg = Rcon_var/rg
                V_cem = 3*Nc*Rcon_rg*Rcon_rg*Rcon_rg*Rcon_rg*V_g/16

                Rnc = Rcon_var**2/(2*(rg-Rcon_var))
                invRnc = 1/Rnc
                theta = np.arctan(rg/(Rcon_var-Rnc))
                Km = invRcon_var - invRnc
                Kg = 2/rg

                # Calculate the radiation time constant
                invtau = ED/(WH2O*nH2O)

                for dzswitch in switches:
                    if dzswitch == '80KmaxEa':
                        lineT = '-.'
                        Ea_surf = Ea_surf_max
                        Ea_bulk = Ea_bulk_max
                        Ea_gb = Ea_gb_max
                        UH2O_surf = UH2O_surf_max
                        UH2O_bulk = UH2O_bulk_max
                        boostie = 1
                    elif dzswitch == '80KminEa':
                        lineT = '-'
                        Ea_surf = Ea_surf_min
                        Ea_bulk = Ea_bulk_min
                        Ea_gb = Ea_gb_min
                        UH2O_surf = UH2O_surf_min
                        UH2O_bulk = UH2O_bulk_min
                        boostie = 1
                    elif dzswitch == '80KBoost':
                        lineT = '-.'
                        Ea_surf = Ea_surf_min
                        Ea_bulk = Ea_bulk_min
                        Ea_gb = Ea_gb_min
                        UH2O_surf = UH2O_surf_min
                        UH2O_bulk = UH2O_bulk_min
                        boostie = 2
                    else:
                        print 'There is nothing to compute here.'

                    # ----- Calculate radiation induced diffusion coefficients ----
                    Dsurf_rad = calc_Drad(Ea_surf, UH2O_surf, Jsurf)
                    Dbulk_rad = calc_Drad(Ea_bulk, UH2O_bulk, Jbulk)
                    Dgb_rad = calc_Drad(Ea_gb, UH2O_bulk, Jbulk)

                    Deff_surfc = Dsurf_rad
                    Deff_bulkc = Dbulk_rad[20]
                    Deff_gbc = Dgb_rad[20]

                    dfig=hfig+1
                    if  body == 'Mimas' and dzswitch=='80KminEa' and rg==5e-4 and phi==0.5 and kice[20]==kxice[20]:
                        a = plot_Diff(T, Dsurf, Dbulk, Dgb, Dvac, Dsurf_rad, Dbulk_rad, Dgb_rad, lineT, dfig)

                    invTc=invT[20]
                    Tc=T[20]
                    Pevc=Pev[20]
                    Gmodc=Gmod[20]

                    Dsurfc=Dsurf[20]
                    Dgbc=Dgb[20]
                    Dbulkc = Dbulk[20]
                    Dvacc=Dvac[20]
                    Nvacc = N_vacH2O[20]

                    A = Dgbc*delta_gb/(Dsurfc*delta_surf)
                    arr = A*Rnc*invRcon_var
                    DKK1n = 3*arr*arr + 4*arr+1
                    DKK1 = 3*arr+1.5-1.5*np.sqrt(DKK1n)
                    K1 = Km/np.sqrt(1+DKK1)
                    K2 = K1*(1+DKK1)
                    d2 = Rnc/(1 + np.sqrt((4/3)*((K2-K1)/K2)))
                    d1 = Rnc-d2

                    Vol_to_rad = 2*pi*theta*Rcon_var*Rnc

                # ----- Calculate the thermal volumetric sintering rates -----
                    dVdtsurf = 3*pi*Rcon_var*Dsurfc*delta_surf*gamma_sv*Omega*(Kg-K2)*invkb*invTc/d2
                    dVdtlat = 3*pi*Rcon_var*Dbulkc*gamma_sv*Omega*(Kg-Km)*invkb*invTc

                    phiVap = Pevc*np.sqrt(Omega*invkb*invTc/(2*pi*rho_ice))*gamma_sv*Omega*invkb*invTc*(Kg-Km)
                    dVdtvap = 2*pi*Rcon_var*Rnc*theta*phiVap

                    dVdtgb = 16*pi*invRcon_var*Dgbc*gamma_sv*delta_gb*Omega*(1-0.5*K1*Rcon_var)*invkb*invTc
                    dVdtlat_bound = 32*pi*invRcon_var*Rnc*theta*Dbulkc*gamma_sv*Omega*(1-0.5*Km*Rcon_var)*invkb*invTc
                    dVdtlat_dloc = 8/9*pi*Rcon_var*Rcon_var*Rnc*theta*Nvacc*Dbulkc*gamma_sv*Omega*(-Km-3*Gmodc*Rcon_var/(2*rg*gamma_sv))*invkb*invTc

                    dVdttherm_tot = dVdtsurf + dVdtlat + dVdtvap + dVdtgb + dVdtlat_bound + dVdtlat_dloc
                    dRdttherm_tot = dVdttherm_tot/Vol_to_rad

                # ----- Calculate Radiation Enhanced volumetric sintering rates -----
                    SurfEF = (delta_surf/d2)*gamma_sv*Omega*invkb*invTc*(Kg-K2)
                    dVdtradsurf = 3*pi*Rcon_var*Deff_surfc*SurfEF

                    BulkEF = gamma_sv*Omega*invkb*invTc*(Kg-Km)
                    dVdtradlat = 3*pi*Rcon_var*Deff_bulkc*BulkEF

                    l_c = 0.25*alat*(beta*DeltaE/UH2O_surf) # cm, length-scale of excitation event
                    phi_sput = c_geo*l_c*invtau*nH2O # sputtered flux
                    SputEF = (phi_sput/nH2O)*l_c*(Kg-Km)
                    dVdtradsput = 2*pi*Rcon_var*Rnc*theta*SputEF

                    GBEF = delta_gb*gamma_sv*Omega*(1-0.5*K1*Rcon_var)*invkb*invTc
                    dVdtradgb = 16*pi*invRcon_var*Deff_gbc*GBEF

                    lat_boundEF = gamma_sv*Omega*(1-0.5*Km*Rcon_var)*invkb*invTc
                    dVdtradlat_bound = 32*pi*invRcon_var*Rnc*theta*Dbulkc*lat_boundEF

                    lat_dlocEF = gamma_sv*Omega*(-Km-3*Gmodc*Rcon_var/(2*rg*gamma_sv))*invkb*invTc
                    dVdtradlat_dloc = 8/9*pi*Rcon_var*Rcon_var*Rnc*theta*Nvacc*Dbulkc*lat_dlocEF

                    dVdtrad_tot = dVdtradsurf + dVdtradlat + dVdtradsput + dVdtradgb + dVdtradlat_bound + dVdtradlat_dloc
                    dRdtrad_tot = dVdtrad_tot/Vol_to_rad

                    rrat = Rcon_var/rg
                    Rmin=RconminSY[20]
                    Rminerr = RminSYerr_m[20]
                    Rmin_p = Rmin + Rminerr
                    Rmin_m = Rmin - Rminerr
                    Rmax=RconmaxSY[20]
                    Rmaxerr = RmaxSYerr_m[20]
                    Rmax_p = Rmax+Rmaxerr
                    Rmax_m = Rmax-Rminerr
                    print "The contact to grain ratio increases from {:.5f} to {:.5f}".format(Rmin/rg, Rmax/rg)

                    vfig=dfig+1
                    if body == 'Mimas' and rg==25e-4 and phi==0.5 and  kice[20]==kxice[20]:
                        colnum=2
                        a = plot_VSR(rrat,V_cem/dVdtrad_tot,V_cem/dVdtradsurf,V_cem/dVdtradlat,V_cem/dVdtradsput,V_cem/dVdtradgb,V_cem/dVdtradlat_bound,vfig,lineT, colnum, Rmin/rg, Rmax/rg)
                    elif body == 'Mimas' and rg==5e-4 and phi==0.5 and  kice[20]==kxice[20]:
                        colnum=1
                        a = plot_VSR(rrat,V_cem/dVdtrad_tot,V_cem/dVdtradsurf,V_cem/dVdtradlat,V_cem/dVdtradsput,V_cem/dVdtradgb,V_cem/dVdtradlat_bound,vfig,lineT, colnum, Rmin/rg, Rmax/rg)

                    ttot = calc_timescale(Rmin, Rmax)
                    ttot_p = calc_timescale(Rmin_p, Rmax_p)
                    ttot_m = calc_timescale(Rmin_m, Rmax_m)

                    print 'The total time for {}, {} is {:.0e}--{:.0e} Myr'.format(body, dzswitch, (ttot_m)/1E6, (ttot_p)/1E6) 

                    fig = plt.figure(vfig+1)
                    ax = fig.add_subplot(1,1,1)
            #        line, = ax.plot(timesteps,R_tarr*1e4, linewidth=2, color='r', label = 'Rad. sint @ %d' %Tc, linestyle = lineT)
            #        ax.annotate('Rmax @ %d K' % Tc,xy=(ttot,Rcur*1e4),xytext=(ttot,Rcur+2.4),arrowprops=dict(facecolor='black',shrink=0.05), horizontalalignment='center')
                    if dzswitch == 1:
                        ax.set_ylabel('Contact Radius [um]')
                        ax.set_xlabel('Timestep [s]')
                    if dzswitch == 3:
                        ax.legend(loc = 3, prop={'size':18})
                #        plt.savefig('./sint_timescales.png',dpi=600)

                #plt.show()
