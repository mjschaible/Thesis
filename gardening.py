'''Program to calculate the gardening rates on the Saturnian satellites''' 

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
matplotlib.rcParams.update({'font.size': 16})

def find_nearest(A,B):
    index = (np.abs(A-B)).argmin()
    return index

# Uniform parameters
vinf = 9.5E5 # cm/s, interplanetary dust velocity
phi = 0.6 # approximate regolith porosity
rho_reg = 0.934*(1-phi) # g/cm^3
stoyear = 3.16888e-8 # year/s, multiply seconds by this to get years
t0 = 1.7E5 # years, mixing timescale

# Saturn and orbital radii
rs = 6.03E9 # cm
rM = 3.2*rs
rT = 4.4*rs
rD = 6.6*rs
#rEu
#rGa

# Escape velocities
vescM = 14.03E5 # cm/s 
vescT = 11.96E5 # cm/s
vescD = 9.77E5 # cm/s
#vescEu = 
#vescGa =

# approximate ejecta yields
YidpM = 1.8E4
YidpT = 1.2E4
YidpD = 9.8E3
YidpJ = 3E4 

YeringM = 33
YeringT = 19
YeringD = 14

# Iterplanetary Dust fluxes
sigma_inf_Sat = 8E-19 #g/cm^2/s
sigma_inf_Jup = 5E-17 #g/cm^2/s
#Circumplanetary Dust fluxes
sigma_Sat_M = 62E-15 #g/cm^2/s 
sigma_Sat_T = 230E-15 #g/cm^2/s 
sigma_Sat_D = 27E-15 #g/cm^2/s

# Define the time scale
treg = np.logspace(3, 8, base=10, num=100) # years, regolith age

bodies = ['Mimas', 'Tethys', 'Dione']
for body in bodies:
    if body == 'Mimas':
        rb = rM
        vesc = vescM
        sigma_Sat = sigma_Sat_M
        Yidp = YidpM
        Yering = YeringM
    if body == 'Tethys':
        rb = rT
        vesc = vescT
        sigma_Sat = sigma_Sat_M
        Yidp = YidpT
        Yering = YeringM
    if body == 'Dione':
        rb = rD
        vesc = vescD
        sigma_Sat = sigma_Sat_M
        Yidp = YidpD
        Yering = YeringM
    sigma_inf = sigma_inf_Sat

    fg = 0.5*np.sqrt(1 + (vesc*vesc)/(vinf*vinf))*(np.sqrt(1 + (vesc*vesc)/(vinf*vinf))+np.sqrt(1 + (vesc*vesc)/(vinf*vinf) - (rs*rs)/(rb*rb)*(1 + (vesc*vesc)/(vinf*vinf))))  # Gravitational focussing factor
#    sigma_ej = fg*sigma_inf*Y0*(rb/(1.8*rs))**(-0.8) # g/cm^2/s, Cuzzi and Estrada ejecta flux 
    sigma_ej_idp = fg*sigma_inf*Yidp
    sigma_ej_ering = sigma_Sat*Yering

    sigma_ej = sigma_ej_idp# + sigma_ej_ering

    g0 = sigma_ej/rho_reg # cm/s, regolith growth rate

    h_t = (g0/stoyear)*treg*(1+treg/t0)**(-0.55) # cm

    index = find_nearest(h_t, 1)
    tovr = treg[index]

    print 'The idp flux at {} is {:.2e} g/cm^2/s'.format(body, sigma_ej_idp/Yidp)
    print 'The regolith growth rate is {:.2e} um/year'.format(g0*1E4/stoyear)
    print 'The overturn timescale at 1cm is {:.2f} Myr'.format(tovr/1E6)
    fig = plt.figure(1)
    ax2 = fig.add_subplot(1,1,1)
    ax2.set_title('Gardening depths for the saturnian icy moons')
    ax2.loglog(treg, h_t, linewidth=2, label = body)
    box = ax2.get_position()
    ax2.set_xlabel('Time [years]', fontsize = 16)
    ax2.set_ylabel('Depth [cm]', fontsize = 16)
    ax2.legend(loc = 1, prop={'size':12})

plt.savefig('./Gardening_depths.png',dpi=600)
plt.show()




