'''Program to calculate the gardening rates on the Saturnian satellites''' 

import numpy as np
import matplotlib.pyplot as plt
import itertools
import matplotlib as mpl

#mpl.use("TkAgg")
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.titlesize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['axes.labelpad'] = 2.5
mpl.rcParams['xtick.labelsize']='large'
mpl.rcParams['ytick.labelsize']='large'
mpl.rcParams['legend.handlelength']=2.5

def find_nearest(A,B):
    index = (np.abs(A-B)).argmin()
    return index

# Uniform parameters
vinf = 9.5E5 # cm/s, interplanetary dust velocity
phi = 0.5 # approximate regolith porosity
rho_reg = 0.934*(1-phi) # g/cm^3
stoyear = 3.16888e-8 # year/s, multiply seconds by this to get years
t0 = 1.7E5 # years, mixing timescale

# Saturn and orbital radii
rs = 6.03E9 # cm
rM = 3.08*rs
rT = 4.89*rs
rD = 6.26*rs
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
treg = np.logspace(4, 7, base=10, num=100) # years, regolith age

bodies = [ 'Mimas', 'Tethys', 'Dione']
color=iter(plt.cm.brg(np.linspace(0,1,len(bodies))))

Body=0
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
        sigma_Sat = sigma_Sat_T
        Yidp = YidpT
        Yering = YeringT
    if body == 'Dione':
        rb = rD
        vesc = vescD
        sigma_Sat = sigma_Sat_D
        Yidp = YidpD
        Yering = YeringD
    sigma_inf = sigma_inf_Sat

    fg = 0.5*np.sqrt(1 + (vesc*vesc)/(vinf*vinf))*(np.sqrt(1 + (vesc*vesc)/(vinf*vinf))+np.sqrt(1 + (vesc*vesc)/(vinf*vinf) - (rs*rs)/(rb*rb)*(1 + (vesc*vesc)/(vinf*vinf))))  # Gravitational focussing factor

#    sigma_ej = fg*sigma_inf*Y0*(rb/(1.8*rs))**(-0.8) # g/cm^2/s, Cuzzi and Estrada ejecta flux 
    sigma_ej_idp = fg*sigma_inf*Yidp
    sigma_ej_ering = sigma_Sat*Yering
#    sigma_ej = sigma_Sat + sigma_ej_idp #

#    g_Er = sigma_ej_ering/rho_reg
    g_Er = sigma_Sat/rho_reg
    g_IDP = sigma_ej_idp/rho_reg
#    g0 = sigma_ej/rho_reg # cm/s, regolith growth rate
    
    h_Er=(g_Er/stoyear)*treg#*(1+treg/t0)**(-0.55)
    h_IDP=(g_IDP/stoyear)*treg*(1+treg/t0)**(-0.55)
    h_tot=h_Er+h_IDP
#    h_t = (g0/stoyear)*treg*(1+treg/t0)**(-0.55) # cm

    reg_dep = 1 # cm
    i_Er = find_nearest(h_Er, reg_dep)
    i_IDP = find_nearest(h_IDP, reg_dep)
#    index = find_nearest(h_t, reg_dep)
    t_Er = treg[i_Er]
    t_IDP = treg[i_IDP]
#    tovr = treg[index]

    print '-------------------------------------------'
    print 'The IDP flux at {} is {:.2e} g/cm^2/s'.format(body, sigma_ej_idp/Yidp)
    print 'The IDP overturn timescale at 1cm is {:.4f} Myr'.format(t_IDP/1E6)
#    print 'The IDP overturn rate is {:.2e} um/year'.format(g_IDP*1E4/stoyear)
    print 'The E-ring growth rate is {:.2e} um/year'.format(g_Er*1E4/stoyear)
    print 'The E-ring growth timescale at 1cm is {:.4f} Myr'.format(t_Er/1E6)
#    print 'The renewal timescale at 1cm is {:.4f} Myr'.format(tovr/1E6)
    fig = plt.figure(1)
    ax2 = fig.add_subplot(1,1,1)
    ax2.set_title('Gardening and growth rates for the icy Saturnian moons')
#    ax2.loglog(treg, h_t, linewidth=2, label = body)
    c=next(color)
    ax2.loglog( treg*1e-6,h_IDP, linewidth=2, ls=':', c=c)
    ax2.loglog( treg*1e-6,h_Er, linewidth=2, ls='--', c=c)
    ax2.loglog( treg*1e-6,h_tot, linewidth=2, ls='-',label = body, c=c)
    box = ax2.get_position()
    ax2.set_xlabel('Time [Myr]', fontsize = 16)
    ax2.set_ylabel('Depth [cm]', fontsize = 16)
    ax2.set_ylim([.01,100])
    ax2.legend(loc = 2, prop={'size':18})
    Body+=1
plt.savefig('./Gardening_depths.png',dpi=600)
plt.show()




