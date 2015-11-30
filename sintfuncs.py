def fix_coeff(T, Fjkr, rg):
    import numpy as np
    pi = np.pi

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

    # ----- Determine the compliance coefficients -----
    ltSij = np.zeros((6,6,len(T)))

    Ymod = np.zeros_like(T)
    Prat = np.zeros_like(T)
    Gmod = np.zeros_like(T)
    Rconjkr = np.zeros_like(T)
    # Calculate Young's modulus in erg/cm^3 from compliance coefficients
    Escale = 1e6 # Convert bar to erg/cm^3 (for mechanical constants)
    jkrpow = float(1)/3
    tempSij = np.zeros_like(Cij)

    for i in range(len(T)):
        ltSij[:,:,i] = np.linalg.inv(Cij[:,:,i])

    for i in range(len(T)):
        Ymod[i] = 1/ltSij[0,0,i]*Escale  # Youngs Modulus
        Prat[i] = -ltSij[0,1,i]/ltSij[0,0,i] # Poisson's ratio
        Gmod[i] = 0.5*Ymod[i]/(1+Prat[i]) # Shear Modulus
        Rconjkr[i] = ((0.75*(1-Prat[i]*Prat[i])*rg*Fjkr)/Ymod[i])**(jkrpow)

    return Ymod, Prat, Gmod, Rconjkr

def plot_JKR(rg, T, Ymod, Prat, Rcon, fignum):
    import matplotlib.pyplot as plt
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.set_title('For {}um water ice grains'.format(rg*1e4))

    ax1.plot(T,Ymod, label = 'Ih water ice (Hobbs, 1974)')
    ax1.set_ylabel('Youngs modulus [erg/cm^3]')
 #   ax3.axis([60,110,5e-3,1])
    ax1.legend(loc = 3, prop={'size':14})

    ax2.plot(T,Prat, label = 'Ih water ice (Hobbs, 1974)')
    ax2.set_ylabel('Poissons ratio')
#    ax2.axis([60,110,5e-3,1])
    ax2.legend(loc = 3, prop={'size':14})

    ax3.plot(T,Rcon*1E4, label = 'Rcon,jkr')
    ax3.set_ylabel('Contact Radius [um]')
    ax3.set_xlabel('Temperature [K]')
#    ax3.axis([60,110,5e-3,1])
    ax3.legend(loc = 3, prop={'size':14})

    plt.savefig('./Hertzian.png',dpi=600)

    return

def plot_Diff(T, Dsurf, Dbulk, Dgb, Dvac, Deffsurf, Deffbulk, Deffgb, fignum):
    import numpy as np
    import matplotlib.pyplot as plt
    fig = plt.figure(fignum)

    ax2= fig.add_subplot(1,1,1)
    #ax2.semilogy(T,DsurfIc, label = 'DsurfIc')
    ax2.semilogy(1000/T,Dsurf, linewidth=2, color='c', linestyle='-.', label = 'Dtherm, surf')
    ax2.semilogy(1000/T,Dbulk, linewidth=2, color='m', linestyle='-.', label = 'Dtherm, bulk')
    ax2.semilogy(1000/T,Dgb, linewidth=2, linestyle='-.', label = 'Dtherm, gb')
    ax2.semilogy(1000/T,Dvac, linewidth=2, linestyle='-.', label = 'Dtherm, vac')

    ax2.semilogy([1000/60,1000/110],[Deffsurf,Deffsurf],label = 'Drad, surf', linestyle='-', color='r')
    ax2.semilogy(1000/T,Deffbulk, linestyle='-', label = 'Drad, bulk', color='k')
#    ax2.semilogy(1000/T,Deffgb, linestyle='-', label = 'Drad, gb', color='g')

    ax2.set_xlabel('1000/T,  [1/K]', fontsize = 16)
    ax2.set_ylabel('Diffusion Coeff. [cm^2/s]', fontsize = 16)
    ax2.set_xlim([1000/70, 1000/105])
    ax2.legend(loc = 2, prop={'size':14})

    ax3= ax2.twiny()
    ax2Ticks = ax2.get_xticks()
    ax3Ticks = ax2Ticks
 
    ax3.set_xticks(ax3Ticks)
    ax3.set_xbound(ax2.get_xbound())
    upperx = 1000/ax3Ticks
    upx = np.around(upperx)
    ax3.set_xticklabels(upx)
    ax3.set_xlim(ax3.get_xlim()[::-1])
    ax3.set_xlabel('Temperature, [K]', fontsize = 16)

    if fignum == 2:
        limit = 'min'
    else:
        limit = 'max'
    plt.savefig('./DiffusionRates{}.png'.format(limit),dpi=600)
    
    return

def plot_VSR(rg,Rcon_var,dVdtrad_tot,dVdtradsurf,dVdtradlat,dVdtradsput,dVdtradgb,dVdtradbound,dVdtraddloc, fignum, lineT):
    import numpy as np
    import matplotlib.pyplot as plt
    fig = plt.figure(fignum)
    ax2 = fig.add_subplot(1,1,1)

#    ax2.semilogy([Rmin2/rg,Rmin2/rg],[1e-36,1e-26], color='k', linewidth=2, linestyle = '-.')
#    ax2.semilogy([Rmax2/rg,Rmax2/rg],[1e-36,1e-26], color='k', linewidth=2, linestyle = '-.')

    ax2.semilogy(Rcon_var/rg,dVdtrad_tot, linewidth=2,linestyle=lineT, color='b', label = 'Rad. Tot')
    
    ax2.semilogy(Rcon_var/rg,dVdtradsurf, linestyle = lineT, label = 'Rad. surface')
    ax2.semilogy(Rcon_var/rg,dVdtradlat, linestyle = lineT, label = 'Rad. lattice')
    ax2.semilogy(Rcon_var/rg,dVdtradsput, linestyle = lineT, label = 'Rad. sputter')
    ax2.semilogy(Rcon_var/rg,dVdtradgb, linestyle = lineT, label = 'Rad. boundary ')
#    ax2.semilogy(Rcon_var/rg,dVdtradbound, linestyle = lineT, label = 'Rad. lat/bound')
    ax2.semilogy(Rcon_var/rg,dVdtraddloc, linestyle = lineT, label = 'Rad. dislocation')
    
    box = ax2.get_position()
    ax2.set_xlabel('Contact/Grain radius (Rcon/rg)', fontsize = 16)
    ax2.set_ylabel('Volumetric Sintering Rate [cm^3/s]', fontsize = 16)
    ax2.set_xlim([0,0.055])
    ax2.set_ylim([1e-33,1e-26])
    ax2.legend(loc = 1, prop={'size':12})
    
    if fignum == 4:
        limit = 'min'
    else:
        limit = 'max'

    plt.savefig('./Rad_Sint_Rates_{}.png'.format(limit),dpi=600)

    return
