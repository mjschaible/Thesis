import numpy as np
import matplotlib.pyplot as plt

def fix_coeff(T, Fjkr, rg):
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

def plot_Diff(T, Dsurf, Dbulk, Dgb, Dvac, Deffsurf, Deffbulk, Deffgb, lineT, fignum):
    fig = plt.figure(fignum)
    ax2= fig.add_subplot(1,1,1)
    lineT = '-'
    if lineT == '-':
        lbl1='Rad. surface'
        lbl2='Rad. bulk'
        lbl7='Rad. GB'
        lbl3='Therm. surface'
        lbl4='Therm, bulk'
        lbl5='Therm, GB'
#        lbl6='Therm, vac'
    else:
        lbl1=None
        lbl2=None
        lbl3=None
        lbl4=None
        lbl7=None
        lbl5=None
        lbl6=None
    lw=3
    ax2.plot([1000/60,1000/110],[Deffsurf,Deffsurf],label = lbl1,lw=lw, linestyle = lineT, color='r')
    ax2.plot(1000/T,Deffbulk, linestyle = lineT, label = lbl2,lw=lw, color='g')
    ax2.semilogy([12.5,12.5],[1e-30,1e-22], color='k', linewidth=1, linestyle = '--')
    ax2.semilogy(1000/T,Deffgb, linestyle=lineT, label =lbl7,lw=lw, color='c')
#    ax2.semilogy(T,DsurfIc, label = 'DsurfIc')
#    

    ax2.set_xlabel('1000/T,  [$1/K$]', fontsize = 16)
    ax2.set_ylabel('Diffusion Coeff. [$cm^2/s$]', fontsize = 16)
    ax2.set_xlim([1000/70, 1000/110])
#    ax2.text(14.7, 2.2e-25, r'(a)', fontsize=20,color='k')
#    if lineT == '-':
    ax2.semilogy(1000/T,Dsurf,lw=lw, linestyle=':', marker='v', color='r',label = lbl3)
    ax2.semilogy(1000/T,Dbulk,lw=lw, linestyle=':', marker='o', color='g',label = lbl4)
    ax2.semilogy(1000/T,Dgb,lw=lw, linestyle=':', marker='>', color='c',label = lbl5)
#        ax2.semilogy(1000/T,Dvac, linewidth=2, linestyle='-.', marker='x',color='m',label = lbl6)
    ax2.set_ylim([Deffsurf/50,Deffsurf*2])
    ax2.legend(loc = 2, prop={'size':14})
        
    ax3= ax2.twiny()
    ax2Ticks = ax2.get_xticks()
 
    ax3.set_xticks(ax2Ticks)
    ax3.set_xbound(ax2.get_xbound())
    upperx = 1000/ax2Ticks
    upx = np.around(upperx)
    ax3.set_xticklabels(upx)
    ax3.set_xlim(ax3.get_xlim()[::-1])
    ax3.set_xlabel('Temperature, [$K$]', fontsize = 16)

    plt.savefig('./DiffusionRates.png',dpi=600)#.format(limit)

    return

def plot_VSR(Rcon,dVdtrad_tot,dVdtradsurf,dVdtradlat,dVdtradsput,dVdtradgb,dVdtradbound,fignum,lineT,colnum, Rmin,Rmax):
    import numpy as np
    import matplotlib.pyplot as plt
    stoyear = 3.16888e-8*1e-6 # convert second to Myr
    fig = plt.figure(fignum)
    ax2 = fig.add_subplot(1,1,colnum)
    if colnum==1 and lineT == '-':
        lbl1='Rad. surface'
        lbl2='Rad. bulk'
        lbl3='Rad. sputter'
        lbl4='Rad. Tot'
        lbl5='Rad. GB'
    else:
        lbl1=None
        lbl2=None
        lbl3=None
        lbl4=None
        lbl5=None
    lw=2

    ax2.semilogy(Rcon,dVdtrad_tot*stoyear,lw=4,ls=lineT, c='k', label = lbl4)
    
    ax2.semilogy(Rcon,dVdtradsurf*stoyear,lw=lw,ls=lineT, c='r',label = lbl1)#, marker='v')
    ax2.semilogy(Rcon,dVdtradlat*stoyear,lw=lw,ls=lineT, c='g',label = lbl2)#, marker='o')
    ax2.semilogy(Rcon,dVdtradsput*stoyear,lw=lw,ls=lineT, c='b',label = lbl3)#, marker='x')
    ax2.semilogy(Rcon,dVdtradgb*stoyear,lw=lw,ls=lineT, c='c',label = lbl5)
#    ax2.semilogy(Rcon,dVdtradbound*stoyear,linestyle = lineT,label = 'Rad. lat/bound')
#    ax2.semilogy(Rcon,dVdtraddloc*stoyear,linestyle = lineT,label = 'Rad. dislocation')
    
    ax2.set_xlim([0,0.04])
#    ax2.xlabel('Contact/Grain radius (Rcon/rg)', fontsize = 16)
    ax2.semilogy([Rmin,Rmin],[1e-3,1e6], color='k', linewidth=1, linestyle = '--')
    ax2.text(0.4*Rmin, 2e4, r'$R_{con,min}$', rotation='vertical', fontsize=20,color='k')
    ax2.semilogy([Rmax,Rmax],[1e-3,1e6], color='k', linewidth=1, linestyle = '--')
    plt.xticks(np.arange(0,0.041,0.01))

    if colnum==1 and lineT == '-':
        ax2.set_title('$r_g = 5 \mu m$', fontsize = 16, y=1.03)
        ax2.set_ylabel('Sintering Timescale $(\\tau_{sint})$ [Myr]', fontsize = 16)
        plt.suptitle('Contact/Grain radius ($R_{con}/r_g$)', fontsize = 16, y=.05)
        box = ax2.get_position()
        ax2.legend(loc = 4, prop={'size':12})
        ax2.set_ylim([1e-2,1e5])
        ax2.text(0.9*Rmax, 2e4, r'$R_{con,max}$', rotation='vertical', fontsize=20,color='k')
#        plt.text(-0.015, 2e5, r'(b)', fontsize=20,color='k')
    if colnum==2 and lineT == '-':
        ax2.set_title('$r_g = 25 \mu m$', fontsize = 16, y=1.03)
        plt.setp(ax2.get_yticklabels(), visible=False)
        ax2.set_ylim([1e-2,1e5])
        ax2.text(1.02*Rmax, 0.2, r'$R_{con,max}$', rotation='vertical', fontsize=20,color='k')

    plt.savefig('./Rad_Sint_Rates.png',dpi=600) #.format(limit)

    return
