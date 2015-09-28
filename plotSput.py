def plot_sputFile(descrip, species, fsteps, specFlu1, specFlu2, specFlu3, av0, av1, av2):
    import numpy as np
    import matplotlib.pyplot as plt
    params = descrip.split()
    simEng = params[0] + 'eV'
    print simEng
    fig = plt.figure(1)
#    fig.suptitle("The simulation run is %s" %descrip)
    ax1 = fig.add_subplot(111)
    #ax1.plot(fsteps,specFlu1, label = species[0])
    #ax1.plot(fsteps,specFlu2, label = species[1])
    ax1.plot(fsteps,specFlu3, label = simEng)
    #ax1.plot(fsteps,av0)
    #ax1.plot(fsteps,av1)
    ax1.plot(fsteps,av2)
    #ax1.set_xlim(.001,20)
    #ax1.set_ylim(.001,0.2)                    
    ax1.set_ylabel('Yield (atom/ion)')
    ax1.set_xlabel('Fluence (x10^16)')
    leg = ax1.legend()
    
'''   fig2 = plt.figure(3)
    fig2.suptitle("Averaged simulation results for %s" %descrip)
    ax2 = fig2.add_subplot(111)
    ax2.set_ylim(0.001,0.2)
    #ax2.set_xlim(0.001,20)
    ax2.set_ylabel('Yield (atom/ion')
    ax2.set_xlabel('Fluence (x10^16)')
    
    stepave=[]
    
    #for i in range(1,len(av0)):
    #    step = i*numave
    #    stepave.append(fsteps[step-1])
    #print len(fsteps), len(av0)
    
    ax2.plot(fsteps,av0, label = species[0])
    ax2.plot(fsteps,av1, label = species[1])
    ax2.plot(fsteps,av2, label = species[2])
    
    leg2 = ax2.legend() '''
    
    #plt.show()

