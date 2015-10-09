def plot_sputFile(descrip, species, fsteps, specFlu1, specFlu2, specFlu3, av0, av1, av2):
    import numpy as np
    import matplotlib.pyplot as plt
    params = descrip.split()
    simEng = params[0] + 'eV'
#    print simEng
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
        
def plot_sputYld(elemName, energies, Flux1, Flux2, Flux3, simCount):
    import numpy as np
    import matplotlib.pyplot as plt
 
    fig = plt.figure(2)
    ax1 = fig.add_subplot(111)
    #ax1.set_title("Sputtered Fluence")    
    ax1.set_xlabel('Energy')
    ax1.set_ylabel('Sputtered Fluence')
    
    totFlu = [x+y for x,y in zip(Flux2,Flux3)]
    print Flux2, Flux3, totFlu
    ax1.loglog(energies,totFlu, label="{0}".format(elemName[1]))
#    ax1.plot(energies,Flux3, c='b', label="{0}".format(elemName[2]))    
    if simCount == 0:
        leg = ax1.legend()

def plot_sputExpt(elemName, energies, sputyld, simCount):
    import numpy as np
    import matplotlib.pyplot as plt
 
    fig = plt.figure(2)
    ax1 = fig.add_subplot(111)
    #ax1.set_title("Sputtered Fluence")    
    ax1.set_xlabel('Energy')
    ax1.set_ylabel('Sputtered Fluence')
    
    ax1.scatter(energies,sputyld, c='r', label="{0}".format(elemName))
#    ax1.plot(energies,Flux3, c='b', label="{0}".format(elemName[2]))    
    if simCount == 0:
        leg = ax1.legend()
        
