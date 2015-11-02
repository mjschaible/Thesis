def plot_sputFile(energy, species, fsteps, specFlu1, specFlu2, specFlu3, av0, av1, av2, nF):
    import numpy as np
    import matplotlib.pyplot as plt

    if species[1]=='Si' and species[2]=='O':
        target = 'SiO2'
    descrip = species[1] + '->' + target
    simEng = energy + 'eV'
#    print simEng
    fig = plt.figure(nF)
    fig.suptitle("The simulation run is {0}".format(descrip))
    ax1 = fig.add_subplot(111)
    
    #ax1.plot(fsteps,specFlu1, label = species[0])
    #ax1.plot(fsteps,specFlu2, label = species[1])
    ax1.plot(fsteps,specFlu3, label = simEng + ', ' + species[2])
    #ax1.plot(fsteps,av0)
    #ax1.plot(fsteps,av1)
    ax1.plot(fsteps,av2)
    #ax1.set_xlim(.001,20)
    #ax1.set_ylim(.001,0.2)                    
    ax1.set_ylabel('Yield (atom/ion)')
    ax1.set_xlabel('Fluence (x10^16)')
    leg = ax1.legend()
        
def plot_sputYld(elemName, energies, Flux1, Flux2, Flux3, sC):
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    #ax1.set_title("Sputtered Fluence")    
    ax1.set_xlabel('Energy')
    ax1.set_ylabel('Sputtered Fluence')
    
    totFlu = [x+y for x,y in zip(Flux2,Flux3)]
    if sC % 2 ==0:
        LT = '-'
        ax1.semilogx(energies,totFlu, label="{0}".format(elemName[1]), linestyle = LT)
    else:
        LT = '-.'
        ax1.semilogx(energies,totFlu, linestyle = LT)
    #ax1.plot(energies,Flux3, c='b', label="{0}".format(elemName[2]))    
    leg = ax1.legend()

def plot_sputExpt(elemName, energies, sputyld, simCount):
    import numpy as np
    import matplotlib.pyplot as plt
 
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    #ax1.set_title("Sputtered Fluence")    
    ax1.set_xlabel('Energy')
    ax1.set_ylabel('Sputtered Fluence')
    
    ax1.scatter(energies,sputyld, c='r', label="{0}".format(elemName))
#    ax1.plot(energies,Flux3, c='b', label="{0}".format(elemName[2]))    
    if simCount == 0:
        leg = ax1.legend()
        
