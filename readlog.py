import numpy as np
import matplotlib.pyplot as plt
kcalToeV = 0.0433641 # convert kcal/mole to eV/atom

class LogData (object):

    def __init__(self,potential,step,pe,ke,etot):
        self.potential=potential
        self.step=step
        self.pe=pe
        self.ke=ke
        self.etot=etot

def generate_plots(n):
    a = n
    b = 3
#    print "a\t=\t%d\nb\t=\t%d\na*b\t=\t%d\nn\t=\t%d" % (a,b,a*b,n)
    fig, axs = plt.subplots(b,a,figsize=(a*b,b*b), sharex='col', sharey='row')
    axs.shape = (b,a)

    for i in range(n):
        axs[0,i].plot(runs[i].step,runs[i].pe*kcalToeV,label = HeadCol[2])
        plt.setp(axs[0,i].get_xticklabels(), visible=False)
        axs[0,i].locator_params(axis='x', tight=True, nbins=2)

        axs[1,i].plot(runs[i].step,runs[i].ke*kcalToeV,label = HeadCol[3])
        plt.setp(axs[1,i].get_xticklabels(), visible=False)
        axs[1,i].locator_params(axis='x', tight=True, nbins=2)
 
        axs[2,i].plot(runs[i].step,runs[i].etot*kcalToeV,label = HeadCol[4])
        axs[2,i].set_xlabel('Step #  (%s fs)' % timesteps[i])
        axs[2,i].ticklabel_format(axis='x', style = 'sci')
        axs[2,i].locator_params(axis='x', tight=True, nbins=2)

        if i == 0:
            axs[0,i].set_ylabel('%s (eV/atom)' % HeadCol[2])
            axs[1,i].set_ylabel('%s (eV/atom)' % HeadCol[3])
            axs[2,i].set_ylabel('%s (eV/atom)' % HeadCol[4])
            plt.suptitle("The pair potential used is %s" % potential)
        if i > 0:
            plt.setp(axs[0,i].get_yticklabels(), visible=False)
            plt.setp(axs[1,i].get_yticklabels(), visible=False)
            plt.setp(axs[2,i].get_yticklabels(), visible=False)
        plt.subplots_adjust(wspace=0.001)
#        axs[2,i].set_ylim(-1700,-1600)

    plt.show()
    return fig


# ------------ Begin main program ----------
filename='log.lammps'

with open(filename,'r') as logfile:
    contents=logfile.read()

lines=contents.split('\n')

data_start=[]
data_end=[]
counter=1
timesteps=[]
markers=[]
for line in lines:
    if 'pair_style' in line:
        pairs=line.split()
        potential=pairs[1]
    if 'timestep' in line:
        column = line.split()
        timesteps.append(float(column[1]))
    if 'PPPM' in line:
        markers.append(counter)
    if 'Step Temp' in line:
        data_start.append(counter)
        if len(data_start)-1==0:
            HeadCol = line.split()
    if 'Loop time' in line:
        data_end.append(counter)
    if len(data_start) > len(timesteps):
        if len(timesteps)==0:
            timesteps.append(0)
        else:
            timesteps.append(timesteps[len(timesteps)-1])
    counter+=1
    if counter == len(lines) and len(data_end) < len(data_start):
        data_end.append(counter)
        
num_runs=len(data_start)

runs=[]
data_vals=['']*num_runs

for n in range(num_runs):
    data_vals[n]=lines[data_start[n]:data_end[n]-1]

    step=[]
    pe=[]
    ke=[]
    etot=[]
    for m in range(len(data_vals[n])):
        vals=data_vals[n][m].split()
#        print vals
        step.append(vals[0])
        pe.append(vals[2])
        ke.append(vals[3])
        etot.append(vals[4])
    steparr=np.array(step,dtype='float')
    pearr=np.array(pe,dtype='float')
    kearr=np.array(ke,dtype='float')
    etotarr=np.array(etot,dtype='float')
#    runs.append(LogData(potential,step,pe,ke,etot))
    runs.append(LogData(potential,steparr,pearr,kearr,etotarr))

#print runs[3].step
a = generate_plots(num_runs)
