import numpy as np
import matplotlib.pyplot as plt
from Tkinter import Tk
from tkFileDialog import askopenfilename
from operator import attrgetter

kcalToeV = 0.0433641 # convert kcal/mole to eV/atom

class LogData (object):

    def __init__(self,potential,step,temp,tempave,pe,peave,ke,keave,etot):
        self.potential=potential
        self.step=step
        self.temp=temp
        self.tempave=tempave
        self.pe=pe
        self.peave=peave
        self.ke=ke
        self.keave=keave
        self.etot=etot

def generate_plots(n):
    a = n
    b = 3
#    print "a\t=\t%d\nb\t=\t%d\na*b\t=\t%d\nn\t=\t%d" % (a,b,a*b,n)
    fig, axs = plt.subplots(b,a,figsize=(a*b,b*b), sharex='col', sharey='row')
    axs.shape = (b,a)

    for i in range(n):
        axs[0,i].plot(runs[i].step,runs[i].pe*kcalToeV,label = HeadCol[3])
        axs[0,i].plot(runs[i].step,runs[i].peave*kcalToeV*num_molec)
        if i == 2:
            pemin=min(runs[i].peave)
            pemax=max(runs[i].peave)
            axs[0,i].set_xscale([pemin-20,pemax+100])
        if i > 2:
            axs[0,i].set_xscale([pemin-20,pemax+100])
        plt.setp(axs[0,i].get_xticklabels(), visible=False)
        axs[0,i].locator_params(axis='x', tight=True, nbins=2)

        axs[1,i].plot(runs[i].step,runs[i].ke*kcalToeV,label = HeadCol[5])
        axs[1,i].plot(runs[i].step,runs[i].keave*kcalToeV*num_molec)
        plt.setp(axs[1,i].get_xticklabels(), visible=False)
        axs[1,i].locator_params(axis='x', tight=True, nbins=2)
 
        axs[2,i].plot(runs[i].step,runs[i].etot*kcalToeV,label = HeadCol[7])
        axs[2,i].set_xlabel('Step #  (%s fs)' % timesteps[i])
        axs[2,i].ticklabel_format(axis='x', style = 'sci')
        axs[2,i].locator_params(axis='x', tight=True, nbins=2)

        if i == 0:
            axs[0,i].set_ylabel('%s (eV)' % HeadCol[3])
            axs[1,i].set_ylabel('%s (eV)' % HeadCol[5])
            axs[2,i].set_ylabel('%s (eV)' % HeadCol[7])
            plt.suptitle("The pair potential used is %s" % potential)
        if i > 0:
            plt.setp(axs[0,i].get_yticklabels(), visible=False)
            plt.setp(axs[1,i].get_yticklabels(), visible=False)
            plt.setp(axs[2,i].get_yticklabels(), visible=False)
        plt.subplots_adjust(wspace=0.001)
        axs[2,i].set_ylim(-1700,-1600)

    plt.show()
    return fig


# ------------ Begin main program ----------
cont = 1
atomline = 0
num_runs=0
runs=[]

data_start=[]
data_end=[]
timesteps=[]
markers=[]

while cont != 0:
#    name = raw_input('Enter filename: ')
#    filename.append(name)

# show an "Open" dialog box and return the path to the selected file
    Tk().withdraw()
    filename = askopenfilename() 
#    filename = 'log.lammps'
    print(filename)
    counter = 1

#filename='log.lammps'

    with open(filename,'r') as logfile:
        contents=logfile.read()

    lines=contents.split('\n')

    for line in lines:
        if '# create groups ###' in line:
            atomline = counter
        if atomline > 0 and counter == atomline+2:
            column = line.split()
            num_molec=float(column[0])
            num_atoms=num_molec*3
            print num_atoms
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

    num_cur=len(data_start)-num_runs

    data_vals=['']*num_cur

    for n in range(num_cur):
        data_vals[n]=lines[data_start[n+num_runs]:data_end[n+num_runs]-1]

        step=[]
        temp=[]
        tempave=[]
        pe=[]
        peave=[]
        ke=[]
        keave=[]
        etot=[]

        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()

            step.append(vals[0])
            temp.append(vals[1])
            tempave.append(vals[2])
            pe.append(vals[3])
            peave.append(vals[4])
            ke.append(vals[5])
            keave.append(vals[6])
            etot.append(vals[7])
        steparr=np.array(step,dtype='float')
        temparr=np.array(temp,dtype='float')
        tempavearr=np.array(tempave,dtype='float')
        pearr=np.array(pe,dtype='float')
        peavearr=np.array(peave,dtype='float')
        kearr=np.array(ke,dtype='float')
        keavearr=np.array(keave,dtype='float')
        etotarr=np.array(etot,dtype='float')
    #    runs.append(LogData(potential,step,pe,ke,etot))
        runs.append(LogData(potential,steparr,temparr,tempavearr,pearr,peavearr,kearr,keavearr,etotarr))

    num_runs+=num_cur
    cont = input('Enter 0 to end: ')

#print runs[3].step
a = generate_plots(num_runs)
