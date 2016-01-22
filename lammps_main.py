import numpy as np
import matplotlib.pyplot as plt
plt.interactive(True)
import glob

from Tkinter import Tk
from tkFileDialog import askopenfilename

import lammps_read
import lammps_plot

# ------------ Begin main program ----------
cont = 1

num_runs=0
run_param=[]
run_thermo=[]

while cont != 0:
# show an "Open" dialog box and return the path to the selected file
    Tk().withdraw()
    filename = askopenfilename()
    filepath = filename.split(".")[0]
    print filepath
    
    num_runs, run_param, run_thermo = lammps_read.log_read(filename)

    filemsd = filepath+'.msd'
    num_msd, run_msd = lammps_read.data_read(filemsd)

    filemsd_gRad = filepath+'_gRad.msd'
#    num_msd_gRad, run_msd_gRad = lammps_read.data_read(filemsd_gRad)

    filemsd_gPKA = filepath+'_gPKA.msd'
#    num_msd_gPKA, run_msd_gPKA = lammps_read.data_read(filemsd_gPKA)

    filedens = filepath+'.dens'
    num_dens, run_dens = lammps_read.data_read(filedens)

    lammps_plot.log_plots(run_thermo, run_param)
    nrows=2
    fig, axs = plt.subplots(nrows,num_msd,sharey='row',sharex='col',squeeze=False)
    lammps_plot.msd_plots(axs[0,:], run_msd, run_param, 'All msd (A^2)')
    lammps_plot.msd_plots(axs[1,:], run_dens, run_param, 'Density (g/cm^3)')
#    lammps_plot.msd_plots(axs[1,:], run_msd_gRad, run_param, 'gRad msd (A^2)')
#    lammps_plot.msd_plots(axs[2,:], run_msd_gPKA, run_param, 'PKA msd (A^2)')
    
    plt.show()

    cont = input('Enter 0 to end: ')
    
raise SystemExit()
