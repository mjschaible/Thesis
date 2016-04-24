import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
plt.interactive(True)
import glob
import itertools

from Tkinter import Tk
from tkFileDialog import askdirectory
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
    if not run_thermo[0].descrip[2]:
        print 'The excited molecule is {}, the position is {}, and the velocity is {:.2f}'\
            .format(run_thermo[0].descrip[2], run_thermo[0].descrip[3], run_thermo[0].descrip[4])
    else:
        print 'No excited molecule specified'

    filedens = filepath+'.dens'
    num_dens, run_dens = lammps_read.data_read(filedens, run_thermo[0].descrip)

    filemsd = filepath+'.msd'
    num_msd, run_msd = lammps_read.data_read(filemsd, run_thermo[0].descrip)

    filemsd_gRad = filepath+'_gRad.msd'
    num_msd_gRad, run_msd_gRad = lammps_read.data_read(filemsd_gRad, run_thermo[0].descrip)

    filemsd_gPKA = filepath+'_gPKA.msd'
    num_msd_gPKA, run_msd_gPKA = lammps_read.data_read(filemsd_gPKA, run_thermo[0].descrip)

    filerdf = filepath+'.rdf'
    num_rdf, run_rdf=lammps_read.rdf_read(filerdf, run_thermo[0].descrip)

    filecom = filepath+'.com'
    num_com, run_com=lammps_read.com_read(filecom, run_thermo[0].descrip)
    msd_com = lammps_read.find_commsd(run_param, run_com, num_com)

    #----- -----
    lammps_plot.log_plots(run_thermo, run_param)

    nrows=3
    fig, axs = plt.subplots(nrows,num_msd,sharey='row',sharex='col',squeeze=False)
#    lammps_plot.msd_plots(axs[0,:], run_dens, run_param, 'Density (g/cm^3)')
    lammps_plot.msd_plots(axs[0,:], run_msd, run_param, 'All msd (A^2)')
    lammps_plot.msd_plots(axs[1,:], run_msd_gRad, run_param, 'gRad msd (A^2)')
    lammps_plot.msd_plots(axs[2,:], run_msd_gPKA, run_param, 'PKA msd (A^2)')
    
    lammps_plot.rdf_plots(run_rdf, run_param)
    lammps_plot.msd2_plots(run_com, msd_com, run_param)
    plt.show()

    cont = input('Enter 0 to end: ')
    
raise SystemExit()

'''    


'''
