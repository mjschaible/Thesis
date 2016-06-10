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

    num_runs, run_thermo = lammps_read.log_read(filename)
    if run_thermo[0].descrip[2] is not None:
        print 'The excited molecule is {}, the position is {}'.format(run_thermo[0].descrip[3], run_thermo[0].descrip[4])
        print 'The molecule velocity is {:.2f}'.format(run_thermo[0].descrip[5])
    else:
        print 'No excited molecule specified', run_thermo[0].descrip[2]

    steplen = run_thermo[0].descrip[0]
    run_info= run_thermo[0].descrip

    filedens = filepath+'.dens'
    run_dens = lammps_read.data_read(filedens, steplen)

    filemsd = filepath+'.msd'
    run_msd = lammps_read.data_read(filemsd, steplen)

#    filemsd_gRad = filepath+'_gRad.msd'
#    num_msd_gRad, run_msd_gRad = lammps_read.data_read(filemsd_gRad, run_thermo[0].descrip)

#    filemsd_gPKA = filepath+'_gPKA.msd'
#    num_msd_gPKA, run_msd_gPKA = lammps_read.data_read(filemsd_gPKA, run_thermo[0].descrip)

    filerdf = filepath+'.rdf'
    run_rdf=lammps_read.rdf_read(filerdf, steplen)

    filecom = filepath+'.com'
    run_com = lammps_read.com_read(filecom, steplen)
    msd_com = lammps_read.find_commsd(run_com, run_info)

    filedump = filepath+'.dump'
    run_dump = lammps_read.dump_read(filedump, steplen)
    #dump_eng = lammps_read.find_dumpeng(run_dump, run_info)
    
    #----- -----
    lammps_plot.log_plots(run_thermo)

    nrows=2
    fig, axs = plt.subplots(nrows,1,sharey='row',squeeze=False)
    #lammps_plot.msd_plots(axs[0,:], run_dens, 'Density (g/cm^3)')
    lammps_plot.msd_plots(axs[0,:], run_msd, 'All msd (A^2)')
    lammps_plot.msd_plots(axs[1,:], msd_com, 'Shell msd (A^2)')

#    lammps_plot.msd_plots(axs[1,:], run_msd_gRad, run_param, 'gRad msd (A^2)')
#    lammps_plot.msd_plots(axs[2,:], run_msd_gPKA, run_param, 'PKA msd (A^2)')
#    lammps_plot.rdf_plots(run_rdf, run_param)
    plt.show()

    cont = input('Enter 0 to end: ')
    
raise SystemExit()
