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
    if run_param[0].mol_pka != 0:
        print 'The excited molecule is {}, the position is {}, and the velocity is {:.2f}'\
            .format(run_param[0].mol_pka, run_param[0].pos_pka, run_param[0].vel_pka)
    else:
        print 'No excited molecule specified'


    filedens = filepath+'.dens'
    num_dens, run_dens = lammps_read.data_read(filedens)

    filemsd = filepath+'.msd'
    num_msd, run_msd = lammps_read.data_read(filemsd)

#    filemsd_gRad = filepath+'_gRad.msd'
#    num_msd_gRad, run_msd_gRad = lammps_read.data_read(filemsd_gRad)

#    filemsd_gPKA = filepath+'_gPKA.msd'
#    num_msd_gPKA, run_msd_gPKA = lammps_read.data_read(filemsd_gPKA)

#    filerdf = filepath+'.rdf'
#    num_rdf, run_rdf=lammps_read.rdf_read(filerdf)

#    filecom = filepath+'.com'
#    num_com, run_com=lammps_read.com_read(filecom)

    lammps_plot.log_plots(run_thermo, run_param)

    nrows=2
    fig, axs = plt.subplots(nrows,num_msd,sharey='row',sharex='col',squeeze=False)
    lammps_plot.msd_plots(axs[0,:], run_dens, run_param, 'Density (g/cm^3)')
    lammps_plot.msd_plots(axs[1,:], run_msd, run_param, 'All msd (A^2)')
#    lammps_plot.msd_plots(axs[1,:], run_msd_gRad, run_param, 'gRad msd (A^2)')
#    lammps_plot.msd_plots(axs[2,:], run_msd_gPKA, run_param, 'PKA msd (A^2)')
    
#    lammps_plot.rdf_plots(run_rdf, run_param)
#    lammps_plot.msd2_plots(run_com, msd, run_param)
    plt.show()

    cont = input('Enter 0 to end: ')
    
raise SystemExit()

'''    
    xpka=run_param[0].pos_pka[0]
    ypka=run_param[0].pos_pka[1]
    zpka=run_param[0].pos_pka[2]
    shell_thickness=5
    num_shells=8
    shell=[[] for y in range(num_shells)]

    for i in range(len(run_com[0].Nmolec)):
        xcom=run_com[0].xpos[i]
        ycom=run_com[0].ypos[i]
        zcom=run_com[0].zpos[i]
        dist_from_pka=np.sqrt((xpka-xcom)**2+(ypka-ycom)**2+(zpka-zcom)**2)
        for n in range(num_shells):
            if dist_from_pka > shell_thickness*n and dist_from_pka<=shell_thickness*(n+1):
                shell[n].append(run_com[0].Nmolec[i])

    msd=[[0 for x in range(1)] for y in range(num_shells)]
    for n in range(1,num_com):
        dist=[[] for y in range(num_shells)]
        for i in range(len(run_com[n].Nmolec)):
            xcom=run_com[n].xpos[i]
            ycom=run_com[n].ypos[i]
            zcom=run_com[n].zpos[i]
            xcom_prev=run_com[n-1].xpos[i]
            ycom_prev=run_com[n-1].ypos[i]
            zcom_prev=run_com[n-1].zpos[i]
            xdisp=(xcom-xcom_prev)**2
            ydisp=(ycom-ycom_prev)**2
            zdisp=(zcom-zcom_prev)**2
            tot_disp=xdisp+ydisp+zdisp
            for j in range(num_shells):
                if run_com[n].Nmolec[i] in shell[j]:
                    dist[j].append(tot_disp)

        for j in range(num_shells):
            msd[j].append(np.mean(dist[j]))

'''
