import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
plt.interactive(True)
import glob
import itertools
import sys

from Tkinter import Tk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename

import lammps_read
import lammps_plot

# ------------ Begin main program ----------
def main():
    cont = 1
    num_runs=0

    while cont != 0:
    # show an "Open" dialog box and return the path to the selected file
        Tk().withdraw()
        filename = askopenfilename()
        filepath = filename.split(".")[0]
        print filepath
        if 'rad' in filepath:
            ftype=1
        else:
            ftype=0
        
        num_runs, run_thermo = lammps_read.log_read(filename)
        if run_thermo[0].descrip[2] is not None:
            print 'The excited molecule is {}, the position is {}'.format(run_thermo[0].descrip[4], run_thermo[0].descrip[5])
            print 'The molecule energy is {:.1f}'.format(run_thermo[0].descrip[6])
        else:
            print 'No excited molecule specified', run_thermo[0].descrip[3]

        steplen = run_thermo[0].descrip[1] # pull out timestep length [fs]
        run_info= run_thermo[0].descrip # pull out entire simulation description

        filedens = filepath+'.dens'
        run_dens = lammps_read.data_read(filedens, steplen, run_info)

        filemsd = filepath+'.msd'
        run_msd = lammps_read.data_read(filemsd, steplen, run_info)

        filerdf = filepath+'.rdf'
        run_rdf=lammps_read.rdf_read(filerdf, steplen)

        if ftype==0:
            nrows=2
            lammps_plot.log_plots(run_thermo, nrows, ftype)

            nrows=2
            lammps_plot.msd_plots(run_msd, 'Avg All msd',nrows,1,ftype)
            lammps_plot.msd_plots(run_dens, 'Density (g/cm^3)',nrows,2,ftype)

        if ftype==1:
            #----- -----
            nrows=2
            #lammps_plot.log_plots(run_thermo, nrows, ftype)

            #----- -----
            filecom = filepath+'.com'
            run_com = lammps_read.com_read(filecom, steplen)
            msd_com, shell = lammps_read.find_commsd(run_com, run_info)

            # create DataFrame and save as csv --> converted dump file
            filedump = filepath+'.dump'
            dataframe, dumpeng=lammps_read.createDataframeFromDump(filedump,shell,run_info,steplen)
            #dataframe.to_csv(filedump+'_conv',sep=' ', index=False)
            # comment lines 59 & 60 and use the function call line 63 instead if reading
            # from converted dump file
            #dataframe=createDataframeFromConvDump(sys.argv[1])        

            nrows=2
            lammps_plot.msd_plots(run_msd, 'Avg All msd', nrows, 1,ftype)
            #lammps_plot.msd_plots(run_dens, 'Density (g/cm^3)', nrows, 2)
            lammps_plot.msd_plots(msd_com, 'Shell Avg. MSD',nrows,3,ftype)
            lammps_plot.msd_plots(dumpeng, 'Shell Avg. KE',nrows,4,ftype)
            #    lammps_plot.rdf_plots(run_rdf, run_param)
            plt.show()

        cont = input('Enter 0 to end: ')

if __name__ == "__main__":
    main()
