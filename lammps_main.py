import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
plt.interactive(True)
import glob
import os.path
import itertools
import sys

from Tkinter import Tk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename

import lammps_read
import lammps_plot

## --- Control what plots you want to output from here ---

def plot_files(run_thermo, run_dens, run_msd, run_rdf, msd_com, dumpeng_avg, dumpeng_max, ftype):
    #fn = lammps_plot.rdf_plots(run_rdf)
    #print 'the rdf file is', fn

    if ftype==0: # normal type e.g. thermostatting
        nrows=2
        lammps_plot.log_plots(run_thermo, nrows, ftype)
        nrows=2
        lammps_plot.msd_plots(run_msd, 'Avg All msd',nrows,ftype)
        lammps_plot.msd_plots(run_dens, 'Density (g/cm^3)',nrows,ftype)

    if ftype==1: # radiation type = 'rad' in filename
        #lammps_plot.log_plots(run_thermo, ftype)
        #lammps_plot.msd_plots(run_msd, 'System avg. MSD',ftype)
        #lammps_plot.msd_plots(msd_com, 'Shell avg. MSD',ftype)
        #lammps_plot.msd_plots(dumpeng_avg, 'Shell avg. KE',ftype)
        #lammps_plot.msd_plots(dumpeng_max, 'Shell max KE',ftype)
        print 'no plots rendered'
    return

# ------------ Begin main program ----------
def main():
    cont = 1
    nf=1

    outfile=[]
    logfile=[]
    ftype=[]
    dumpfile=[]
    msdfile=[]
    comfile=[]
    densfile=[]
    rdffile=[]
    for root, dirs, files in os.walk("."):
        #print files
        for filename in [f for f in files if f.endswith(".out")]:
            fn_out = os.path.join(root, filename)
            outfile.append(fn_out)
        for filename in [f for f in files if f.endswith(".log")]:
            fn_log = os.path.join(root, filename)
            if 'rad' in fn_log:
                ftype.append(1)
            else:
                ftype.append(0)
            logfile.append(fn_log)
        for filename in [f for f in files if f.endswith(".dump")]:
            fn_dump = os.path.join(root, filename)
            dumpfile.append(fn_dump)
            #print fn_dump
        for filename in [f for f in files if f.endswith(".rdf")]:
            fn_rdf = os.path.join(root, filename)
            rdffile.append(fn_rdf)
        for filename in [f for f in files if f.endswith(".msd")]:
            fn_msd = os.path.join(root, filename)
            msdfile.append(fn_msd)
        for filename in [f for f in files if f.endswith(".com")]:
            fn_com = os.path.join(root, filename)
            comfile.append(fn_com)
        for filename in [f for f in files if f.endswith(".dens")]:
            fn_dens = os.path.join(root, filename)
            densfile.append(fn_dens)

    print 'logfile'
    print logfile
    #print 'ftype'
    #print  ftype
    #print dumpfile
    #print comfile
    #print ftype
    msd_avg=[]
    for n,lf in enumerate(logfile):
        run_thermo = lammps_read.log_read(lf)
        run_info=run_thermo[0].descrip # pull out entire simulation description
        steplen=run_info[1] # pull out timestep length [fs]
        run_dens = lammps_read.data_read(densfile[n], steplen, run_info)
        run_msd = lammps_read.data_read(msdfile[n], steplen, run_info)
        run_rdf=lammps_read.rdf_read(rdffile[n], steplen, run_info)
        run_com = lammps_read.com_read(comfile[n], steplen)
        #print "{} [fs] step length".format(steplen)

        # If the files being considered are radiation files...
        # calculate an array vs. time of the mean squared displacement for shells around the PKA 
        # Additionally, create DataFrame and save as csv --> converted dump file
        if ftype[n] == 1:
            msd_com = lammps_read.find_commsd(run_com, run_info) # find <MSD> for shells around the PKA
            shell = msd_com.descrip[8]
            #print shell
            dataframe, dumpeng_avg, dumpeng_max=lammps_read.createDataframeFromDump(dumpfile[n],shell,run_info,steplen)

            #dataframe.to_csv(filedump+'_conv',sep=' ', index=False)
            # if reading from converted dump file
            #dataframe=createDataframeFromConvDump(sys.argv[1])
        else:
            msd_com=[]
            dumpeng=[]

        # Create function to calculate permanent displacements based on analysis of first and last dump files

        if ftype[n] == 1:
            #ke_max.append(lamms_read.kemaxf(dumpeng))
            # Calculate the average MSD of a shell after steady state has been reached
            msd_avg.append(lammps_read.msd_avgf(msd_com))

        Plot1 = plot_files(run_thermo, run_dens, run_msd, run_rdf, msd_com, dumpeng_avg, dumpeng_max, ftype[n])
        
    if ftype[n] == 1:
        print "Plot the average MSD for the runs analyzed"
        lammps_plot.com_plots(msd_avg)

    cont = input('Enter 0 to end: ')

if __name__ == "__main__":
    main()



'''mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.titlesize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'large'
#mpl.rcParams['axes.labelpad'] = 2.5
mpl.rcParams['xtick.labelsize']='large'
mpl.rcParams['ytick.labelsize']='large'
mpl.rcParams['legend.handlelength']=2.5 '''
