import numpy as np
import matplotlib
matplotlib.use("TkAgg")
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

## --- Define several functions ---
def open_files(path, cur_dir, steplen, run_info):
    pdens = path+'/'+cur_dir+"/*.dens"
    for filedens in glob.glob(pdens): 
        run_dens = lammps_read.data_read(filedens, steplen, run_info)
        
    pmsd = path+'/'+cur_dir+"/*.msd"
    for filemsd in glob.glob(pmsd): 
        run_msd = lammps_read.data_read(filemsd, steplen, run_info)

    prdf = path+'/'+cur_dir+"/*.rdf"
    for filerdf in glob.glob(prdf): 
        run_rdf=lammps_read.rdf_read(filerdf, steplen)

    pcom = path+'/'+cur_dir+"/*.com"
    for filecom in glob.glob(pcom): 
        run_com = lammps_read.com_read(filecom, steplen) # read the center of mass file
        msd_com = lammps_read.find_commsd(run_com, run_info) # find <MSD> for shells around the PKA

    # create DataFrame and save as csv --> converted dump file
    pdump = path+'/'+cur_dir+"/*.dump"
    for filedump in glob.glob(pdump):
        print
        shell = msd_com.descrip[8]
        dataframe, dumpeng=lammps_read.createDataframeFromDump(filedump,shell,run_info,steplen)
        #dataframe.to_csv(filedump+'_conv',sep=' ', index=False)
        # if reading from converted dump file
        #dataframe=createDataframeFromConvDump(sys.argv[1])

    return run_dens, run_msd, run_rdf, msd_com, dumpeng
    
# ------------ Begin main program ----------
def main():
    cont = 1
    num_runs=0
    nf=0
    root = Tk()
    root.withdraw() 
    path = askdirectory()
    # show an "Open" dialog box and return the path to the selected file
    #Tk().withdraw()
    #filename = askopenfilename()
    #filepath = filename.split(".")[0]
    all_dir = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    my_dir = [x for x in all_dir if 'rad00' in x]

    msd_avg=[]
    for cur_dir in my_dir:
        print cur_dir
        logfile = path+'/'+cur_dir+'/log.lammps'
        if 'rad' in cur_dir:
            ftype=1
        else:
            ftype=0
        
        num_runs, run_thermo = lammps_read.log_read(logfile)
        run_info= run_thermo[0].descrip # pull out entire simulation description
        steplen=run_thermo[0].descrip[1] # pull out timestep length [fs]

        run_dens, run_msd, run_rdf, msd_com, dumpeng = open_files(path, cur_dir, steplen, run_info) #

        msd_avg.append(lammps_read.msd_avgf(msd_com))

        if ftype==0:
            nrows=2
            lammps_plot.log_plots(run_thermo, nrows, ftype)

            nrows=1
            #lammps_plot.msd_plots(run_msd, 'Avg All msd',nrows,1,ftype)
            lammps_plot.msd_plots(run_dens, 'Density (g/cm^3)',nrows,1,ftype, nf)

    
        if ftype==1:
            #----- -----
            #lammps_plot.rdf_plots(run_rdf, run_param)
            nrows=2
            lammps_plot.msd_plots(run_msd, 'Avg All msd', nrows, 1,ftype, nf)
            lammps_plot.msd_plots(msd_com, 'Shell Avg. MSD',nrows,2,ftype, nf)

            nrows=2
            lammps_plot.log_plots(run_thermo, nrows, ftype, nf+1)
            lammps_plot.msd_plots(dumpeng, 'Shell Avg. KE',nrows,2,ftype, nf+1)
            
            nf+=2

        plt.show()
        #cont = input('Enter 0 to end: ')
        
    lammps_plot.com_plots(msd_avg)
    cont = input('Enter 0 to end: ')

if __name__ == "__main__":
    main()
