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

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.titlesize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'large'
#mpl.rcParams['axes.labelpad'] = 2.5
mpl.rcParams['xtick.labelsize']='large'
mpl.rcParams['ytick.labelsize']='large'
mpl.rcParams['legend.handlelength']=2.5

## --- Define several functions ---
def open_files(logfile,ftype):
    run_thermo = lammps_read.log_read(logfile)
    run_info= run_thermo[0].descrip # pull out entire simulation description
    steplen=run_thermo[0].descrip[1] # pull out timestep length [fs]
    cpath = '/'.join(logfile.split('/')[0:-1])
    
    pdens = cpath+"/*.dens"
    for filedens in glob.glob(pdens): 
        run_dens = lammps_read.data_read(filedens, steplen, run_info)
        
    pmsd = cpath+"/*.msd"
    for filemsd in glob.glob(pmsd): 
        run_msd = lammps_read.data_read(filemsd, steplen, run_info)

    prdf = cpath+"/*.rdf"
    for filerdf in glob.glob(prdf): 
        run_rdf=lammps_read.rdf_read(filerdf, steplen, run_info)

    pcom = cpath+"/*.com"
    for filecom in glob.glob(pcom): 
        run_com = lammps_read.com_read(filecom, steplen) # read the center of mass file

    # create DataFrame and save as csv --> converted dump file
    pdump = cpath+"/*.dump"
    if ftype == 1:
        msd_com = lammps_read.find_commsd(run_com, run_info) # find <MSD> for shells around the PKA
        for filedump in glob.glob(pdump):
            shell = msd_com.descrip[8]
            dataframe, dumpeng=lammps_read.createDataframeFromDump(filedump,shell,run_info,steplen)
            #dataframe.to_csv(filedump+'_conv',sep=' ', index=False)
            # if reading from converted dump file
            #dataframe=createDataframeFromConvDump(sys.argv[1])
    else:
        msd_com=[]
        dumpeng=[]

    return run_thermo, run_dens, run_msd, run_rdf, msd_com, dumpeng

def plot_files(run_thermo, run_dens, run_msd, run_rdf, msd_com, dumpeng, ftype, nf):
    if nf==1:
        lammps_plot.rdf_plots(run_rdf, nf)
        nf+=1
        #print 'the rdf file is', nf
        
        if ftype==0:
            nrows=2
            lammps_plot.log_plots(run_thermo, nrows, ftype)
            nrows=2
            lammps_plot.msd_plots(run_msd, 'Avg All msd',nrows,ftype)
            lammps_plot.msd_plots(run_dens, 'Density (g/cm^3)',nrows,ftype)

        if ftype==1:
            #----- -----
            nrows=2
            lammps_plot.msd_plots(run_msd, 'System avg. MSD', nrows,ftype, nf)
            lammps_plot.msd_plots(msd_com, 'Shell avg. MSD',nrows,ftype, nf)

            lammps_plot.log_plots(run_thermo, nrows, ftype, nf+1)
            lammps_plot.msd_plots(dumpeng, 'Shell avg. KE',nrows,ftype, nf+1)
            nf+=2

    return nf

# ------------ Begin main program ----------
def main():
    cont = 1
    nf=1
    root = Tk()
    root.withdraw() 
    path = askdirectory()
    # show an "Open" dialog box and return the path to the selected file
    #Tk().withdraw()
    #filename = askopenfilename()
    #filepath = filename.split(".")[0]
    all_dir = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    print all_dir
    if 'rad00' in all_dir[0]:
        my_dir = [x for x in all_dir if 'rad00' in x]
    elif 'EQ' in all_dir[0]:
        my_dir = [x for x in all_dir if 'EQ' in x]
        
    msd_avg=[]
    for cur_dir in my_dir:
        logfile = path+'/'+cur_dir+'/log.lammps'
        if 'rad' in cur_dir:
            ftype=1
        else:
            ftype=0

        run_thermo, run_dens, run_msd, run_rdf, msd_com, dumpeng = open_files(logfile,ftype)
        if ftype == 1:
            #ke_max.append(lamms_read.kemaxf(dumpeng))
            msd_avg.append(lammps_read.msd_avgf(msd_com))
        
        nf = plot_files(run_thermo, run_dens, run_msd, run_rdf, msd_com, dumpeng, ftype, nf)

        #plt.show()
        #cont = input('Enter 0 to end: ')

    if ftype == 1:
        lammps_plot.com_plots(msd_avg)
    cont = input('Enter 0 to end: ')

if __name__ == "__main__":
    main()
