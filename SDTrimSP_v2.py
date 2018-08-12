import numpy as np
import os, csv
import itertools
from operator import add

import SDTrimSP_readSput
import SDTrimSP_plot
import SDTrimSP_output
import SDTrimSP_comps

def set_region(cur_dir):
    if 'CC' in cur_dir:
        title = 'Carbonaceous Chondrites'
    elif 'Mars' in cur_dir:
        title = 'Martian Meteorites'

    if 'H-He' in cur_dir:
        label = 'Solar wind'
    elif 'Lobe' in cur_dir:
        label = 'Plasma Lobe'
    elif 'Sheath' in cur_dir:
        label = 'Sheath Region'
    elif 'pSheet' in cur_dir:
        label ='Plasma Sheet'

    region = [title, label]
    return region


targets = ['Mars','LobeO','SheathO', 'pSheetO', 'H-He','LobeO2','SheathO2', 'pSheetO2']
#targets = ['Lunar','H-He']

log_yld=[]
out_yld=[]
cur_dir=[] # Arrary of the directory locations of each *.out file 
tcount=[]
shift = -0.031

for root, dirs, files in os.walk("."):
    dirs[:]=[d for d in dirs if d in targets]
    lcount=0
    for filename in [f for f in files if f.endswith(".out")]:
        fn = os.path.join(root, filename)
        #print fn
        # Create 'output' file containing sputtering data for each meteorite class.
        # The energy, fluence step, nElem x yld, and tot_yld arrays are given for
        # nIncident ions.
        out_yld.append(SDTrimSP_readSput.read_outfile(fn))
        lcount+=1
    if lcount>0:
        # Label should contain the names of the meteorite classes, e.g. "Mars
        for l in range(lcount):
            cur_dir.append(root)
            tcount.append(lcount)

    for filename in [f for f in files if f.endswith(".log")]:
        fn = os.path.join(root, filename)
        #print fn
        # Create independent 'log' files for each simulation file 
        # containing meta data (specific meteorite + ion combination) 
        log_yld.append(SDTrimSP_readSput.read_logfile(fn)) 

print cur_dir
# Append meta data to [out_yld] file for book keeping
out_yld = SDTrimSP_output.appendmeta(log_yld, out_yld)

# Output
#-------------------
# Define color and marker schemes for plotting purposes
# Each color/marker value corresponds to a single meteorite class

class_avgY=[]
class_avgF=[]
for j, out in list(enumerate(out_yld)):
    # Output each meteorite yld vs. fluence data to .csv file
    # *Note that [lbl] contains the meteorite class directory
    # Naming format = Region_composition.csv
    Print_csv = SDTrimSP_output.write_sputvfluence(out, cur_dir[j])

    # Computer the average (steady state) elemental and total yields for each 
    # individual meteorite composition. The n_out are simply the raw fluence and
    # and yields after sorting wrt 'elemlist', while the i_out yields are the 
    # same after correcting by the relative ion yield fractions.
    # 'outfile' contains a single line for the elemental averages of each
    # meteorite composition. 
    # naming format = Class_region.csv
    n_out,i_out,IndvAvgs=SDTrimSP_comps.yld_comp(out)
    Print_csv = SDTrimSP_output.write_indvAvgs(IndvAvgs, cur_dir[j])

    # Compute the mean SI YIELD per meteorite class for a given irradiation condition
    #nE, n_out, i_out=SDTrimSP_readSput.comp_yield(out_yld,tar)
    # yields = neutY,neutY_std,ionY,ionY_std
    yields,elems=SDTrimSP_comps.comp_iavg(IndvAvgs,cur_dir[j])
    yld_info=[cur_dir[j],elems,yields[0],yields[1],yields[2],yields[3]]
    class_avgY.append(yld_info)

    # Compute the mean SI FLUX per meteorite class
    # fluxes = neutF,neutF_std,ionF,ionF_std
    fluxes = SDTrimSP_comps.comp_flux(yields, cur_dir[j])
    flux_info=[cur_dir[j],elems,fluxes[0],fluxes[1],fluxes[2],fluxes[3]]
    class_avgF.append(flux_info)

## I don't totally understand what is going on in the next couple functions
# Determine solar wind fluxes for 95% H and 5% He
cf = iter(enumerate(class_avgF))
fn_next=0
fn_nnext=0
for j, flux in cf:
    fn_next=0
    fn_nnext=0
    fn=class_avgF[j][0]
    if j<len(class_avgF)-1:
        fn_next=class_avgF[j+1][0]
    if j<len(class_avgF)-2:
        fn_nnext=class_avgF[j+2][0]
    if j<len(class_avgF)-1 and fn == fn_next and 'H' in fn:
        sw_fluxes=[class_avgF[j],class_avgF[j+1]]
        for flux in sw_fluxes:
            if 'H' == flux[1][0]:
                H_flux = [0.95*f for i,f in enumerate(flux) if i > 1]
            elif 'He'==flux[1][0]:
                He_flux = [0.05*f for i,f in enumerate(flux) if i > 1]
        sw_f=H_flux+He_flux
        flux[1][0]='SW'
        region = set_region(flux[0])
        sw_flux=[region,flux[1],sw_f[0],sw_f[1],sw_f[2],sw_f[3]]
        next(itertools.islice(cf,1,1),None)

pos=next(i for i,d in enumerate(cur_dir) if 'H-He' in d)
class_avgF[pos]=sw_flux
del class_avgF[pos+1]

# Determine surface fluxes in magnetosphere regions O, O2, (and H can be neglected due to low fluxes)
cf = iter(enumerate(class_avgF))
Lobe=[inc_cond for d,inc_cond in cf if 'Lobe' in inc_cond[0]]
region=set_region(Lobe[0][0])
Lobe=[region,Lobe[0][1],Lobe[0][2]+Lobe[1][2],Lobe[0][3]+Lobe[1][3],Lobe[0][4]+Lobe[1][4],Lobe[0][5]+Lobe[1][5]]

cf = iter(enumerate(class_avgF))
Sheath=[inc_cond for d,inc_cond in cf if 'Sheath' in inc_cond[0]]
region=set_region(Sheath[0][0])
Sheath=[region,Sheath[0][1],Sheath[0][2]+Sheath[1][2],Sheath[0][3]+Sheath[1][3],Sheath[0][4]+Sheath[1][4],Sheath[0][5]+Sheath[1][5]]

cf = iter(enumerate(class_avgF))
pSheet=[inc_cond for d,inc_cond in cf if 'pSheet' in inc_cond[0]]
region=set_region(pSheet[0][0])
pSheet=[region,pSheet[0][1],pSheet[0][2]+pSheet[1][2],pSheet[0][3:]+pSheet[1][3],pSheet[0][4]+pSheet[1][4],pSheet[0][5]+pSheet[1][5]]


class_avgF=[sw_flux,Lobe,Sheath,pSheet]
# Plotting
# Plot the elemental ratios based on yield calculations (these shouldn't change based on radiation)
nf=1
#Print_SiRatio=SDTrimSP_plot.plot_siratio(class_avgY, nf)

# For Neutral Yields
nf=1
#Plot_ClassAvgs = SDTrimSP_plot.plot_clavg(elems,neut_avgs,cur_dir[j])
#ploty_sput=SDTrimSP_plotSput.plot_iavg(n_out,nf,c,shift,mk)

# For Ion Yields
nf=2
#Print_csv=SDTrimSP_output.write_iavg(i_out,nf,c,shift,mk,cur_dir[j])
#Plot_ClassAvgs = SDTrimSP_plot.plot_clavg(class_avgY, nf)

# For Ion Fluxes
nf=3
Plot_ClassAvgs = SDTrimSP_plot.plot_clavg(class_avgF, nf)

# For Neut Fluxes
nf=4
#Plot_ClassAvgs = SDTrimSP_plot.plot_clavg(class_avgF, nf)

#ER_plot=SDTrimSP_plotSput.plot_ER(ER,c,met_class,adir,nf,mk)
nf=4
#flux_plot=SDTrimSP_plotSput.plot_flux(mean,std,adir,nf)


#cont = input('Enter 0 to end: ')
