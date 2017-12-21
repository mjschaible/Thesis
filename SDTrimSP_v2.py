import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import itertools

import SDTrimSP_readSput
import SDTrimSP_plotSput

targets = ['Mars','Aubrites','LobeO']
#targets = ['Lunar','H-He']

log_yld=[]
out_yld=[]
lbl=[]
tcount=[]
shift = -0.03

for root, dirs, files in os.walk("."):
    dirs[:]=[d for d in dirs if d in targets]
    lcount=0
    for filename in [f for f in files if f.endswith(".log")]:
        fn = os.path.join(root, filename)
        print fn
        log_yld.append(SDTrimSP_readSput.read_logfile(fn))
        lcount+=1
    if lcount>0:
        lbl.append(root.split('/')[1])
        tcount.append(lcount)
        
    for filename in [f for f in files if f.endswith(".out")]:
        fn = os.path.join(root, filename)
        print fn
        out_yld.append(SDTrimSP_readSput.read_outfile(fn))

#tar=[i.label[0].split('->')[1] for i in log_yld]
tar=[i.label[0] for i in log_yld]
#print tar
#print tcount

for i in range(len(out_yld)):
    for j in range(len(out_yld[i])):
        #ctar = out_yld[i][j].label[1].split('->')[1]
        ctar = out_yld[i][j].label[0]
        for n, t in enumerate(tar):
            if t == ctar:
                index=n
        if len(out_yld[i][j].Flux) == len(log_yld[index].label[4]):
            out_yld[i][j].label.append(log_yld[index].label[4]) # append element names
            out_yld[i][j].label.append(log_yld[index].label[5]) # append atomic masses
            out_yld[i][j].label.append(log_yld[index].label[6]) # append SBE
        else:
            print "the number of elements does not match"
            print len(out_yld[i][j].Flux), len(log_yld[nlog].label[4])

color=iter(['red','green','purple','orange','grey','blue','black'])
marker=iter(['D', '^', 'v', '<', 'o', 's',  '>'])

for j, out in reversed(list(enumerate(out_yld))):
    c=next(color)
    mk=next(marker)

    # Determine average yield
    out_yld_avg=SDTrimSP_readSput.average(out)
    nsim = out_yld[0][0].label[0]

    n_out,i_out=SDTrimSP_readSput.yld_comp(out,tar)
    #nE, n_out, i_out=SDTrimSP_readSput.comp_yield(out_yld,tar)

    shift+=0.1
    nf=1
    #ploty_sput=SDTrimSP_plotSput.plot_iavg(sw_out,nf,c,shift,mk)
    nf=2
    mean, std=SDTrimSP_plotSput.plot_iavg(i_out,nf,c,shift,mk,lbl[j])
    nf=3
    #ER_plot=SDTrimSP_plotSput.plot_ER(ER,c,met_class,adir,nf,mk)
    nf=4
    #flux_plot=SDTrimSP_plotSput.plot_flux(mean,std,adir,nf)

plt.show()
#cont = input('Enter 0 to end: ')
