class CompData (object):

    def __init__(self, name,  position, FWHM, area, areamod, atcon, fit):
        self.name=name # 
        self.position=position # 
        self.FWHM=FWHM # 
        self.area=area # 
        self.areamod=areamod #
        self.atcon=atcon #
        self.fit=fit # 

import os
import numpy as np
import matplotlib.pyplot as plt
import itertools

dirs_c=[] # array for 'comp' files containing position and FWHM
dirs_a=[] # array for 'area' files containing component areas
element='o1s'
samples='cramp'
#samples='grblank'
for root, dirs, files in os.walk("."):
    dirs_c[:]=[d for d in files if element in d and samples in d]
    dirs_a[:]=[d for d in files if 'areas' in d and samples in d]

print dirs_c
print dirs_a

for d, fn in enumerate(dirs_c):
    peak_names = []
    peak_dat = []
    plop=0
    
    with open(fn,'r') as logfile:
        contents = logfile.read()
    logfile.close()
    
    lines=contents.split('\n')
    for n, line in enumerate(lines):
        lin = line.strip('\r')
        if n<1:
            col = lin.split('\\')
            fname = col[-1]
            print fname
        elif lin=='':
            plop+=1
        elif 'Name' in lin and plop==1:
            headers = lin.split('\t')
            #print headers
        elif plop<2:
            col = lin.split('\t')
            if col[0] not in peak_names:
                peak_names.append(col[0])
                peak_dat.append([[0] for i in range(len(headers))])
                for i,h in enumerate(col):
                    peak_dat[-1][i][0]=col[i]
                    #print  n, len(peak_dat)-1, peak_dat[-1][i], col[i]
            elif col[0] in peak_names:
                j = peak_names.index(col[0])
                for i,h in enumerate(col):
                    peak_dat[j][i].append(col[i])
                    #print  n, j, i, peak_dat[j][i], col[i]
    print peak_names
    
    color=iter(plt.cm.rainbow(np.linspace(0,1,len(peak_dat))))
    for j, peak in enumerate(peak_dat):
        c=next(color)
        fig = plt.figure(d)
        plt.suptitle(fname)
        for i, head in enumerate(peak):
            if i ==0:
                lbl=head[0]
                print lbl
            if i == 1:
                ax1 = fig.add_subplot(1,4,i)
                ax1.scatter(range(len(peak[i])),peak[i],s=80,color=c)
                ax1.set_title('Position')
            if i == 2:
                ax2 = fig.add_subplot(1,4,i)
                ax2.scatter(range(len(peak[i])),peak[i],s=80,color=c)
                ax2.set_title('FWHM')
            if i == 3:
                ax3 = fig.add_subplot(1,4,i)
                ax3.scatter(range(len(peak[i])),peak[i],s=80,color=c)
                ax3.set_title('Area')
            if i == 5:
                ax4 = fig.add_subplot(1,4,i-1)
                ax4.scatter(range(len(peak[i])),peak[i],s=80,color=c,label=lbl)
                ax4.set_title('%At Conc')
                ax4.legend(bbox_to_anchor=(1,1), loc='upper left', fontsize=10)
            

        
    plt.show()
