import os

import XPSread
import XPSout

comp_files=[] # array for 'comp' files containing position and FWHM
current_files=[] 
area_files=[] # array for 'area' files containing component areas

element='p2p'
samples='ampcu'
#samples='grblank'

for root, dirs, files in os.walk("."):
    comp_files[:]=[d for d in files if element in d and samples in d and 'txt' in d]
    area_files[:]=[d for d in files if element in d and samples in d and 'txt' in d and 'areas' in d]
    current_files[:]=[d for d in files if element in d and samples in d and 'cur' in d]

print comp_files
print current_files

for d, fn in enumerate(comp_files):
    if current_files==[]:
        cf=None
    else:
        cf=current_files[d]
    fname, peak_names, peak_dat = XPSread.read_comps(fn)

    csv_photo_series = XPSout.csv_series(fname, peak_dat, cf)

    #plot_photo_series = XPSout.plot_series(d, fname, peak_dat)

