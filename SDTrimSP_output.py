def plot_log(log,nf,ls,c,lbl=None,h=None,mc=None):
    e_sort = [[] for i in range(len(isbv_pos))]
    yld_sort = [[] for i in range(len(isbv_pos))]
    for i in range(len(log)):
        isbv=log[i].label[1]
        for j in range(len(isbv_pos)):
            match="isbv={}".format(isbv_pos[j])
            if match==isbv:
                e_sort[j].append(log[i].energy)
                yld_sort[j].append(log[i].totYld)

    for j in range(len(isbv_pos)):
        yld_sort[j] = [x for (y,x) in sorted(zip(e_sort[j],yld_sort[j]), key=lambda pair:pair[0])]
        e_sort[j].sort()

    return yld_sort, e_sort

def plot_out1(sput, nf, c, mk,mfc, lbl=None):
    e_sort = [[] for i in range(len(isbv_pos))]
    yld_sort = [[] for i in range(len(isbv_pos))]
    yld_yerr = [[] for i in range(len(isbv_pos))]
    for sputy in sput:
        for i in range(len(sputy)):
            if sputy[i].label[1] != iontar:
                print 'Different ions used'
                exit

            isbv=sputy[i].label[3]

            energy = np.mean(sputy[i].energy)

            totYld = np.mean(sputy[i].totYld[len(sputy[i].totYld)-10:])
            totYld_yerr = np.std(sputy[i].totYld[len(sputy[i].totYld)-10:])
            for j in range(len(isbv_pos)):
                if isbv_pos[j]==isbv:
                    e_sort[j].append(energy)
                    yld_sort[j].append(totYld)
                    yld_yerr[j].append(totYld_yerr)
    for j in range(len(isbv_pos)):
        yld_sort[j] = [x for (y,x) in sorted(zip(e_sort[j],yld_sort[j]), key=lambda pair:pair[0])]
        yld_yerr[j] = [x for (y,x) in sorted(zip(e_sort[j],yld_yerr[j]), key=lambda pair:pair[0])]
        e_sort[j].sort()

    return 


def plot_iavg(sput, nf, ct,shift, mk, lbl=None): 
    for n,i in enumerate(sput):
        for y,elem in enumerate(i.label[4]):
            if not isinstance(i.Flux[y], (int,long)) and y>1 and elem in plotElem:
                ni=plotElem.index(elem)
                dd=len(i.Flux[y])
                navg=10
                elem_iyld[n][ni]=np.mean(i.Flux[y][dd-navg:])
                elem_iyld_var[n][ni]=np.std(i.Flux[y][dd-navg:])
        tot_iyld=np.sum(elem_iyld)

    sw_iavg=np.zeros(len(elem_iyld[0]))
    myarray = np.asarray(elem_iyld)
    metclass_mean=np.mean(myarray, axis=0)
    metclass_std=np.std(myarray, axis=0)
    elem_vsSi_mean=metclass_mean/metclass_mean[nSi]
    #elem_vsSi_std=elem_vsSi_mean*((metclass_std/metclass_mean)**2+(metclass_std[5]/metclass_mean[5])**2)**0.5
    
    print 'Class Averages, {}'.format(t)
    for y, elem in enumerate(plotElem):
        print '{}, {}$\pm${}: vs. Si = {}'.format(elem, metclass_mean[y], metclass_std[y], elem_vsSi_mean[y])

    tot_iyld=np.sum(metclass_mean)
    print 'Total {} yield = {:.3}'.format(t, tot_iyld)

    return


def plot_flux(mean_yld,tld_std,met,nf)

        phiSW=1e8 # solar wind flux, ion/cm^2/s
        P=1./3 # porosity reduction factor
        theta=3.78 # accounting for cosine distribution of incident ion angle
        theta_p=0.68404 # accounting for flux ejected within 20deg from surface normal
        r_sb = 10 # radius of small body, km
        dist = np.arange(5.,100.,5)
        dist_fac=(1-(1-(r_sb/(r_sb+dist))**2)**0.5)

        for i, yld in enumerate(mean_yld):
            c=next(color)
            mk=next(marker)
            surf_flux=phiSW*yld*P*theta*theta_p
            dist_flux=surf_flux*dist_fac
            print 'The {} flux at {} km is {} and at {} km is {}'.format(elem[i], dist[0], dist_flux[0], dist[7], dist_flux[7])

