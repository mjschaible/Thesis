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
    iontar=sput[0][0].label[0]
    ttar=iontar.replace('->', r'$\rightarrow$')
    ntar=iontar.replace('->', '')
    ion=sput[0][0].label[1].split('->')[0]
    ax1.text(.9,.25,"{}".format(ttar), fontsize=14, horizontalalignment='right', transform=ax1.transAxes)
    
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
