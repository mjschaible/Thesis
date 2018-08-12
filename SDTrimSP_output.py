import csv
import os

def appendmeta(log_yld, out_yld):
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

    return out_yld

def write_sputvfluence(out,met_dir):
    for n, met_comp in enumerate(out):
        # Obtain arrays for the meteorite composition elements (+masses and surface binding energies)
        ion_met=met_comp.label[0]
        ion=ion_met.split()[0]+ion_met.split()[1]
        met=ion_met.split()[3]
        comp_elems = met_comp.label[1]
        comp_masses = met_comp.label[2]
        comp_Esbe = met_comp.label[3]
        indv_met = (ion+'_'+met+'.csv')
        Header1 = ('energy', 'fluence', comp_elems, 'TotYld')
        #Header2 = ('', '', 'masses', '')
        #Header3 = ('', '', 'Esbe', '')
        #        return repr((self.label, self.energy, self.fluence, self,Flux, self.totYld))
        #print indv_met
        with open(os.path.join(met_dir,indv_met), "w+") as output:
            output.truncate
            writer = csv.writer(output)
            writer.writerow(Header1)

            for i, flu in enumerate(met_comp.fluence):
                row=[]
                row.append(met_comp.energy[i])
                row.append(met_comp.fluence[i])
                for j,val in enumerate(met_comp.Flux):
                    row.append(met_comp.Flux[j][i])
                row.append(met_comp.totYld[i])
                writer.writerow(row)
    return

def write_indvAvgs(mets,met_dir):
    situation=met_dir.split('/')
    met_class=situation[1]
    ion=situation[-1]
    indv_met = (met_class+'_'+ion+'.csv')
    with open(os.path.join(met_class,indv_met), "w+") as output:
        output.truncate
        #writer = csv.writer(output)
        met_l=list(mets[0][0])
        comp_elems = ','.join(met_l[1])
        comp_masses = met_l[2]
        comp_Esbe = met_l[3]
        Header1 = 'Met Comp,' + comp_elems + ',TotYld,'+ comp_elems + ',TotIYld' + '\n'
        output.write(Header1)

        for n, met in enumerate(mets):
            # Obtain arrays for the meteorite composition elements (+masses and surface binding energies)
            met_l=list(met[0])
            ion_met=met_l[0]
            neut_ylds=','.join(str(val) for val in met[1])
            ion_ylds=','.join(str(val) for val in met[2])
            row = ion_met+','+neut_ylds+','+ion_ylds+'\n'
            #print row
            #print ion_met
            output.write(row)

    return


def write_iavg(sput, lbl=None):

    return

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


def plot_flux(mean_yld,tld_std,met,nf):

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

