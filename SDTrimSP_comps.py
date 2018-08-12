import numpy as np
from operator import add
from SDTrimSP_readSput import LogData

elemlist = ['H', 'He', 'C_g','Na','Mg','Al','Si','S','K','Ca','Ti','Mn','Fe','Ni','O']
ifrac = 0.01
relcorr = [0.00,0.00,0.01,42.13,8.80,3.85,1.000,0.01,229.488,7.00,5.71,8.81,1.78,2.886,0.01]
navg=100

def find_sputvar(sput, i):
    sput_data=[]
    for j in range(len(sput)):
        iyld=[]
        fs=sput[j].fluence[i]

        for k in range(len(sput[j].Flux)):
            iyld.append(sput[j].Flux[k][i])
        totYld=sput[j].totYld[i]

        sput_data.append(LogData(sput[j].label, sput[j].energy, fs, iyld, totYld))

    return sput_data

def yld_comp(out_yld, met_class=None):
    # Array of neut and ion yields vs. fluence for all elements
    n_met = len(out_yld)
    nfyld=[[0]*len(elemlist) for i in range(n_met)] 
    ifyld=[[0]*len(elemlist) for i in range(n_met)]
    # Arrays for returned LogData objects
    tot_nyld=[] 
    tot_iyld=[]
    dd = [] 
    navg=100

    k=-1
    # out_yld = files
    # cycle over ion/target combinations
    outfile=[]
    for i in range(len(out_yld)):
        outline=[]
        outline.append(out_yld[i].label)
        comp_lbl=[]
        ion = out_yld[i].label[0].split('->')[0]
        target = out_yld[i].label[0].split('->')[1]

        # Take the average of the last 'navg' values to get steady state mean
        tyld=np.mean(out_yld[i].totYld[len(out_yld[i].totYld)-navg:])
        dd.append(len(out_yld[i].totYld)) # number of fluence steps

        k+=1
        for y, elem in enumerate(elemlist):
            # Iterate over the elements in the composition and organize according to elemlist
            match = [l for l, x in enumerate(out_yld[i].label[1]) if x == elem]
            if match:
                yld=out_yld[i].Flux[match[0]]
                #print k, y, elem, yld
                nfyld[k][y]=yld
                ifyld[k][y]=yld*relcorr[y]*ifrac

        # Array of average neut and ion yields at steady state for all elements
        avg_nyld=[]
        avg_iyld=[]
        for y, elem in enumerate(elemlist):
            if isinstance(nfyld[k][y], (int,long)):
                pass
            else:
                avg=np.mean(nfyld[k][y][dd[k]-navg:])
                if avg>1: avg-=1
                avg_nyld.append(avg)
                avg_iyld.append(avg*relcorr[y]*ifrac)
                comp_lbl.append(elem)
                #print '{}: Y^elem={:.5f}, Y^i={:.6f}'.format(elem,avg_nyld[k][y],avg_iyld[k][y])
        outline[0][1]=comp_lbl
        ttot=np.sum(avg_nyld)
        itot=np.sum(avg_iyld)
        avg_nyld.append(ttot)
        outline.append(avg_nyld)
        avg_iyld.append(itot)
        outline.append(avg_iyld)
        
        outfile.append(outline)
        #print '{},{}: Y^tot={:.4f}, Y^i={:.6f}'.format(ion,target,ttot,itot)

        tyld=np.zeros(len(nfyld[k][-1]))
        iyld=np.zeros(len(ifyld[k][-1]))
        # Determine the total yield as a function of fluence
        # len(tyld) = len(fluence)
        for y, elem in enumerate(elemlist):
            if not isinstance(nfyld[k][y],(int,long)) and y>1:
                tyld=map(add,nfyld[k][y],tyld) 
                iyld=map(add,ifyld[k][y],iyld)

        lbl=out_yld[i].label[:]
        lbl[3]=elemlist
        tot_nyld.append(LogData(lbl,out_yld[i].energy,out_yld[i].fluence,nfyld[k],tyld))
        tot_iyld.append(LogData(lbl,out_yld[i].energy,out_yld[i].fluence,ifyld[k],iyld))
    return tot_nyld, tot_iyld, outfile

def comp_iavg(IndvAvgs,cur_dir):
    elems=IndvAvgs[0][0][1]
    elem_yld=[[]*len(IndvAvgs[0][0][1]) for i in xrange(len(IndvAvgs[0][0][1]))]
    i_yld=[[]*len(IndvAvgs[0][0][1]) for i in xrange(len(IndvAvgs[0][0][1]))]
    for n, indvMet in enumerate(IndvAvgs):
        for y,elem in enumerate(indvMet[0][1]):
            elem_yld[y].append(indvMet[1][y])
            i_yld[y].append(indvMet[2][y])

    neut_avgs=np.average(elem_yld, axis=1)
    neut_std=np.std(elem_yld, axis=1)
    ion_avgs=np.average(i_yld, axis=1)
    ion_std=np.std(i_yld, axis=1)
    yields = [neut_avgs,neut_std,ion_avgs,ion_std]
    return yields,elems

def comp_flux(yields, cd):
    #[SW Flux, Lobe Flux, Sheath Flux, pSheet Flux], ion/cm^2/s
    fluxes = [ 1e8, 1.21e7, 8.37e5, 1.11e7]
    P=1./3 # porosity reduction factor
    theta=3.78 # accounting for cosine distribution of incident ion angle
    theta_p=0.68404 # accounting for flux ejected within 20deg from surface normal
    surf_flux=[]

    if 'H' in cd:
        Phi=fluxes[0]
    elif 'Lobe' in cd:
        Phi=fluxes[1]
    elif 'Sheath' in cd:
        Phi=fluxes[2]
    elif 'pSheet' in cd:
        Phi=fluxes[3]

    for Y in yields:
        surf_flux.append(Phi*Y*P*theta*theta_p)

    return surf_flux

def f_v_dist():
    return

def comp_yield(out_yld, tar):
    hyld=[[0]*len(elemlist) for i in range(len(tar))]
    h_iyld=[[0]*len(elemlist) for i in range(len(tar))]
    heyld=[[0]*len(elemlist) for i in range(len(tar))]
    he_iyld=[[0]*len(elemlist) for i in range(len(tar))]
    sw_avg=[[0]*len(elemlist) for i in range(len(tar))]
    sw_iavg=[[0]*len(elemlist) for i in range(len(tar))]
    sw_yld=[]
    sw_iyld=[]
    totyld=[]
    totiyld=[]
    swtot_yld = []
    swion_yld = []
    dd = []
    Elem_Ratio=[[0]*4 for i in range(len(tar))]
    #elemlist = ['H', 'He', 'C_g','Na','Mg','Al','Si','S','K','Ca','Ti','Mn','Fe','Ni','O']
    nE=[4, 5, 9, 12]

    for i in range(len(out_yld)): # cycle over files
        for j in range(len(out_yld[i])): # cycle over ion/target combinations
            ion = out_yld[i][j].label[0].split()[1]
            target = out_yld[i][j].label[0].split()[3]
            #print ion, target
            tyld=np.mean(out_yld[i][j].totYld[len(out_yld[i][j].totYld)-navg:])
            dd.append(len(out_yld[i][j].totYld))
            
            for k in range(len(tar)):
                if ion == 'H' and target in tar[k]:
                    #print 'The {} total yield is {}'.format(out_yld[i][j].label[0], tyld)
                    for y, elem in enumerate(elemlist):
                        match = [l for l, x in enumerate(out_yld[i][j].label[3]) if x == elem]
                        if match:
                            yld=0.95*out_yld[i][j].Flux[match[0]]
                            hyld[k][y]=yld
                            h_iyld[k][y]=yld*relcorr[y]*ifrac
    
                elif ion == 'He' and target in tar[k]:
                    #print 'The {} total yield is {}'.format(out_yld[i][j].label[0], tyld)
                    for y, elem in enumerate(elemlist):
                        match = [l for l, x in enumerate(out_yld[i][j].label[3]) if x == elem]
                        if match:
                            yld=0.05*out_yld[i][j].Flux[match[0]]
                            heyld[k][y]=yld
                            he_iyld[k][y]=yld*relcorr[y]*ifrac

    for k in range(len(tar)):
        print tar[k]
        for y, elem in enumerate(elemlist):
            if isinstance(hyld[k][y], (int,long)) or isinstance(heyld[k][y], (int,long)):
                if y<2:
                    havg=np.mean(hyld[k][0][dd[k]-navg:])
                    heavg=np.mean(heyld[k][1][dd[k]-navg:])
                else:
                    havg=0
                    heavg=0
                #print 'Y{}={} and Y{}={}'.format(elem, havg, elem, heavg)
            else:
                #print elemlist[match[0]], elem
                #print match[0]
                #print hyld[k][y]
                havg=np.mean(hyld[k][y][dd[k]-navg:])
                heavg=np.mean(heyld[k][y][dd[k]-navg:])
                sw_avg[k][y]=havg+heavg
                sw_iavg[k][y]=sw_avg[k][y]*relcorr[y]*ifrac
                print 'Y^tot({})={:.5f}, Y^i({})={:.6f}'.format(elem,sw_avg[k][y],elem,sw_iavg[k][y])
        ttot=np.sum(sw_avg[k][2:])
        itot=np.sum(sw_iavg[k][2:])
        print 'Y^tot({})={:.4f}, Y^i_({})={:.6f}'.format(tar[k],ttot ,tar[k],itot)
        
        for l in range(len(nE)):
            Elem_Ratio[k][l]=sw_iavg[k][nE[l]]/sw_iavg[k][6]

        sw_yld.append(map(add,hyld[k],heyld[k]))
        sw_iyld.append(map(add,h_iyld[k],he_iyld[k]))
        swtyld=np.zeros(len(sw_iyld[k][0]))
        swiyld=np.zeros(len(sw_iyld[k][0]))
        # Determine the total SW yield as a function of fluence
        # len(swtyld) = len(fluence)
        for y, elem in enumerate(elemlist):
            if not isinstance(sw_iyld[k][y], (int,long)) and y>1:
                swtyld=map(add,sw_yld[k][y],swtyld) 
                swiyld=map(add,sw_iyld[k][y],swiyld)
        totyld.append(swtyld)
        totiyld.append(swiyld)
        
    for k, ttar in enumerate(tar):
        tot_iyld=0
        tot_yld=0
        sw_iavg_c=0
        for y,elem in enumerate(elemlist):
            if not isinstance(sw_iyld[k][y], (int,long)) and y>1:
                havg=np.mean(hyld[k][y][dd[k]-navg:])
                heavg=np.mean(heyld[k][y][dd[k]-navg:])
                sw_iavg_c+=(havg+heavg)*relcorr[y]*ifrac

                elem_yld= np.mean(sw_yld[k][y][dd[k]-navg:])
                elem_iyld = np.mean(sw_iyld[k][y][dd[k]-navg:])
                
                tot_iyld+=elem_iyld
                tot_yld+=elem_yld
        #print ttar
        #print tot_yld, tot_iyld, sw_iavg_c
        sw_yld_lbl=out_yld[i][k].label[:]
        sw_yld_lbl[4]=elemlist
        asdf=sw_yld_lbl[1].replace('He','SW')
        sw_yld_lbl[1]=asdf
        asdf=sw_yld_lbl[1].replace('H','SW')
        sw_yld_lbl[1]=asdf
        swtot_yld.append(LogData(sw_yld_lbl,out_yld[i][k].energy,out_yld[i][k].fluence,sw_yld[k],totyld[k]))
        swion_yld.append(LogData(sw_yld_lbl,out_yld[i][k].energy,out_yld[i][k].fluence,sw_iyld[k],totiyld[k]))
        
    return Elem_Ratio, swtot_yld, swion_yld

def log_comp(yld, tar):
    hyld=[[] for i in range(len(tar))]
    h_iyld=[[] for i in range(len(tar))]
    heyld=[[] for i in range(len(tar))]
    he_iyld=[[] for i in range(len(tar))]
    sw_yld=[]
    sw_iyld=[]
    Elem_Ratio=[[0]*4 for i in range(len(tar))]
    nE=[4, 5, 9, 12]
    #print out_yld[i].label[4]
    #print out_yld[i].Flux
    #print tar
    neut_data=[]
    ion_data=[]

    for i in range(len(yld)):
        for k in range(len(tar)):
            ion = yld[i].label[0].split('->')[0]
            target = yld[i].label[0].split('->')[1]
            tyld=yld[i].totYld
            if ion == 'H' and target == tar[k]:
                #print 'The {} yield is {}'.format(yld[i].label[0], tyld)
                #print yld[i].label[4]
                for y, elem in enumerate(elemlist):
                    match = [l for l, x in enumerate(yld[i].label[4]) if x == elem]
                    if match:
                        ny=yld[i].Flux[match[0]]
                        hyld[k].append(ny)
                        #print elem, yld[i].label[4][match[0]], relcorr[y], ny
                        h_iyld[k].append(ny*relcorr[y]*ifrac)

            elif ion == 'He' and target == tar[k]:
                #print 'The {} yield is {}'.format(yld[i].label[0], tyld)
                for y, elem in enumerate(elemlist):
                    match = [l for l, x in enumerate(yld[i].label[4]) if x == elem]
                    if match:
                        cy=yld[i].Flux[match[0]]
                        heyld[k].append(cy)
                        he_iyld[k].append(cy*relcorr[y]*ifrac)

    for k in range(len(tar)):
        sw_yld.append(map(add,hyld[k],heyld[k]))
        sw_iyld.append(map(add,h_iyld[k],he_iyld[k]))
        neut_data.append(LogData(yld[k].label,yld[k].energy,yld[k].fluence,sw_yld[k],np.sum(sw_yld[k])))
        ion_data.append(LogData(yld[k].label,yld[k].energy,yld[k].fluence,sw_iyld[k],np.sum(sw_iyld[k])))
        npR=0
        if npR==1:
            print 'H Yield'
            for l in range(len(yld[0].label[4])):
                #print '{}, {:.4f}, {:.6f}'.format(out_yld[i][j].label[4][l],sw_yld[k][l],sw_iyld[k][l])
                print '{} {:.2e}'.format(yld[0].label[4][l], hyld[k][l])
            print 'He Yield'
            for l in range(len(yld[i].label[4])):
                print '{} {:.2e}'.format(yld[0].label[4][l], heyld[k][l])
        #print 'The {} total SW yield is {:.3f}'.format(tar[k], np.sum(sw_yld[k]))
        #print 'The {} total SW ion yield is {:.5f}'.format(tar[k], np.sum(sw_iyld[k]))

    return neut_data, ion_data
                        
