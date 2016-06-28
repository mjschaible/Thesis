import numpy as np
import pandas as pd
import re
import csv
import sys
import string

# Define a class for the standard lammps log output file
class LogData (object):

    def __init__(self,descrip,thermoCol,step,temp,tempave,pe,peave,ke,keave,etot,press):
        self.descrip=descrip
        self.thermoCol=thermoCol
        self.step=step
        self.temp=temp
        self.tempave=tempave
        self.pe=pe
        self.peave=peave
        self.ke=ke
        self.keave=keave
        self.etot=etot
        self.press=press

# Define a class for mean-squared-displacement(msd) and density files
class runData (object):

    def __init__(self,Head,step,data):
        self.Head=Head
        self.step=step
        self.data=data

# Define a class for radial distribution function files
class runRDF (object):

    def __init__(self,timestep,Nbin,pos,RDF,coordN):
        self.timestep=timestep
        self.Nbin=Nbin
        self.pos=pos
        self.RDF=RDF
        self.coordN=coordN

# Define a class for the center-of-mass file
class runCOM (object):

    def __init__(self,timestep,Nmolec,xpos,ypos,zpos):
        self.timestep=timestep
        self.Nmolec=Nmolec
        self.xpos=xpos
        self.ypos=ypos
        self.zpos=zpos

# Define a class for the dump file
class runDump (object):

    def __init__(self,natoms,timestep,atype,atid,molid,xpos,ypos,zpos,peatom,keatom):
        self.natoms=natoms
        self.timestep=timestep
        self.atype=atype
        self.atid=atid
        self.molid=molid
        self.xpos=xpos
        self.ypos=ypos
        self.zpos=zpos
        self.peatom=peatom
        self.keatom=keatom
        
class prettyfloat(float):
    def __repr__(self):
        return "{:0.2f}".format(self)
    
def log_read(filename):
    counter = 1
    atomline = 0

    mol_pka=0
    pos_pka=0
    vel_pka=0    
    
    data_start=[]
    data_end=[]
    timesteps=[]
    markers=[]
    
    num_runs=0
    run_param = []
    run_thermo = []

    with open(filename,'r') as logfile:
        contents=logfile.read()

    lines=contents.split('\n')

    for line in lines:
        column = line.split()
        if '# create groups ###' in line:
            atomline = counter
        if atomline > 0 and counter == atomline+2:
            num_molec=float(column[0])
            num_atoms=num_molec*3
        if 'pair_style' in line:
            pot_full = line
            potential = column[1]
        if 'timestep' in line and not 'reset' in line and not '$' in line:
            timesteps.append(float(column[1]))
        if 'PPPM' in line:
            markers.append(counter)
        if 'Step Temp' in line:
            data_start.append(counter)
            HeadCol = column
        if 'group gPKA molecule' in line:
            mol_pka = column[3]        
        if 'PKA position' in line:
            pos_x = column[5]
            pos_y = column[6]
            pos_z = column[7]
            if pos_z != '${ypos}':
                pos_x=float(pos_x)
                pos_y=float(pos_y)
                pos_z=float(pos_z)
                pos_pka=[pos_x, pos_y, pos_z]
                pos_pka = map(prettyfloat,pos_pka)
        if 'velocity gPKA' in line:
            vel_x = column[3]
            vel_y = column[4]
            vel_z = column[5]
            if vel_z != '${vz}':
                vel_x=float(vel_x)
                vel_y=float(vel_y)
                vel_z=float(vel_z)
                vel_pka = np.sqrt(vel_x*vel_x+vel_y*vel_y+vel_z*vel_z)
        if 'Loop time' in line:
            data_end.append(counter)
        if len(data_start) > len(timesteps):
            if len(timesteps)==0:
                timesteps.append(0)
            else:
                timesteps.append(timesteps[len(timesteps)-1])
        counter+=1
        if counter == len(lines): 
            if len(data_end) < len(data_start):
                data_end.append(counter)
        
    num_cur=len(data_start)-num_runs

    data_vals=['']*num_cur

    for n in range(num_cur):
        descrip=[]
        data_vals[n]=lines[data_start[n+num_runs]:data_end[n+num_runs]-1]

        step=[]
        temp=[]
        tempave=[]
        pe=[]
        peave=[]
        ke=[]
        keave=[]
        etot=[]
        press=[]

        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()

            step.append(vals[0])
            temp.append(vals[1])
            tempave.append(vals[2])
            pe.append(vals[3])
            peave.append(vals[4])
            ke.append(vals[5])
            keave.append(vals[6])
            etot.append(vals[7])
            press.append(vals[8])
        steparr=np.array(step,dtype='float')
        temparr=np.array(temp,dtype='float')
        tempavearr=np.array(tempave,dtype='float')
        pearr=np.array(pe,dtype='float')
        peavearr=np.array(peave,dtype='float')
        kearr=np.array(ke,dtype='float')
        keavearr=np.array(keave,dtype='float')
        etotarr=np.array(etot,dtype='float')
        pressarr=np.array(press,dtype='float')

        descrip.append(timesteps[n]) # timestep in femtoseconds
        descrip.append('pot={}'.format(potential)) 
        descrip.append(num_molec)
        if 'mol_pka' in locals():
            descrip.append(mol_pka)
            descrip.append(pos_pka)
            descrip.append(float(vel_pka))
        
        run_thermo.append(LogData(descrip,HeadCol,steparr,temparr,tempavearr,pearr,peavearr,kearr,keavearr,etotarr,pressarr))

    num_runs+=num_cur

    return num_runs, run_thermo

def data_read(filename, ts):
    counter = 1    
    data_start=[]
    data_end=[]
    run_data = []

    with open(filename,'r') as logfile:
        contents=logfile.read()

    lines=contents.split('\n')

    for line in lines:
        column = line.split()
        if 'TimeStep' in line:
            dataHead = column
            data_start.append(counter+1)
        elif '#' in line:
            title = line
        
        counter+=1
        if counter == len(lines) and len(data_end) < len(data_start):
            data_end.append(counter)

    num_cur=len(data_start)
    data_vals=['']*num_cur
    
    for n in range(num_cur):
        data_vals[n]=lines[data_start[n]:data_end[n]-1]
        step=[]
        data_val=[]

        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()
            step.append(float(vals[0])*float(ts)/1000) # append timesteps converted to picoseconds
            data_val.append(vals[1]) # append data values (system MSD or density)
        steparr=np.array(step,dtype='float')
        dataarr=np.array(data_val,dtype='float')
        run_data.append(runData(dataHead,steparr,dataarr))

    return run_data

def rdf_read(filename, ts):
    counter=1
    timestep=[]
    data_start=[]
    data_end=[]

    rdf_data=[]

    with open(filename,'r') as logfile:
        contents=logfile.read()

    lines=contents.split('\n')

    for line in lines:
        column=line.split()
        if 'Row' in line:
            Header=column
        elif 'TimeStep' in line:
            Header2=column
        elif len(column)==2: # the timestep + Nbins
            timestep.append(float(column[0])*float(ts)/1000) # Read in timestep and convert to picoseconds
            data_start.append(counter+1)
            if len(data_start)>1:
                data_end.append(counter-1)
        
        counter+=1
        if counter == len(lines) and len(data_end) < len(data_start):
            data_end.append(counter)

    num_cur=len(data_start)
    data_vals=['']*num_cur

    for n in range(num_cur):
        data_vals[n]=lines[data_start[n]:data_end[n]-1]
        nbin=[]
        pos=[]
        rdf=[]
        coordn=[]

        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()
            nbin.append(vals[0])
            pos.append(vals[1])
            rdf.append(vals[2])
            coordn.append(vals[3])

        nbinarr=np.array(nbin,dtype='float')
        posarr=np.array(pos,dtype='float')
        rdfarr=np.array(rdf,dtype='float')
        coordnarr=np.array(coordn,dtype='float')
        
        rdf_data.append(runRDF(timestep,nbinarr,posarr,rdfarr,coordnarr))

    return rdf_data

def com_read(filename, ts):
    
    counter=1
    timestep=[]
    data_start=[]
    data_end=[]

    com_data=[]

    with open(filename,'r') as logfile:
        contents=logfile.read()

    lines=contents.split('\n')

    for line in lines:
        column=line.split()
        if 'Row' in line:
            Header=column
        elif 'TimeStep' in line:
            Header2=column
        elif len(column)==2:
            timestep.append(float(column[0])*float(ts)/1000) # Read in timestep and convert to picoseconds
            data_start.append(counter)
            if len(data_start)>1:
                data_end.append(counter-1)
        
        counter+=1
        if counter == len(lines) and len(data_end) < len(data_start):
            data_end.append(counter-1)

    num_cur=len(data_start)
    data_vals=['']*num_cur

    for n in range(num_cur):
        data_vals[n]=lines[data_start[n]:data_end[n]]
        xpos=[]
        ypos=[]
        zpos=[]
        Nmolec=[]

        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()
            Nmolec.append(vals[0])
            xpos.append(vals[1])
            ypos.append(vals[2])
            zpos.append(vals[3])
        Nmolecarr=np.array(Nmolec,dtype='float')
        xposarr=np.array(xpos,dtype='float')
        yposarr=np.array(ypos,dtype='float')
        zposarr=np.array(zpos,dtype='float')
        
        com_data.append(runCOM(timestep[n],Nmolecarr,xposarr,yposarr,zposarr))

    return com_data

def find_commsd(run_com, descrip):
    # Read the PKA position from the log data file
    xpka=descrip[4][0]
    ypka=descrip[4][1]
    zpka=descrip[4][2]
    
    shell_thickness=5
    num_shells=3
    shell_dat=[]
    shell_dat.append(shell_thickness)
    shell_dat.append(num_shells)
    # This array of arrays will contain the chunk ID for molecules within a shell extending from:
    # dist_from_pka > shell_thickness*n to dist_from_pka<=shell_thickness*(n+1)
    # for n=0 to num_shells
    shell=[[] for y in range(num_shells)] 

    # Define an array that will contain each timestep the COM data is printed out at and
    # the MSD of each molecule with reference to it's initial position at that timestep
    # divided into the shells defined above.
    #->Could add functionality to group all COMs with radius greater than shell_thickness*num_shells
    run_data = []

    # Here I loop though all the chunk ID's and record their position at step=0 
    for i in range(len(run_com[0].Nmolec)):
        xcom_init=run_com[0].xpos[i]
        ycom_init=run_com[0].ypos[i]
        zcom_init=run_com[0].zpos[i]
        #print xcom, ycom, zcom
        # Calculate each molecule's distance from the PKA above
        dist_from_pka=np.sqrt((xpka-xcom_init)**2+(ypka-ycom_init)**2+(zpka-zcom_init)**2)
        for n in range(num_shells):
            # and if that COM is in a given shell it's chunk/molecule ID is added to the appropriate shell array
            if dist_from_pka > shell_thickness*n and dist_from_pka<=shell_thickness*(n+1):
                shell[n].append(run_com[0].Nmolec[i])
    # print the molecule ID's that are within shell_thickness of the PKA
    # print shell[0]

    # Define an array to contain the AVERAGE MSD of all the molecules in a given shell for ALL timesteps
    msd=[[0 for x in range(1)] for y in range(num_shells)]

    # Define an array of ALL the timesteps
    step=[run_com[x].timestep for x in range(len(run_com))]
    steparr=np.array(step,dtype='float')

    # For each timestep
    for n in range(1,len(run_com)):
        # define an array that will have the msd of all molecules in a given shell
        dist=[[] for y in range(num_shells)]
        # loop through all molecules (which should be constant...)
        for i in range(len(run_com[n].Nmolec)):
            # determine the current (x,y,z) position of each molecule
            xcom_cur=run_com[n].xpos[i]
            ycom_cur=run_com[n].ypos[i]
            zcom_cur=run_com[n].zpos[i]
            xcom_init=run_com[0].xpos[i]
            ycom_init=run_com[0].ypos[i]
            zcom_init=run_com[0].zpos[i]
            # calculate the (x,y,z) dispacement from it's initial position
            xdisp=(xcom_cur-xcom_init)**2
            ydisp=(ycom_cur-ycom_init)**2
            zdisp=(zcom_cur-zcom_init)**2
            tot_disp=xdisp+ydisp+zdisp
            # and add the total displacement to it's occording shell
            for j in range(num_shells):
                if run_com[n].Nmolec[i] in shell[j]:
                    dist[j].append(tot_disp)
        # For each cell calculate the mean MSD of all molecules in that shell for this timestep
        for j in range(num_shells):
            msd[j].append(np.mean(dist[j]))

    run_data.append(runData(shell_dat,steparr,msd))
    #print run_data.Head, len(run_data.step), len(run_data.data[0])
    return run_data, shell

def createDataframeFromDump(dumpfile, shell, info, tslen):
    # create dummy headers for 14 columns. 
    # Chis is needed because rows are of unequal length.
    # Missing values will be assigned as NaN
    my_cols=string.ascii_uppercase[0:14]
    dataframe = pd.read_table(dumpfile, delimiter=' ', names=my_cols) #, engine='python')
    # Pull out all row indices that contain 'TIMESTEP' in 2nd column ('B')
    timestep_startrow=dataframe[dataframe.B == 'TIMESTEP'].index
    # Get list of TIMESTEP values
    timesteps = dataframe['A'].loc[timestep_startrow + 1]
    ts=timesteps.values
    ts = [ float(x) * tslen/1000 for x in ts]
    print ts
    # Create Series for TIMESTEPS
    rowsPerTimestep=timestep_startrow[1]-timestep_startrow[0]
    matrix = [ ([i] * rowsPerTimestep) for i in timesteps ]
    timeSeriesData=np.hstack(np.asarray(matrix)) # flatten to 1D array
    # Add TIMESTEP series as new column 
    dataframe['timestep']=pd.Series(timeSeriesData,index=dataframe.index)
    
    # Rename columns
    dataframe.columns=['type','id', 'mol', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'c_pe_atom', 'c_ke_atom', 'L','M','N','timestep']
    
    # Drop all rows that are headers. Headers contain NaN in column 'vx' 
    dataframe=dataframe.dropna(subset=['vx'])
    dataframe=dataframe[dataframe.type != 'ITEM:']
    # Remove additional columns
    dataframe=dataframe.dropna(axis=1)
    # Convert entire dataframe to numeric data type
    dataframe=dataframe.apply(pd.to_numeric)

    run_data=[]
    shellke=[]
    for i in range(len(shell)):
        #print shell[i]
        selectMol=shell[i] # select molecules in shell 'i'
        df=dataframe.loc[dataframe['mol'].isin(selectMol)]
        # aggregate some numbers for data grouped by TIMESTEPS and mol
        groups=df.groupby(['mol', 'timestep'])
        sums=groups['c_pe_atom','c_ke_atom'].sum()
        #print(sums.describe())
        ggroup=df.groupby(['timestep'])
        gsum=ggroup['c_ke_atom'].sum()
        gmean=ggroup['c_ke_atom'].mean()
        shellke.append(gmean.values)

    run_data.append(runData(info,ts,shellke))    
    return dataframe, run_data
    
def createDataframeFromConvDump(convDumpfile):    
    dataframe = pd.read_table(convDumpfile, delimiter=' ')
    return dataframe

def dump_read(filename, ts):
    
    counter=1
    timestep=[]
    natoms=[]
    data_start=[]
    data_end=[]

    nal=0
    dump_data=[]

    with open(filename,'r') as logfile:
        contents=logfile.read()

    lines=contents.split('\n')

    for line in lines:
        column=line.split()
        if 'ITEM: TIMESTEP' in line:
            tsl=counter
            if len(data_start)>0:
                data_end.append(counter-1)
        if counter==tsl+1:
            timestep.append(float(column[0])*float(ts)/1000)
        if 'ITEM: NUMBER OF ATOMS' in line:
            nal=counter
        if counter==nal+1 and nal!=0:
            natoms.append(column[0])
        if 'ITEM: ATOMS' in line:
            data_start.append(counter)
            headers=column[2:]
        counter+=1
        if counter == len(lines) and len(data_end) < len(data_start):
            data_end.append(counter-1)

    num_cur=len(data_start)
    data_vals=['']*num_cur

    for n in range(num_cur):
        data_vals[n]=lines[data_start[n]:data_end[n]]
        atype=[]
        atom_id=[]
        mol_id=[]
        xpos=[]
        ypos=[]
        zpos=[]
        pe_atom=[]
        ke_atom=[]
        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()
            #print vals
            atype.append(vals[0])
            atom_id.append(vals[1])
            mol_id.append(vals[2])
            xpos.append(vals[3])
            ypos.append(vals[4])
            zpos.append(vals[5])
            pe_atom.append(vals[9])
            ke_atom.append(vals[10])
        atypearr=np.array(atype,dtype='float')
        atidarr=np.array(atom_id,dtype='float')
        molidarr=np.array(mol_id,dtype='float')
        xposarr=np.array(xpos,dtype='float')
        yposarr=np.array(ypos,dtype='float')
        zposarr=np.array(zpos,dtype='float')
        pearr=np.array(pe_atom,dtype='float')
        kearr=np.array(ke_atom,dtype='float')
        
        dump_data.append(runDump(natoms[n],timestep[n],atype,atidarr,molidarr,xposarr,yposarr,zposarr,pearr,kearr))

    return dump_data

def find_dumpeng(dump, info):
    xpka=info[4][0]
    ypka=info[4][1]
    zpka=info[4][2]
    shell_thickness=5
    num_shells=3
    shell_dat=[]
    shell_dat.append(shell_thickness)
    shell_dat.append(num_shells)
    shell=[[] for y in range(num_shells)]
    run_data = []
    step=[dump[x].timestep for x in range(len(dump))]
    #print step
    steparr=np.array(step,dtype='float')
    molids=list(set(dump[0].molid)) # Determine a unique set of molecule ID's
    kearr=[]

    for key,value1,xpos,ypos,zpos in zip(dump[0].molid,dump[0].atype,dump[0].xpos,dump[0].ypos,dump[0].zpos):
            if int(value1)==2:
                #print xpos, ypos, zpos
                dist=np.sqrt((xpka-xpos)**2+(ypka-ypos)**2+(zpka-zpos)**2)
                for t in range(num_shells):
                    if dist > shell_thickness*t and dist<=shell_thickness*(t+1):
                        shell[t].append(int(key))
    # print the molecule ID's that are within shell_thickness of the PKA
    #print shell[0]
    # Define an array to contain the AVERAGE MSD of all the molecules in a given shell for ALL timesteps
    keng=[[0 for x in range(1)] for y in range(num_shells)]

    molids=[]
    kemol=[]
    values = defaultdict(int)
    for n in dump: # Loop through each of the output steps in the dump file
        for key, value in zip(n.molid,n.keatom):
            values[key]+=value
        molids_cur, kemol_cur = zip(*sorted(values.items()))
        molids.append(molids_cur)
        kemol.append(kemol_cur)
    #print len(molids[0])
    #print len(kemol[0])

    for n in range(len(dump)-1):
        ke=[[] for y in range(num_shells)]
        for i in range(len(molids[n])):
            for j in range(num_shells):
                if molids[n][i] in shell[j]:
                    ke[j].append(kemol[n][i])
        for j in range(num_shells):
            keng[j].append(np.mean(ke[j]))
    run_data.append(runData(shell_dat,steparr,keng))
    #print run_data[0].Head, len(run_data[0].step), run_data[0].data[0]
    return run_data
