import numpy as np

class SimParam (object):

    def __init__(self,potential,num_molec,timesteps,mol_pka,pos_pka,vel_pka):
        self.potential=potential
        self.num_molec=num_molec
        self.timesteps=timesteps
        self.mol_pka=mol_pka
        self.pos_pka=pos_pka
        self.vel_pka=vel_pka

class LogData (object):

    def __init__(self,thermoCol,step,temp,tempave,pe,peave,ke,keave,etot,press):
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
        if 'timestep' in line:
            timesteps.append(float(column[1]))
        if 'PPPM' in line:
            markers.append(counter)
        if 'Step Temp' in line:
            data_start.append(counter)
            HeadCol = column
        if 'group gPKA molecule' in line:
            mol_pka = column[3]        
        if 'region rRad sphere' in line:
            pos_x = column[3]
            pos_y = column[4]
            pos_z = column[5]
            if pos_z != '${zpos}':
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
    #    runs.append(LogData(potential,step,pe,ke,etot))
        run_param.append(SimParam(potential,num_molec,timesteps[n],mol_pka,pos_pka,vel_pka))
        run_thermo.append(LogData(HeadCol,steparr,temparr,tempavearr,pearr,peavearr,kearr,keavearr,etotarr,pressarr))

    num_runs+=num_cur

    return num_runs, run_param, run_thermo


def data_read(filename):
    counter = 1    
    num_data=0
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
            if len(data_start)>1:
                data_end.append(counter-1)
        elif '#' in line:
            title = line
        
        counter+=1
        if counter == len(lines) and len(data_end) < len(data_start):
            data_end.append(counter)

    num_cur=len(data_start)-num_data
    data_vals=['']*num_cur

    for n in range(num_cur):
        data_vals[n]=lines[data_start[n+num_data]:data_end[n+num_data]-1]
        step=[]
        data_val=[]

        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()
            step.append(vals[0])
            data_val.append(vals[1])
        steparr=np.array(step,dtype='float')
        dataarr=np.array(data_val,dtype='float')
        run_data.append(runData(dataHead,steparr,dataarr))
    num_data+=num_cur

    return num_data, run_data

def rdf_read(filename):
    
    counter=1
    num_rdf=0
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
        elif len(column)==2:
            timestep.append(column[0])
            data_start.append(counter+1)
            if len(data_start)>1:
                data_end.append(counter-1)
        
        counter+=1
        if counter == len(lines) and len(data_end) < len(data_start):
            data_end.append(counter)

    num_cur=len(data_start)-num_rdf
    data_vals=['']*num_cur

    for n in range(num_cur):
        data_vals[n]=lines[data_start[n+num_rdf]:data_end[n+num_rdf]-1]
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
    num_rdf+=num_cur

    return num_rdf, rdf_data

def com_read(filename):
    
    counter=1
    num_com=0
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
            timestep.append(column[0])
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
        Nmolec=[]
        xpos=[]
        ypos=[]
        zpos=[]

        for m in range(len(data_vals[n])):
            vals=data_vals[n][m].split()
            Nmolec.append(vals[0])
            xpos.append(vals[1])
            ypos.append(vals[2])
            zpos.append(vals[3])

        nmolecarr=np.array(Nmolec,dtype='float')
        xposarr=np.array(xpos,dtype='float')
        yposarr=np.array(ypos,dtype='float')
        zposarr=np.array(zpos,dtype='float')
        
        com_data.append(runCOM(float(timestep[n]),nmolecarr,xposarr,yposarr,zposarr))
    num_com+=num_cur

    return num_com, com_data
