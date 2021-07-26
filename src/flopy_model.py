import numpy as np
import flopy
import flopy.utils.binaryfile as bf
import time

#################################################################################
# Setting up flopy model
#################################################################################

# def flopy_model(modelname,t,wave,Lx=Lx,Ly=Ly,ztop=ztop,zbot=zbot,nlay=nlay,
#                 nrow=nrow,ncol=ncol,delr=delr,delc=delc,delv=delv,botm=botm,
#                 hk=hk,vka=vka,ss=ss,laytyp=laytyp,Plot_length=Plot_length,
#                 headstart=headstart):

def flopy_model(dt,
                wave,
                task_name = 'wave_in_aquifer',
                task_root = './',
                exe_name = "mf2005",
                Lx = 10000,                                # x domain length [m]
                Ly = 10,                                   # y domain length [m]
                ztop = 10,                                 # Top level [m]
                zbot = 0.0,                                # Bottom level [m]
                nlay = 1,                                  # Number of layers
                nrow = 1,                                  # Number of rows
                ncol = 15000,                              # Number of columns
                hk = 25,                                   # Horizontal conductivity
                vka = 1.0,                                 # Vertical conductivity
                ss = 1e-05,                                # Specific storage [m-1]
                laytyp = 0,                                # Unconfined (laytyp > 0)
                headstart = 15,                
                **settings):

    modelname = '{}{}'.format(task_root,task_name)
    # Start timer
    start = time.time()

    # Calculate x series
    xi = np.exp(np.linspace(start=0, stop=np.log(Lx), num=ncol+1)) - 1
    xi[ncol] += 1

    delr = np.diff(xi)  # Voxel length x direction
    delc= Ly / nrow     # Voxel length y direction
    
    Plot_length = np.argmax(delr>2)
    # print(xi[Plot_length])
    
    # extinction = -np.log(0.01/2) / np.sqrt(12*ss/(2*hk))
    # extinction
 
    
    # Define the Stress Periods
    nper   = len(wave)                                  # Number of stress periods
    perlen = dt*np.ones(len(wave))                       # Stress period lengths
    nstp   = np.ones((len(wave)), dtype=np.int32)       # Number of timesteps for each stress period
    steady = np.zeros((len(wave)), dtype=bool)
    
    # Create hydraulic head array
    head_matrix=np.ones((nper, nrow, ncol))
    
    # Create grid 
    mf = flopy.modflow.Modflow(modelname, exe_name=exe_name)
    dis = flopy.modflow.ModflowDis(
        mf,
        nlay,
        nrow,
        ncol,
        delr=delr,     
        delc=delc, 
        top=ztop,
        botm=np.linspace(ztop, zbot, nlay + 1)[1:],
        nper=nper,
        perlen=perlen,
        nstp=nstp,
        steady=steady,
    )
    
    gridx = dis.get_node_coordinates()[1]
    
    # Assign boundary conditions
    strt = headstart * np.ones((nlay, nrow, ncol), dtype=np.float32)
    flopy.modflow.ModflowBas(
        mf, 
        ibound=np.ones((nlay, nrow, ncol), dtype=np.int32), 
        strt=strt,
        )
 
    # Assign hydraulic parameters
    flopy.modflow.ModflowLpf(
        mf, 
        hk=hk, 
        vka=vka, 
        ss=ss, 
        laytyp=laytyp, 
        ipakcb=53,
        )
 
    # Load Modflow solver
    flopy.modflow.ModflowPcg(mf)
    
    # Define stress period data    
    stress_period_data = {}
    
    for i in range(len(wave)):
        stageleft = wave[i] + strt[0,0][0]
        bound_sp = []
        for il in range(nlay):
            condleft = hk * (stageleft - zbot) * delc * 1000
            for ir in range(nrow):
                bound_sp.append([il, ir, 0, stageleft, condleft])
                #bound_sp.append([il, ir, ncol - 1, stageright, condright])
        stress_period_data[i] = bound_sp
    
    # Load stress period data
    flopy.modflow.ModflowGhb(
        mf, 
        stress_period_data=stress_period_data,
        )
    
    # Define output
    stress_period_data = {}
    for kper in range(nper):
        for kstp in range(nstp[kper]):
            stress_period_data[(kper, kstp)] = [
                "save head",
                "save budget",
                "print head",
                "print budget",
            ]
    # Load output file
    flopy.modflow.ModflowOc(
        mf, 
        stress_period_data=stress_period_data, 
        compact=True,
        )
      
    # Write model input files
    mf.write_input()
    
    # Run model
    success, mfoutput = mf.run_model(silent=True, pause=False)
    if not success:
        raise Exception("MODFLOW did not terminate normally.")
       
    # Create the headfile and budget output objects
    headobj = bf.HeadFile(modelname + ".hds") 
    times = headobj.get_times()
    # cbb = bf.CellBudgetFile(modelname + ".cbc")
    head_matrix = np.ones(shape=(nper,Plot_length))
    
    # Extract head files in a readable form
    for i in range(nper):
        head_matrix[i,:] = headobj.get_data(totim=times[i])[0][0,0:Plot_length] - headstart

    # stop timer
    end = time.time()

    # total time taken
    print("Simulation finished successfully")
    print(f" Runtime [s]: {end - start}")
   
    return([times,gridx[0:Plot_length],head_matrix])


