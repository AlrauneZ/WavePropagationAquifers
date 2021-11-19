import numpy as np
import flopy
import flopy.utils.binaryfile as bf
import time

#################################################################################
# Setting up flopy model
#################################################################################

def model_confined(dt,
                wave,
                task_name = 'wave_in_aquifer_confined',
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
#                vka = 25,                                  # Vertical conductivity
                ss = 5e-05,                                # Specific storage [m-1]
                dx_max = 2,
                **settings):

    modelname = '{}{}'.format(task_root,task_name)
    # Start timer
    start = time.time()

    # Calculate x series
    xi= np.logspace(0,np.log10(Lx),ncol+1)-1
    # xi = np.exp(np.linspace(start=0, stop=np.log(Lx), num=ncol+1)) - 1
    xi[ncol] += 1

    delr = np.diff(xi)  # Voxel length x direction
    delc= Ly / nrow     # Voxel length y direction
    
    ix_max = np.argmax(delr>dx_max)
    # print(xi[ix_max])
    
    # Define the Stress Periods
    nper   = len(wave)                             # Number of stress periods 
   
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
        perlen=dt*np.ones(nper),                     # Stress period lengths,
        nstp= np.ones((nper), dtype=np.int32),       # Number of timesteps for each stress period
        steady= np.zeros((nper), dtype=bool),
    )
    
    gridx = dis.get_node_coordinates()[1]
    headstart = ztop - min(wave-np.mean(wave)) + 1
    
    # Assign boundary conditions
    flopy.modflow.ModflowBas(
        mf, 
        ibound=np.ones((nlay, nrow, ncol), dtype=np.int32), 
        strt=headstart * np.ones((nlay, nrow, ncol), dtype=np.float32),
        )
 
    # Assign hydraulic parameters
    flopy.modflow.ModflowLpf(
        mf, 
        hk=hk, 
        vka=hk, 
        ss=ss, 
        ipakcb=53,
        )
 
    # Load Modflow solver
    flopy.modflow.ModflowPcg(mf)
    
    # Define stress period data    
    stress_period_data = {}
    for i in range(nper):
        stageleft = wave[i] + headstart 
        # condleft = hk * (stageleft - zbot) * delc * 1000
        condleft = 10000	
        stress_period_data[i] = [0, 0, 0, stageleft, condleft]
       
    # Load stress period data
    flopy.modflow.ModflowGhb(
        mf, 
        stress_period_data=stress_period_data,
        )
    
    # Define output
    stress_period_out = {}
    for kper in range(nper):
        stress_period_out[(kper, 0)] = [
            "save head",
            "print head",
        ]
    
    # Load output file
    flopy.modflow.ModflowOc(
        mf, 
        stress_period_data=stress_period_out, 
        compact=True,
        )
      
    # Write model input files
    mf.write_input()
    
    # Run model
    success, mfoutput = mf.run_model(silent=True, pause=False)
    if not success:
        raise Exception("MODFLOW did not terminate normally.")
       
    # Create the headfile output objects
    headobj = bf.HeadFile(modelname + ".hds") 
    times = headobj.get_times() # end time points of stress periods
   
    # Create hydraulic head array
    head_matrix = np.ones(shape=(nper,ix_max))
    
    # Extract head files in a readable form
    for i in range(nper):
        head_matrix[i,:] = headobj.get_data(totim=times[i])[0][0,0:ix_max] - headstart

    # stop timer
    end = time.time()

    # total time taken
    print("Simulation finished successfully")
    print(f" Runtime [s]: {end - start}")
   
    return([times,gridx[0:ix_max],head_matrix])

def model_leakage(dt,
                wave,
                task_name = 'wave_in_aquifer_leakage',
                task_root = './',
                exe_name = "mf2005",
                Lx=10000,
                Ly=10,
                zbot=0.0,
                nlay=3,
                nrow=1,
                ncol=15000,
                thickness_lay1 = 10,
                thickness_lay2 = 1, # need to be one to be consistant with the leakage coefficient
                thickness_lay3 = 10,
                hk_unconfined = 25,
                hk_confining = 0.01,
                hk_confined = 25,
#                vka=25,
                ss_confined = 5e-5,
                ss_unconfined = 0.025,
                fluctuation_unconfined_aquifer=False, ### regulating BC in unconfined aquifer
                dx_max = 2,
                **settings):
    
    modelname = '{}{}'.format(task_root,task_name)
    # Start timer
    start = time.time()
    
    # Calculate x series
    xi= np.logspace(0,np.log10(Lx),ncol+1)-1
    # xi = np.exp(np.linspace(start=0, stop=np.log(Lx), num=ncol+1)) - 1
    xi[ncol] += 1
    
    delr = np.diff(xi)  # Voxel length x direction
    delc= Ly / nrow     # Voxel length y direction

    ix_max = np.argmax(delr>dx_max)
    
    # Define the Stress Periods

    nper   = len(wave)                                     
    ztop = thickness_lay1 + thickness_lay2 + thickness_lay3

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
        botm=np.array([thickness_lay3+thickness_lay2,thickness_lay3,0.0]),
        nper=nper,
        perlen=dt*np.ones(nper),                 # Number of stress periods
        nstp=np.ones((nper), dtype=np.int32),   # Number of timesteps for each stress period
        steady = np.zeros((nper), dtype=bool),  # Stress period lengths
    )
    
    gridx = dis.get_node_coordinates()[1]
    
    # Assign boundary conditions
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    if fluctuation_unconfined_aquifer is False:  
        ibound[0,:,:] = -1

    headstart = ztop - min(wave-np.mean(wave)) + 1
       
    flopy.modflow.ModflowBas(mf, 
                             ibound=ibound, 
                             strt = headstart * np.ones((nlay, nrow, ncol), dtype=np.float32),
                             )
    
    # Assign hydraulic conductivities
    hk = np.ones((nlay,nrow,ncol))
    hk[0,:,:] = hk_unconfined
    hk[1,:,:] = hk_confining
    hk[2,:,:] = hk_confined
    
    vka = hk
       
    ss = np.ones((nlay,nrow,ncol))
    ss[0,:,:] = ss_unconfined
    ss[1,:,:] = ss_confined
    ss[2,:,:] = ss_confined
    # ss = ss_confined * np.ones((nlay,nrow,ncol))
       
    flopy.modflow.ModflowLpf(mf, 
                             hk=hk, 
                             vka=vka, 
                             ss=ss, 
                             # laytyp=laytyp, 
                             ipakcb=53
                             )
    
    # Load Modflow solver
    flopy.modflow.ModflowPcg(mf)
    
    # Define stress period data    
    stress_period_data = {}
    
    for i in range(len(wave)):
        stageleft = wave[i] + headstart
        bound_sp = []
        for il in range(nlay):
            bound_sp.append([il, 0, 0, stageleft, 10000])
        stress_period_data[i] = bound_sp
    
    # Load stress period data
    flopy.modflow.ModflowGhb(mf, 
                             stress_period_data=stress_period_data,
                             )
    
    # Define output
    stress_period_out = {}
    for kper in range(nper):
        stress_period_out[(kper, 0)] = [
            "save head",
            "print head",
        ]
    
    # Load output file
    flopy.modflow.ModflowOc(mf, 
                            stress_period_data=stress_period_out, 
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
    
    # Create hydraulic head array
    head_matrix = np.ones(shape=(nper,ix_max))
    #head_matrix2 = np.ones(shape=(nper,ix_max))
    
    # Extract head files in a readable form
    for i in range(nper):
        head_matrix[i,:] = headobj.get_data(totim=times[i])[2][0,0:ix_max] - headstart # results in confined aquifer (3rd layer)
        #head_matrix2[i,:] = headobj.get_data(totim=times[i])[0][0,0:ix_max] - headstart    # results in unconfined aquifer
    
    # stop timer
    end = time.time()
    
    # total time taken
    print("Simulation finished successfully")
    print(f" Runtime [s]: {end - start}")
    
    return([times,gridx[0:ix_max],head_matrix])

def model_leakage_old(dt,
                wave,
                task_name = 'wave_in_aquifer_leakage',
                task_root = './',
                exe_name = "mf2005",
                Lx=10000,
                Ly=10,
                zbot=0.0,
                nlay=3,
                nrow=1,
                ncol=15000,
                thickness_lay1 = 10,
                thickness_lay2 = 1,
                thickness_lay3 = 9,
                hk_unconfined=25,
                hk_confining=0.01,
                hk_confined=25,
#                vka=25,
                ss_confined=5E-5,
                lay1_unconfined = True,
                sy=0.25,
                laytyp=0,
                headstart=25,
                fluctuation_1st_aquifer=False,
                **settings):
    
    modelname = '{}{}'.format(task_root,task_name)
    # Start timer
    start = time.time()
    
    # Calculate x series
    xi = np.exp(np.linspace(start=0, stop=np.log(Lx), num=ncol+1)) - 1
    xi[ncol] += 1
    
    delr = np.diff(xi)  # Voxel length x direction
    delc= Ly / nrow     # Voxel length y direction
    ztop = thickness_lay1 + thickness_lay2 + thickness_lay3

    ix_max = np.argmax(delr>2)
    
    # Define the Stress Periods
    nper   = len(wave)                                  # Number of stress periods
    perlen = dt*np.ones(len(wave))                       # Stress period lengths
    nstp   = np.ones((len(wave)), dtype=np.int32)       # Number of timesteps for each stress period
    steady = np.zeros((len(wave)), dtype=bool)
    
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
        botm=np.array([thickness_lay3+thickness_lay2,thickness_lay3,0.0]),
        nper=nper,
        perlen=perlen,
        nstp=nstp,
        steady=steady,
    )
    
    gridx = dis.get_node_coordinates()[1]
    
    # Assign boundary conditions
    if fluctuation_1st_aquifer:  
        ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    else:
        ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
        ibound[0,:,:] = -1
    
    headstart = thickness_lay1 + thickness_lay2 + thickness_lay3 - min(wave) + 1
    strt = headstart * np.ones((nlay, nrow, ncol), dtype=np.float32)
    
    flopy.modflow.ModflowBas(mf, 
                              ibound=ibound, 
                              strt=strt,
                              )
    
    # Assign hydraulic conductivities
    hk = np.ones((nlay,nrow,ncol))
    hk[0,:,:] = hk_unconfined
    hk[1,:,:] = hk_confining
    hk[2,:,:] = hk_confined
    
    vka = hk
    
    # Assign storage values    
    if lay1_unconfined:
        ss_unconfined = sy / thickness_lay1
    else:
        ss_unconfined = ss_confined
    
    ss = np.ones((nlay,nrow,ncol))
    ss[0,:,:] = ss_unconfined
    ss[1,:,:] = ss_confined
    ss[2,:,:] = ss_confined
       
    flopy.modflow.ModflowLpf(mf, 
                              hk=hk, 
                              vka=vka, 
                              ss=ss, 
                              sy=sy, 
                              laytyp=laytyp, 
                              ipakcb=53
                              )
    
    # Load Modflow solver
    flopy.modflow.ModflowPcg(mf)
    
    # Define stress period data    
    stress_period_data = {}
    
    for i in range(len(wave)):
        stageleft = wave[i] + headstart
        bound_sp = []
        for il in range(nlay):
            for ir in range(nrow):
                bound_sp.append([il, ir, 0, stageleft, 10000])
        stress_period_data[i] = bound_sp
    
    # Load stress period data
    flopy.modflow.ModflowGhb(mf, 
                              stress_period_data=stress_period_data,
                              )
    
    # Define output
    stress_period_data = {}
    for kper in range(nper):
        for kstp in range(nstp[kper]):
            stress_period_data[(kper, kstp)] = [
                "save head",
                "print head",
            ]
    
    # Load output file
    flopy.modflow.ModflowOc(mf, 
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
    
    # Create hydraulic head array
    head_matrix = np.ones(shape=(nper,ix_max))
    #head_matrix2 = np.ones(shape=(nper,ix_max))
    
    # Extract head files in a readable form
    for i in range(nper):
        head_matrix[i,:] = headobj.get_data(totim=times[i])[2][0,0:ix_max] - headstart
        #head_matrix2[i,:] = headobj.get_data(totim=times[i])[0][0,0:ix_max] - headstart
    
    # stop timer
    end = time.time()
    
    # total time taken
    print("Simulation finished successfully")
    print(f" Runtime [s]: {end - start}")
    
    return([times,gridx[0:ix_max],head_matrix])
   
def model_barrier(dt,
                wave,
                task_name = 'wave_in_aquifer_barrier',
                task_root = './',
                exe_name = "mf2005",
                Lx = 10000,         # x domain length [m]
                Ly = 10,            # y domain length [m]
                ztop = 10,          # Top level [m]
                zbot = 0.0,         # Bottom level [m]
                nlay = 1,           # Number of layers
                nrow = 1,           # Number of rows
                ncol = 15000,       # Number of columns
                hk = 25,            # Horizontal conductivity
#                vka = 25,          # Vertical conductivity
                ss = 5e-05,         # Specific storage [m-1]
                # hk_barrier = 0.01,  # hydraulic conductivity barrier
                d_barrier = 0.1,    # Thickness sheetpiles [m]
                c_barrier = 10.,    # resistance barrier [d]
                dx_max = 2,         # max 
                # sy=0.33,          # specific yield
                # laytyp = -1       # Unconfined (laytyp > 0)
                **settings):

    modelname = '{}{}'.format(task_root,task_name)
    # Start timer
    start = time.time()

    # Calculate x series
    xi= np.logspace(0,np.log10(Lx),ncol+1)-1
    xi = np.exp(np.linspace(start=0, stop=np.log(Lx), num=ncol+1)) - 1
    xi[ncol] = xi[ncol] + 1
    delr = np.diff(xi)                        # Voxel length x direction
    delc = Ly / nrow                          # Voxel length y direction

    ix_max = np.argmax(delr>dx_max)

    # Define the Stress Periods
    nper   = len(wave)                                  # Number of stress periods
       
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
        perlen=dt*np.ones(nper) ,
        nstp=np.ones((nper), dtype=np.int32),
        steady=np.zeros((nper), dtype=bool),
        )
    
    gridx = dis.get_node_coordinates()[1]
    headstart = ztop - min(wave-np.mean(wave)) + 1
    
    # Assign boundary conditions
    flopy.modflow.ModflowBas(
        mf, 
        ibound=np.ones((nlay, nrow, ncol), dtype=np.int32),
        strt= headstart * np.ones((nlay, nrow, ncol), dtype=np.float32),
        )

    ### Assign horizontal hydraulic conductivities
    hks = hk * np.ones((nlay,nrow,ncol))                #  conductivity
    hk_barrier = d_barrier/c_barrier                    # conductivity of barrier [m/d]
    i_barrier  = (np.abs(xi - d_barrier)).argmin()      # Index of x-location where barrier ends
    hks[:,:,0:i_barrier] = hk_barrier      
    # print('index:', hks.shape[2],i_barrier)
    # print('K_values:', hk,hk_barrier)
    # print(hks[:,:,0:i_barrier+20])

    # Assign hydraulic parameters
    flopy.modflow.ModflowLpf(
        mf, 
        hk=hks, 
        vka=hks, 
        ss=ss, 
        # sy=sy, 
        # laytyp=laytyp, 
        ipakcb=53,
    )
    
    # Load Modflow solver
    flopy.modflow.ModflowPcg(mf)
    
    # Define stress period data    
    stress_period_data = {}
    
    for i in range(nper):
        stageleft = wave[i] + headstart
        condleft = 10000	
        stress_period_data[i] = [0, 0, 0, stageleft, condleft]
   
    # Load stress period data
    flopy.modflow.ModflowGhb(
        mf, 
        stress_period_data=stress_period_data
        )
    
    # Define output
    stress_period_out = {}
    for kper in range(nper):
        stress_period_out[(kper, 0)] = [
            "save head",
            "print head",
        ]

    # Load output file
    flopy.modflow.ModflowOc(
        mf, 
        stress_period_data=stress_period_out, 
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

    # Create hydraulic head array
    head_matrix = np.ones(shape=(nper,ix_max))   

    # Extract head files in a readable form
    for i in range(nper):
        head_matrix[i,:] = headobj.get_data(totim=times[i])[0][0,0:ix_max] - headstart

    # stop timer
    end = time.time()

    # total time taken
    print("Simulation finished successfully")
    print(f" Runtime [s]: {end - start}")

    return([times,gridx[0:ix_max],head_matrix])
