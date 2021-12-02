import numpy as np
# import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers, save_heads_cL, read_heads_cL
import os

###############################################################################
### Model settings
###############################################################################

BC1 = dict(
    BC_setting =  'tide_wave',   
    normalize_time = True,
    normalize_hs = True,
    )
BC2 = dict(
    BC_setting = 'square_wave',    
    normalize_time = False,
    normalize_hs = False
    )
BC3 = dict(
    BC_setting = 'river_levels', 
    normalize_time = True,
    normalize_hs = True,
    )

##############################################################################
### Settings
##############################################################################
### input wave BC
BC = BC1 

### flow setting of numerical simulation
flow_setting = 'leakage'
# flow_setting = 'barrier'

### run simulation for individual resistance values
run_sims = True #False #

### select heads at x_piez for all resistance values
select_heads = True #False #

### perform inverse estimation of numerical results with confined solution
inverse_est_confined = True #False #

### perform inverse estimation of numerical results with leakage solution
inverse_est_leakage =  True # False # 
### When inverse estimation with leakage solution: fix resistance

### observation location
x_piez = 400

### range of tested resistances
if flow_setting == 'leakage':
    cL_range = np.logspace(0,3,10 * 3 + 1,endpoint = True)
elif flow_setting == 'barrier':
    cL_range = np.logspace(-3,1,10 * 4 + 1,endpoint = True)

##############################################################################
### Calculate/Load Results of numerical and analytical model 
##############################################################################

Ex = WavePropagationAquifers(
    flow_setting = flow_setting,
    **BC)
Ex.read_wave(**BC) 
Ex.x_piez =  x_piez   

### Specify separate directory for all numerical simulation results
dir_cL = '{}/{}_cL/'.format(Ex.task_root,Ex.task_name)
if not os.path.exists(dir_cL):
    os.makedirs(dir_cL)

### files with simulation results for individual resistance values
file_results = '{}{}_cL'.format(dir_cL,Ex.task_name)+'_{:0.1f}.p'   

### files with heads at selected x for various resistance values
file_heads_cL = '{}{}_heads_cL_x{:.0f}.txt'.format(dir_cL,Ex.task_name,x_piez)      

### files containing inverse estimation results: for different combination of 
### aquifer simulation setting and analytical solution

### file containing fitting results assuming confined solution
file_fit_cL_confined = '{}_fit_cL_x{:.0f}_confined.txt'.format(Ex.task_name,x_piez)    

### file containing fitting results assuming leakage solution
file_fit_cL_leakage = '{}_fit_cL_x{:.0f}_leakage.txt'.format(Ex.task_name,x_piez)      

##############################################################################
### run simulation for individual resistance values
##############################################################################

if run_sims:   
    print('### Running Simulation for various resistance values ###\n')
    ### run all numerical simulation to get observational data
    for ic, c_L in enumerate(cL_range):
        if flow_setting == 'leakage':
            Ex.c_L = c_L    
        elif flow_setting == 'barrier':
            Ex.c_barrier  = c_L

        print('##################')
        Ex.run_numerical_model(
            dir_sim = '{}/Modflow_cL/'.format(dir_cL),
            write_to_file = True,
            file_results=file_results.format(np.log10(c_L)),   ### simulation results for individual resistance values
            coarsen_x=1.,
            )

    print('### Simulations finished ###\n')
    
##############################################################################
### Postprocess numerical simulation data 
### --> store all heads at x_piez for all resistance values
##############################################################################

if select_heads:
    print('### Postprocess numerical simulation data ###\n')

    h_piez_range = []
    for ic, c_L in enumerate(cL_range):
        Ex.c_L = c_L    
    
        ### Postprocess numerical simulation data --> store all heads at x_piez for all resistance values
        Ex.read_num_from_pickle(
            file_results = file_results.format(np.log10(c_L))) 
        Ex.prepare_piezometric_data(x = x_piez, write_to_file=False)
        h_piez_range.append(Ex.h_piez)
   
    ### save postprocessed numerical simulation data at x_piez to one file
    save_heads_cL(
        np.array(h_piez_range).T,
        x_piez,
        Ex.t_piez,
        cL_range,
        file_results_cL = file_heads_cL,
        )
    print('\n### Postprocessing finished ###\n')

##############################################################################
### Inverse parameter estimation for analytical solution(s)
##############################################################################

if inverse_est_confined:   
    print('\n### Inverse estimation for resistance values ###\n')

    fit_results_confined = np.zeros((3,len(cL_range)))
    fit_results_confined[0,:] = cL_range

    ### load numerical simulation data at x_piez from file
    h_piez,x_piez,t_piez,cL_range = read_heads_cL(file_heads_cL) 
    Ex.t_piez =  t_piez   
    
    for ic, c_L in enumerate(cL_range):
        Ex.h_piez = h_piez[:,ic] 
        Ex.fit_data_to_analytical_confined()
        fit_results_confined[1,ic] = Ex.diff_fit
        fit_results_confined[2,ic] = Ex.eps_diff
                     
    np.savetxt(file_fit_cL_confined,fit_results_confined,delimiter =',')
    print('\nInverse estimation results saved to file: \n', file_fit_cL_confined)
    print('\n### Inverse estimation finished ###\n')


if inverse_est_leakage:   
    print('\n### Inverse estimation for resistance values ###\n')

    if flow_setting != 'leakage':
        raise ValueError('Flow setting needs to be set to leakage')

    fit_results_leakage = np.zeros((5,len(cL_range)))
    fit_results_leakage[0,:] = cL_range

    ### load numerical simulation data at x_piez from file
    h_piez,x_piez,t_piez,cL_range = read_heads_cL(file_heads_cL) 
    Ex.t_piez =  t_piez   
    
    for ic, c_L in enumerate(cL_range):
        Ex.h_piez = h_piez[:,ic] 
        Ex.c_L = c_L

        cS_init = Ex.c_L * Ex.ss* Ex.d_conf
        diff_init = Ex.ss/Ex.hk

        Ex.fit_data_to_analytical_leakage(
            p0 = [0.5*diff_init,0.5*cS_init],
            bounds = ([1e-9,1e-6], [1e-3,10]), # lower bounds, upper bounds
            )
        fit_results_leakage[1,ic] = Ex.diff_fit
        fit_results_leakage[2,ic] = Ex.eps_diff
        fit_results_leakage[3,ic] = Ex.cS_fit
        fit_results_leakage[4,ic] = Ex.eps_cS
        np.savetxt(file_fit_cL_leakage,fit_results_leakage,delimiter =',')
        print('\nInverse estimation results saved to file: \n', file_fit_cL_leakage)
                           
    print('\n### Inverse estimation finished ###\n')


