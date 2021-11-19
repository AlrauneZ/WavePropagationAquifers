import numpy as np
import matplotlib.pyplot as plt
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
run_sims = False #True #

### select heads at x_piez for all resistance values
select_heads = True #False #

### perform inverse estimation
inverse_est = False #True #

### range of tested resistances

if flow_setting == 'leakage':
    # cL_range = np.logspace(2,5,10 * 3 + 1,endpoint = True)
    # cL_range = np.logspace(4,5,10 + 1,endpoint = True)
    cL_range = np.logspace(0,3,10 * 3 + 1,endpoint = True)
elif flow_setting == 'barrier':
    cL_range = np.logspace(-3,1,10 * 4 + 1,endpoint = True)
    # cL_range = np.logspace(-2,1,10 * 3 + 1,endpoint = True)
    # cL_range = np.logspace(-2.9,-2.1,5 ,endpoint = True)
# print((np.log10(cL_range)))

### observation location
x_piez = 50

Ex = WavePropagationAquifers(
    flow_setting = flow_setting,
    **BC)
Ex.read_wave(**BC) 
Ex.x_piez =  x_piez   

#cL_range = [0.1,0.01,0.001]
# cL_range = np.logspace(-4,0, 4 + 1,endpoint = True)
# cL_range = np.logspace(-4,0,10 * 4 + 1,endpoint = True)

##############################################################################
### Calculate/Load Results of numerical and analytical model 
##############################################################################

##############################################################################
### Specify setting of simulation result files

dir_cL = '{}/{}_cL/'.format(Ex.task_root,Ex.task_name)
if not os.path.exists(dir_cL):
    os.makedirs(dir_cL)

file_results = '{}{}_cL'.format(dir_cL,Ex.task_name)+'_{:0.1f}.p'   ### simulation results for individual resistance values
file_heads_cL = '{}{}_heads_cL_x{:.0f}.txt'.format(dir_cL,Ex.task_name,x_piez)      ### heads at selected x for various resistance values
# file_fit_cL_leakage = '{}{}_fit_cL_x{:.0f}_leakage.txt'.format(dir_cL,Ex.task_name,x_piez)      ### file containing fitting results assuming leakage solution
file_fit_cL_confined = '{}{}_fit_cL_x{:.0f}_confined.txt'.format(dir_cL,Ex.task_name,x_piez)    ### file containing fitting results assuming confined solution

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
            # file_results=file_results.format(abs(np.log10(c_L))),   ### simulation results for individual resistance values
            # file_results = '{}/{}_cL_{:0.1f}.p'.format(dir_cL,Ex.task_name,abs(np.log10(c_L))),
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

if inverse_est:   
    print('\n### Inverse estimation for resistance values ###\n')

    fit_results_confined = np.zeros((3,len(cL_range)))
    fit_results_confined[0,:] = cL_range

    # if flow_setting == 'leakage':
    #     fit_results_leakage = np.zeros((5,len(cL_range)))
    #     fit_results_leakage[0,:] = cL_range

    ### load numerical simulation data at x_piez from file
    h_piez,x_piez,t_piez,cL_range = read_heads_cL(file_heads_cL) 
    Ex.t_piez =  t_piez   
    
    for ic, c_L in enumerate(cL_range):
        Ex.h_piez = h_piez[:,ic] 
        Ex.fit_data_to_analytical_confined()
        fit_results_confined[1,ic] = Ex.diff_fit
        fit_results_confined[2,ic] = Ex.eps_diff

        # if flow_setting == 'leakage':
        #     cS_init = Ex.c_L * Ex.ss*Ex.d_conf
        #     Ex.fit_data_to_analytical_leakage(
        #         p0 = [1e-6,0.5*cS_init],
        #         bounds = ([1e-8,1e-1], [1e-4,1e7]))
        #     print('c_L ', Ex.c_L)
        #     print('cS_init ', cS_init)
        #     print('cS_fit ' , Ex.cS_fit)
        #     # print(relative_difference(Ex.cT_fit,cT_init))
            
        #     fit_results_leakage[1,ic] = Ex.diff_fit
        #     fit_results_leakage[2,ic] = Ex.eps_diff
        #     fit_results_leakage[3,ic] = Ex.cS_fit
        #     fit_results_leakage[4,ic] = Ex.eps_cS
                     
    np.savetxt(file_fit_cL_confined,fit_results_confined,delimiter =',')
    # if flow_setting == 'leakage':
    #     np.savetxt(file_fit_cL_leakage,fit_results_leakage,delimiter =',')

    print('\nInverse estimation results saved to file: \n', file_fit_cL_confined)
    print('\n### Inverse estimation finished ###\n')

    ###############################################################################
    ### Plot Fitting Results
    ###############################################################################

    # textsize = 8    
    # plt.close('all')
    # plt.figure(figsize=[3.75,2.5])
    # plt.plot(cL_range,fit_results_confined[2,:],'-o',c = 'C1',label = 'fit confined solution')
    # if flow_setting == 'leakage':
    #     plt.plot(cL_range,fit_results_leakage[2,:],'-o',c = 'C0',label = 'fit leakage solution')
    
    # plt.title("Diffusivity fit at $x = {:.0f}$m as function of resistance value".format(x_piez),fontsize = textsize)
    # plt.xlabel("Resistance $c$ [d]",fontsize = textsize)
    # plt.ylabel("Relative difference [%]",fontsize = textsize)
    # plt.xscale('log')
    # plt.ylim([0,10])
    # plt.grid(True)
    # plt.legend(fontsize = textsize)
    # plt.savefig('{}/{}_Fit4Resistance_x{:.0f}.png'.format(dir_cL,Ex.task_name,x_piez),dpi=300)

# inverse_est_leakage = False
# file_fit_cL_leakage = '{}{}_fit_cL_x{:.0f}_leakage.txt'.format(dir_cL,Ex.task_name,x_piez)      ### file containing fitting results assuming leakage solution

# if inverse_est_leakage:   
#     print('\n### Inverse estimation for resistance values ###\n')

#     if flow_setting != 'leakage':
#         raise ValueError('Flow setting needs to be leakage')

#     fit_results_leakage = np.zeros((5,len(cL_range)))
#     fit_results_leakage[0,:] = cL_range

#     ### load numerical simulation data at x_piez from file
#     h_piez,x_piez,t_piez,cL_range = read_heads_cL(file_heads_cL) 
#     Ex.t_piez =  t_piez   
    
#     for ic, c_L in enumerate(cL_range):
#         Ex.h_piez = h_piez[:,ic] 
#         Ex.fit_data_to_analytical_leakage()

#         cS_init = Ex.c_L * Ex.ss* Ex.d_conf

#         Ex.fit_data_to_analytical_leakage(
#             p0 = [1e-6,0.5*cS_init],
#             bounds = ([1e-8,1e-1], [1e-4,1e7]))

#         print('c_L ', Ex.c_L)
#         print('cS_init ', cS_init)
#         print('cS_fit ' , Ex.cS_fit)
#         # print(relative_difference(Ex.cT_fit,cT_init))
        
#         fit_results_leakage[1,ic] = Ex.diff_fit
#         fit_results_leakage[2,ic] = Ex.eps_diff
#         fit_results_leakage[3,ic] = Ex.cS_fit
#         fit_results_leakage[4,ic] = Ex.eps_cS
                     
#     np.savetxt(file_fit_cL_leakage,fit_results_leakage,delimiter =',')
#     print('\nInverse estimation results saved to file: \n', file_fit_cL_leakage)
#     print('\n### Inverse estimation finished ###\n')
