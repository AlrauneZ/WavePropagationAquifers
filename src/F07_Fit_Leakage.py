import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers,read_heads_cL

###############################################################################
### Model settings
###############################################################################
BC1 = dict(
    BC_setting =  'tide_wave',   
    normalize_time = True,
    normalize_hs = True,
    )

Ex = WavePropagationAquifers(
    flow_setting = 'leakage',
    **BC1)

Ex.read_wave(**BC1) 

dir_cL = '{}/{}_cL/'.format(Ex.task_root,Ex.task_name)
cL_range = np.logspace(0,3,10 * 3 + 1,endpoint = True)
x_piez = 200 

# file_results = '{}{}_cL'.format(dir_cL,Ex.task_name)+'_{:0.1f}.p'   ### simulation results for individual resistance values
file_heads_cL = '{}{}_heads_cL_x{:.0f}.txt'.format(dir_cL,Ex.task_name,x_piez)      ### heads at selected x for various resistance values
file_fit_cL_leakage = '{}{}_fit_cL_x{:.0f}_leakage.txt'.format(dir_cL,Ex.task_name,x_piez)      ### file containing fitting results assuming leakage solution

##############################################################################
### Plot and File settings
##############################################################################

print('\n### Inverse estimation for resistance values ###\n')

fit_results_leakage = np.zeros((5,len(cL_range)))
fit_results_leakage[0,:] = cL_range

### load numerical simulation data at x_piez from file
h_piez,x_piez,t_piez,cL_range = read_heads_cL(file_heads_cL) 
Ex.t_piez =  t_piez   
Ex.x_piez =  x_piez   
    
for ic, c_L in enumerate(cL_range):
    Ex.h_piez = h_piez[:,ic] 
    Ex.c_L = c_L
    cS_init = Ex.c_L * Ex.ss* Ex.d_conf

    diff_init = Ex.ss/Ex.hk
    Ex.fit_data_to_analytical_leakage(
        fix_cS = cS_init,
        p0 = [2*diff_init,cS_init],
        # p0 = [0.5*diff_init,0.5*cS_init],
        bounds = ([1e-9,1e-3], [1e-5,10]),
        # bounds = ([1e-8,1e-1], [1e-4,1e7])
        )

    print('c_L ', Ex.c_L)
    print('cS_init ', cS_init)
    print('cS_fit ' , Ex.cS_fit)
    # print(relative_difference(Ex.cT_fit,cT_init))
    
    fit_results_leakage[1,ic] = Ex.diff_fit
    fit_results_leakage[2,ic] = Ex.eps_diff
    fit_results_leakage[3,ic] = Ex.cS_fit
    fit_results_leakage[4,ic] = Ex.eps_cS
                 
# np.savetxt(file_fit_cL_leakage,fit_results_leakage,delimiter =',')
# print('\nInverse estimation results saved to file: \n', file_fit_cL_leakage)
# print('\n### Inverse estimation finished ###\n')

d_barrier = 1.
hk = 25

# # # ###############################################################################
# # # ### Plot Results
# # # ###############################################################################

plt.close('all')
textsize = 8
# text_id = ['a','b','c']

plt.figure(figsize=[3.75,2.5])

ax = plt.subplot(1,1,1)

# fit_results_leakage = np.loadtxt(file_fit_cL_leakage.format(dir_cL,x_piez),delimiter = ',')
# fit_results_confined = np.loadtxt(file_fit_cL_confined.format(dir_cL,x_piez),delimiter = ',')
# cL_range = fit_results_leakage[0,:]

ax.plot(cL_range,fit_results_leakage[2,:],'-o',c = 'C0',label = 'eps diff')
ax.plot(cL_range,fit_results_leakage[4,:],'-o',c = 'C1',label = 'eps cS')

# x_rel = d_barrier/cL_range/hk
# ax.plot(x_rel,fit_results_leakage[2,:],'-o',c = 'C0',label = 'fit leakage')
# ax.set_xlabel("Relative resistance $K_{barrier}/K_{aquifer}$",fontsize=textsize)

# ax.set_title("Diffusivity fit at $x = {:.0f}$m as function of resistance value".format(x_piez),fontsize=textsize)
ax.set_xscale('log')
ax.set_ylim([0,10])
# ax.set_xlim([0.9*cL_range[0],1.1*cL_range[-1]])
# ax.grid(True)

ax.set_xlabel("Resistance $c$ [d]",fontsize=textsize)
ax.set_ylabel(r'$\varepsilon_{rel}$ [%] of diffusivity',fontsize=textsize)
ax.legend(loc = 'upper left',fontsize=textsize)
ax.tick_params(axis="both",which="major",labelsize=textsize)

# ax.text(-0.13,-0.15,text_id[ix], bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)
ax.text(0.7,0.9,'$x = {:.0f}$m'.format(x_piez), bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)

# plt.tight_layout()
# # plt.savefig('../results/Fig07_Fit_Leakage.pdf')


