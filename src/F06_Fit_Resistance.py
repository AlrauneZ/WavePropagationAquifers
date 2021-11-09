import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers, read_heads_cL

###############################################################################
### Model settings
###############################################################################
flow_setting = 'leakage'
# flow_setting = 'barrier'

BC1 = dict(
    BC_setting =  'tide_wave',   
    normalize_time = True,
    normalize_hs = True,
    )

Ex = WavePropagationAquifers(
    flow_setting = flow_setting,
    **BC1)
Ex.read_wave(**BC1) 

##############################################################################
### Plot and File settings
##############################################################################

plt.close('all')
textsize = 8
text_id = ['a','b','c']

x_piez = 200
file_fit_cL_leakage = '{}/fit_cL_x{:.0f}_leakage.txt'.format(Ex.task_root,x_piez)
file_fit_cL_confined = '{}/fit_cL_x{:.0f}_confined.txt'.format(Ex.task_root,x_piez)
fit_results_leakage = np.loadtxt(file_fit_cL_leakage,delimiter = ',')
fit_results_confined = np.loadtxt(file_fit_cL_confined,delimiter = ',')

# file_heads_cL = '{}{}_heads_cL_x{:.0f}.txt'.format(dir_cL,Ex.task_name,x_piez)      ### heads at selected x for various resistance values
# file_fit_cL_leakage = '{}{}_fit_cL_x{:.0f}_leakage.txt'.format(dir_cL,Ex.task_name,x_piez)      ### file containing fitting results assuming leakage solution
# file_fit_cL_confined = '{}{}_fit_cL_x{:.0f}_confined.txt'.format(dir_cL,Ex.task_name,x_piez)    ### file containing fitting results assuming confined solution

# # ###############################################################################
# # ### Plot Results
# # ###############################################################################

plt.figure(figsize=[3.75,2.5])

ax = plt.subplot(1,1,1)
# fit_results_leakage = np.loadtxt(file_fit_cL_leakage.format(Ex.task_root,x_piez),delimiter = ',')
# fit_results_confined = np.loadtxt(file_fit_cL_confined.format(Ex.task_root,x_piez),delimiter = ',')
cL_range = fit_results_leakage[0,:]

ax.plot(cL_range,fit_results_leakage[2,:],'-o',c = 'C0',label = 'fit leakage')
ax.plot(cL_range,fit_results_confined[2,:],'-o',c = 'C1',label = 'fit confined')

# ax.set_title("Diffusivity fit at $x = {:.0f}$m as function of resistance value".format(x_piez),fontsize=textsize)
ax.set_xlabel("Resistance $c$ [d]",fontsize=textsize)
ax.set_xscale('log')
ax.set_ylim([0,10])
ax.set_xlim([0.9*cL_range[0],1.1*cL_range[-1]])
ax.grid(True)
ax.set_ylabel(r'$\varepsilon_{rel}$ [%] of diffusivity',fontsize=textsize)
ax.legend(loc = 'upper left',fontsize=textsize)
ax.tick_params(axis="both",which="major",labelsize=textsize)
ax.text(0.8,0.9,'$x = {:.0f}$m'.format(x_piez), bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)

plt.tight_layout()
plt.savefig('../results/Fig06_Fit_Resistance.pdf')


