import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers #, read_heads_cL

###############################################################################
### Model settings
###############################################################################
BC1 = dict(
    BC_setting =  'tide_wave',   
    normalize_time = True,
    normalize_hs = True,
    )

Ex = WavePropagationAquifers(
    flow_setting = 'confined',
    **BC1)
Ex.read_wave(**BC1) 


file_fit_cL_leakage = '{}/{}_leakage_fit_cL_x{:.0f}_confined.txt'#.format(Ex.task_root,Ex.BC_setting,x_piez)
file_fit_cL_barrier = '{}/{}_barrier_fit_cL_x{:.0f}_confined.txt'#.format(Ex.task_root,Ex.BC_setting,x_piez)
x_piez_range = [50,100,200,400]

###############################################################################
### Plot Results
###############################################################################

plt.close('all')
textsize = 8
text_id = ['a','b','c']

plt.figure(figsize=[7.5,2.5])
ax = plt.subplot(1,2,1)
ax2 = plt.subplot(1,2,2)

for ix,x_piez in enumerate(x_piez_range):
    
    fit_results_leakage = np.loadtxt(file_fit_cL_leakage.format(Ex.task_root,Ex.BC_setting,x_piez),delimiter = ',')
    fit_results_barrier = np.loadtxt(file_fit_cL_barrier.format(Ex.task_root,Ex.BC_setting,x_piez),delimiter = ',')
    
    rel_leakage = 1./fit_results_leakage[0,:]/Ex.hk
    rel_barrier = 0.1/fit_results_barrier[0,:]/Ex.hk
    ax.plot(rel_leakage,fit_results_leakage[2,:],'-o',c = 'C{}'.format(ix),label = 'x = {:.0f}m'.format(x_piez))
    ax2.plot(rel_barrier,fit_results_barrier[2,:],'-o',c = 'C{}'.format(ix),label = 'x = {:.0f}m'.format(x_piez))

ax.set_xlabel("$K_{L}/K$",fontsize=textsize) 
ax.set_xlim([1e-6,2e-4])
ax2.set_xlabel("$K_{B}/K$",fontsize=textsize) 
ax2.set_xlim([0.003,1])

ax.set_title('Leaking aquifer with tide wave BC',fontsize=textsize)
ax.set_xscale('log')
ax.set_ylim([0,10])
ax.grid(True)
ax.set_ylabel(r'$\varepsilon_{rel}$ [%] of diffusivity',fontsize=textsize)
ax.legend(loc = 'upper left',fontsize=textsize)
ax.tick_params(axis="both",which="major",labelsize=textsize)
ax.text(-0.13,-0.15,text_id[0], bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)

ax2.set_title('Flow barrier with tide wave BC',fontsize=textsize)
ax2.set_xscale('log')
ax2.set_ylim([0,10])
ax2.grid(True)
ax2.legend(loc = 'upper right',fontsize=textsize)
ax2.tick_params(axis="both",which="major",labelsize=textsize)
ax2.text(-0.13,-0.15,text_id[1], bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax2.transAxes)

plt.tight_layout()
plt.savefig('../results/Fig06_Fit_Resistance.pdf')

    
