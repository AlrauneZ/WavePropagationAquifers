import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers #,read_heads_cL

###############################################################################
### Model settings
###############################################################################
BC1 = dict(
    BC_setting =  'tide_wave',   
    normalize_time = True,
    normalize_hs = True,
    )

### file containing fitting results assuming leakage solution
file_fit_cL_leakage = '../results/{}_fit_cL_x{:.0f}_leakage.txt'      
x_piez_range = [50,200]

Ex = WavePropagationAquifers(
    flow_setting = 'leakage',
    **BC1)
Ex.read_wave(**BC1) 

###############################################################################
### Plot Results
###############################################################################

plt.close('all')
textsize = 8

plt.figure(figsize=[3.75,2.5])
ax = plt.subplot(1,1,1)

for ix,x_piez in enumerate(x_piez_range):
    fit_results_leakage = np.loadtxt(file_fit_cL_leakage.format(Ex.task_name,x_piez), delimiter = ',') 
    cL_range = fit_results_leakage[0,:] 
    eps_diff = fit_results_leakage[2,:]
    eps_cS = fit_results_leakage[4,:]
       
    ax.plot(cL_range,eps_diff,'-o',ms = 4,c = 'C{}'.format(2*ix),label = r'$\varepsilon$ (diff) at x = {:.0f} m'.format(x_piez))
    ax.plot(cL_range,eps_cS,'-s',ms = 4,c = 'C{}'.format(2*ix+1),label = r'$\varepsilon$ (cS) at x = {:.0f} m'.format(x_piez))
 
ax.set_xscale('log')
ax.set_ylim([0,10])
ax.grid(True)
ax.set_xlim([cL_range[0],cL_range[-1]])

ax.set_xlabel("Resistance of confining layer $c_L$  [d]",fontsize=textsize)
ax.set_ylabel(r'Relative differences $\varepsilon_{rel}$ [%]',fontsize=textsize)
ax.legend(loc = 'upper right',fontsize=textsize)
ax.tick_params(axis="both",which="major",labelsize=textsize)

plt.tight_layout()
plt.savefig('../results/Fig07_Fit_Leakage.pdf')


