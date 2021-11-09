import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

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

##############################################################################
### Plot and File settings
##############################################################################

plt.close('all')
textsize = 8
text_id = ['a','b','c']
dir_cL = '{}/{}_cL/'.format(Ex.task_root,Ex.task_name)

x_range = [50,200,400]

file_fit_cL_leakage = '{}fit_cL_x{:.0f}_leakage.txt'#.format(dir_cL,x_piez)
file_fit_cL_confined = '{}fit_cL_x{:.0f}_confined.txt'#.format(dir_cL,x_piez)

# # ###############################################################################
# # ### Plot Results
# # ###############################################################################

plt.figure(figsize=[7.5,2.5])

for ix,x_piez in enumerate(x_range):
    ax = plt.subplot(1,len(x_range),ix+1)
    fit_results_leakage = np.loadtxt(file_fit_cL_leakage.format(dir_cL,x_piez),delimiter = ',')
    fit_results_confined = np.loadtxt(file_fit_cL_confined.format(dir_cL,x_piez),delimiter = ',')
    cL_range = fit_results_leakage[0,:]

    ax.plot(cL_range,fit_results_leakage[2,:],'-o',c = 'C0',label = 'fit leakage')
    ax.plot(cL_range,fit_results_confined[2,:],'-o',c = 'C1',label = 'fit confined')

    # ax.set_title("Diffusivity fit at $x = {:.0f}$m as function of resistance value".format(x_piez),fontsize=textsize)
    ax.set_xlabel("Resistance $c$ [d]",fontsize=textsize)
    ax.set_xscale('log')
    ax.set_ylim([0,10])
    ax.set_xlim([0.9*cL_range[0],1.1*cL_range[-1]])
    ax.grid(True)

    if ix == 0:
        ax.set_ylabel(r'$\varepsilon_{rel}$ [%] of diffusivity',fontsize=textsize)
    elif ix == 1:
        ax.legend(loc = 'upper left',fontsize=textsize)
    ax.tick_params(axis="both",which="major",labelsize=textsize)
    ax.text(-0.13,-0.15,text_id[ix], bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)
    ax.text(0.7,0.9,'$x = {:.0f}$m'.format(x_piez), bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)

plt.tight_layout()
# plt.savefig('../results/SI_Fig06_Fit_Resistance.pdf')


