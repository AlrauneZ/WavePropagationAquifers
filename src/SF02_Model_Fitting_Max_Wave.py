# import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

###############################################################################
### Model settings
###############################################################################
flow_setting = 'confined'
# flow_setting = 'leakage'

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

if flow_setting == 'confined':
    x_piezo = 400
elif flow_setting == 'leakage':
    x_piezo = 200

###############################################################################
### Plot Results
###############################################################################

plt.close('all')
textsize = 8
lw = 2
text_id = ['a','b','c','d','e','f']
text_legend = ['Tidal wave', 'Square wave', 'River levels']

fig = plt.figure(figsize=[7.5,4.5])

for ii, BC in enumerate([BC1,BC2,BC3]):
    ax = fig.add_subplot(2,3,ii+1)
    ax2 = fig.add_subplot(2,3,ii+4)

    ##############################################################################
    ### Calculate/Load Results of numerical and analytical model 
    ##############################################################################
    Ex = WavePropagationAquifers(**BC, flow_setting = flow_setting)
    Ex.read_wave(**BC)
    Ex.prepare_piezometric_data(x = x_piezo, write_to_file=False)
    Ex.fit_data_to_analytical()

    ax.plot(Ex.t_piez,Ex.h_piez,lw = lw,c = 'C0',label = 'observed')
    ax.plot(Ex.t_piez,Ex.head_ana_fit,ls = '--',lw = lw-.5,c = 'C1',label = 'model fit')

    ax2.plot(Ex.t_piez,Ex.wave_piez_max,lw = lw,c = 'C2',label = 'observed')
    ax2.plot(Ex.t_piez,Ex.wave_ana_max,ls = '--',lw = lw-.5,c = 'C3',label = 'model fit')

    if ii == 0:
        ax.set_ylabel(r"Water table $h(t)$ [m]",fontsize=textsize)
        ax2.set_ylabel(r"$w_{max}(t)$ [m]",fontsize=textsize)
    elif ii == 2:
        ax.legend(loc='center right',fontsize=textsize,ncol=1)
        ax2.set_ylim([ax2.get_ylim()[0],ax2.get_ylim()[1]*1.3])
        ax2.legend(loc='upper center',fontsize=textsize,ncol=2)
    ax.text(0.7,0.9,r'$\varepsilon = {:.1f}\%$'.format(Ex.eps_diff), bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)

    ax.grid(True)
    ax.tick_params(axis="both",which="major",labelsize=textsize)
    ax.set_title(text_legend[ii],fontsize=textsize)
    ax.text(-0.13,-0.15,text_id[ii], bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)

    ax2.set_xlabel(r"Time $t$ [days]",fontsize=textsize)
    ax2.grid(True)
    ax2.tick_params(axis="both",which="major",labelsize=textsize)
    ax2.text(-0.13,-0.15,text_id[ii+3], bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax2.transAxes)

plt.tight_layout()
# plt.savefig('../results/SI_Fig02_Model_Fitting_Max_Wave_{}.png',dpi=300.format(flow_setting))   
plt.savefig('../results/SI_Fig02_Model_Fitting_Max_Wave_{}.pdf'.format(flow_setting)) 


