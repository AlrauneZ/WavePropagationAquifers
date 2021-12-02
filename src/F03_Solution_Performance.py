import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

###############################################################################
### Model settings
###############################################################################

flow_setting = 'confined'   
flow_setting = 'leakage'   

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

###############################################################################
### Plot Results
###############################################################################

range_BCs = [BC1,BC2,BC3]

plt.close('all')
textsize = 8
text_legend = ['Tidal wave', 'Square wave', 'River levels']
text_id = ['a','b','c']

if flow_setting == 'confined':
    times_rel_compare = [0.3,0.6,0.9] # relative simulation time points for comparison
    x_max = [2000,2000,6000]
    fig_name  = '../results/Fig03_Solution_Performance_{}.pdf'.format(flow_setting)
elif flow_setting == 'leakage':
    x_max = 800.*np.ones(len(range_BCs))
    times_rel_compare = [0.3,0.5,0.7,0.9] # relative simulation time points for comparison
    fig_name  = '../results/SI_Fig03_Solution_Performance_{}.pdf'.format(flow_setting)

fig = plt.figure(figsize=[7.5,2.5])
for ii, BC in enumerate(range_BCs):
    ax = fig.add_subplot(1,3,ii+1)

    ### Calculate/Load Results of numerical and analytical model 
    Ex = WavePropagationAquifers(
        flow_setting = flow_setting,
        **BC)

    Ex.read_wave(**BC)
    Ex.read_num_from_pickle() 
    Ex.select_sim_times(times_select = times_rel_compare)
    Ex.run_analytical_model(     # calculate solution at specified time points
                            x_ana = np.linspace(0,x_max[ii],20),
                            t_ana= times_rel_compare,
                            t_rel = True)

    for it,t_rel in enumerate(times_rel_compare):
        ax.plot(Ex.x_ana, Ex.head_tx_ana[it,:], 'o',ms = 4, color='C{}'.format(it), label="T={:.1f} ana".format(t_rel))
        ax.plot(Ex.x_num, Ex.head_num_select[it,:], '-', color='C{}'.format(it), label="T={:.1f} num".format(t_rel))

    ax.set_xlabel(r"Distance from interface $x$ [m]",fontsize=textsize)
    if ii ==0:
        ax.set_ylabel(r"Water table $h(x)$ [m]",fontsize=textsize)
    # elif ii ==1:
        ax.legend(loc='upper right',fontsize=textsize-1,ncol=1)
    ax.set_xlim([-20,x_max[ii]+20])

    ax.grid(True)
    ax.tick_params(axis="both",which="major",labelsize=textsize)
    ax.set_title(text_legend[ii],fontsize=textsize)
    ax.text(-0.13,-0.15,text_id[ii], bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)

plt.tight_layout()
# plt.savefig(fig_name,dpi=300)   
plt.savefig(fig_name)   


