import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

###############################################################################
### Model settings
###############################################################################
flow_setting = 'confined'   

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

file_name = '../results/{}_diff_fit_locs.txt'

###############################################################################
### Plot Results
###############################################################################

plt.close('all')
textsize = 8
lw = 2
text_legend = ['Tidal wave', 'Square wave', 'River levels']

fig = plt.figure(figsize=[7.5,2.5])
ax = plt.subplot(121)
ax2 = plt.subplot(122)

for ii, BC in enumerate([BC1,BC2,BC3]):

    Ex = WavePropagationAquifers(
        flow_setting = flow_setting,   
        **BC)
    
    Ex.read_wave(**BC)   
    fit_locs = np.loadtxt(file_name.format(Ex.task_name),delimiter = ',',skiprows=1)

    ax.plot(fit_locs[:,0],fit_locs[:,2],'-',lw = 2,label = text_legend[ii])
    ax2.plot(fit_locs[:,3],fit_locs[:,2],'-',lw = 2)#,label = text_legend[ii])

ax.set_xlabel('Distance to Interface $x$ [m]',fontsize=textsize)
ax.set_ylabel(r'Relative Difference $\varepsilon(x)$ [%]',fontsize=textsize)
ax.grid(True)
ax.set_ylim([0,10])
ax.set_xlim([0,3000])
ax.legend(loc='best',fontsize=textsize,ncol=1)
ax.tick_params(axis="both",which="major",labelsize=textsize)

ax2.set_xlabel('Damping Coefficient [-]',fontsize=textsize)
ax2.set_ylim([0,10])
ax2.set_xlim([0,1.0])
ax2.grid(True)
ax2.tick_params(axis="both",which="major",labelsize=textsize)

plt.tight_layout()
# plt.savefig('../results/Fig05_Impact_Location_confined.png',dpi=300)   
plt.savefig('../results/Fig05_Impact_Location_confined.pdf')   


