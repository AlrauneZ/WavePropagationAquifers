#import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

###############################################################################
### Model settings
###############################################################################

BC1 = dict(
    BC_setting = 'tide_wave', 
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

plt.close('all')
lw = 1
textsize = 8
text_id = ['a','b','c']
text_legend = ['Tidal wave', 'Square wave', 'River levels']
plt.figure(figsize=[7.5,2.2])

for ii, BC in enumerate([BC1,BC2,BC3]):
    ax=plt.subplot(1,3,ii+1)
    
    Ex = WavePropagationAquifers(**BC)
    Ex.read_wave(**BC)

    ax.plot(Ex.wave_time,Ex.wave,'-',lw=lw,label=text_legend[ii])        
    ax.set_xlabel(r"Time $t$ [days]",fontsize=textsize)
    if ii ==0:
        ax.set_ylabel(r"Water table $h(t)$ [m]",fontsize=textsize)
    ax.grid(True)
    ax.tick_params(axis="both",which="major",labelsize=textsize)
    # ax.legend(loc='upper right',fontsize=textsize)
    ax.set_title(text_legend[ii],fontsize=textsize)
    ax.text(-0.13,-0.15,text_id[ii], bbox=dict(facecolor='w', alpha=1,boxstyle='round'),fontsize=textsize, transform=ax.transAxes)
    
plt.tight_layout()
# plt.savefig('../results/Fig02_Input_BC_Waves.png',dpi=300)   
plt.savefig('../results/Fig02_Input_BC_Waves.pdf')   
