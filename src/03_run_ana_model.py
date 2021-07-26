import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

##############################################################################
### Specify Model settings
##############################################################################

BC1 = dict(
    task_name =  'tide_wave',    
    cut_input = 1310,
    normalize_time = False,
    # normalize_time = True,
    normalize_hs = True
    )
BC2 = dict(
    task_name = 'square_wave',    
    cut_input = 753,
    normalize_time = False,
    normalize_hs = False
    )
BC3 = dict(
    task_name = 'river', 
    normalize_time = True,
    normalize_hs = True,
    )

Ex = WavePropagationAquifers(**BC1)
Ex.read_wave(**BC1) 

# Ex = WavePropagationAquifers(**BC2)
# Ex.read_wave(**BC2) 

Ex = WavePropagationAquifers(**BC3)
Ex.read_wave(**BC3) 

##############################################################################
### Calculate/Load Results of numerical and analytical model 
##############################################################################

Ex.run_analytical_model()#write_to_file = False)
# Ex.read_ana_from_pickle()

###############################################################################
### Plot Results
###############################################################################

plt.close('all')
plt.figure(1)

""" These are the indices according to Wout's relative time specification
BUT they are not correctly representing the relative simulation time!
They are refering to the end point of the analytical time range which is 
shorter then the wave_time length (due to time array specifications)
A corrected version is implemented in the class, but not yet used to allow for
comparison on results with plots in the manuscript
"""

ti_range = [3,10,14]
# ti_range = [1,8,14]   ### time range Wout for square wave

plt.plot(Ex.x_ana, Ex.head_tx_ana[ti_range[0],:], 'o:', color='blue', label="T=0.21 Semi-analytic")
plt.plot(Ex.x_ana, Ex.head_tx_ana[ti_range[1],:], 'o:', color='red', label="T=0.71 Semi-analytic")
plt.plot(Ex.x_ana, Ex.head_tx_ana[ti_range[2],:], 'o:', color='orange', label="T=1.00 Semi-analytic")

#plt.ylim(-0.2,1.2)

plt.title("Analytical result (spatial distribution)")
plt.xlabel("Distance from coast $x$ [m]")
plt.ylabel("Hydraulic head h(x) [m]")
plt.grid(True)
plt.legend()
# plt.savefig('{}/Analytical_Results_hx_{}.png'.format(Ex.task_root,Ex.task_name),png = 300)

plt.figure(2)
xloc = 850
ix =  (np.abs(Ex.x_ana - xloc)).argmin()

plt.plot(Ex.times_ana, Ex.head_tx_ana[:,ix], '-', color='C0',label='h(t) at x = {}m'.format(xloc))
plt.plot(Ex.wave_time,Ex.wave, color='C1',ls = '--',label = 'h0 (input wave)')
plt.title("Analytical result (time series at x = {}m)".format(xloc))
plt.xlabel("Time $t$ [?]")
plt.ylabel("Hydraulic head h(t) [m]")
plt.grid(True)
plt.legend()
# plt.savefig('{}/Analytical_Results_ht_{}.png'.format(Ex.task_root,Ex.task_name),png = 300)
