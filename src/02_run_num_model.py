import numpy as np
import matplotlib.pyplot as plt

from WavePropagationAquifers import WavePropagationAquifers

###############################################################################
### Select Model settings
###############################################################################

BC1 = dict(
    task_name =  'tide_wave',    
    cut_input = 1310,
    # normalize_time = True,
    normalize_time = False,
    normalize_hs = True,
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
    coarsen_time= 6,
    )

Ex = WavePropagationAquifers(**BC1)
Ex.read_wave(**BC1) 

# Ex = WavePropagationAquifers(**BC2)
# Ex.read_wave(**BC2)

# Ex = WavePropagationAquifers(**BC3)
# Ex.read_wave(**BC3)

##############################################################################
### Calculate/Load Results of numerical and analytical model 
##############################################################################

Ex.run_numerical_model()#write_to_file = False)
# Ex.read_num_from_pickle() 
Ex.select_times()
# indices = Ex.t_indices
# t_select = Ex.t_select
# t_rel = Ex.t_rel

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

i0 = Ex.t_indices[14]
i1 = Ex.t_indices[10]
i2 = Ex.t_indices[3]

# i0 = -1
# i1 = int(0.21*len(Ex.times_num))
# i2 = int(0.71*len(Ex.times_num))

plt.plot(Ex.x_num, Ex.head_tx_num[i0,:], '-', color='C0', label="T = {:.2f} (final time)".format(Ex.times_num[i0]))
plt.plot(Ex.x_num, Ex.head_tx_num[i1,:], '-', color='C1', label="T = {:.2f} (0.71 sim time)".format(Ex.times_num[i1]))
plt.plot(Ex.x_num, Ex.head_tx_num[i2,:], '-', color='C2', label="T = {:.2f} (0.21 sim time)".format(Ex.times_num[i2]))

plt.title("Numerical result (spatial distribution)")
plt.xlabel("Distance from coast $x$ [m]")
plt.ylabel("Hydraulic head h(x) [m]")
plt.grid(True)
plt.legend()
# plt.savefig('{}/Numerical_Results_hx_{}.png'.format(Ex.task_root,Ex.task_name),png = 300)


plt.figure(2)
xloc = 850
ix =  (np.abs(Ex.x_num - xloc)).argmin()

plt.plot(Ex.times_num, Ex.head_tx_num[:,ix], '-', color='C0',label='h(t) at x = {}m'.format(xloc))
plt.plot(Ex.wave_time,Ex.wave, color='C1',ls = '--',label = 'h0 (input wave)')
plt.title("Numerical result (time series at x = {}m)".format(xloc))
plt.xlabel("Time $t$ [?]")
plt.ylabel("Hydraulic head h(t) [m]")
plt.grid(True)
plt.legend()
# plt.savefig('{}/Numerical_Results_ht_{}.png'.format(Ex.task_root,Ex.task_name),png = 300)

