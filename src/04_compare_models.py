import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

##############################################################################
### Specify Model settings
##############################################################################

BC1 = dict(
    task_name =  'tide_wave',    
    cut_input = 1310,
    normalize_time = False, ### chose this option to allow comparison with Wout's results
    # normalize_time = True, ### proper option to use later
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

# Ex = WavePropagationAquifers(**BC3)
# Ex.read_wave(**BC3) 

##############################################################################
### Calculate/Load Results of numerical and analytical model 
##############################################################################
Ex.select_times()
t_rel = Ex.t_rel

# Ex.run_numerical_model()
Ex.read_num_from_pickle()
x_num = Ex.x_num
t_num = Ex.times_num[Ex.t_indices] 
head_num = Ex.head_tx_num[Ex.t_indices,:]

Ex.run_analytical_model()
# Ex.read_ana_from_pickle()
x_ana = Ex.x_ana
t_ana = Ex.times_ana
head_ana = Ex.head_tx_ana

###############################################################################
### Plot Results
###############################################################################

plt.close('all')

ti_range = [3,10,14] ### time ranges for tide_wave and river
# ti_range = [1,8,14]   ### time range Wout for square wave

for i,ti in enumerate(ti_range):
    plt.plot(x_num, head_num[ti,:], '-', color='C{}'.format(i), label="T={:.2f} Numerical".format(t_num[ti]))
    plt.plot(Ex.x_ana, Ex.head_tx_ana[ti,:], 'o', color='C{}'.format(i), label="T={:.2f} Analytical".format(t_num[ti]))
    # plt.plot(x_num, head_num[ti,:], '-', color='C{}'.format(i), label="T={:.2f} Numerical".format(t_rel[ti]))
    # plt.plot(Ex.x_ana, Ex.head_tx_ana[ti,:], 'o', color='C{}'.format(i), label="T={:.2f} Analytical".format(t_rel[ti]))

# plt.ylim(-0.2,1.2)
plt.xlim([-20,Ex.x_ana[-1]+20])
plt.title("Comparisson of Numerical vs Analytical result")
plt.xlabel("Distance from coast [m]")
plt.ylabel("Hydraulic head [m]")
plt.grid(True)
plt.legend()
# plt.savefig('{}/Compare_Results_hx_{}.png'.format(Ex.task_root,Ex.task_name),png = 300)
