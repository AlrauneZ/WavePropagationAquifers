import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

###############################################################################
### Select Model settings
###############################################################################

BC3 = dict(
    task_name = 'river', 
    normalize_time = True,
    normalize_hs = True,
    coarsen_time= 6,
    )
BC3b = dict(
    task_name = 'river', 
    normalize_time = True,
    normalize_hs = True,
    )

### coarsened river data
Ex = WavePropagationAquifers(**BC3)
Ex.read_wave(**BC3)

### original fine time resolution of river data
Ex2 = WavePropagationAquifers(**BC3b)
Ex2.read_wave(**BC3b)

##############################################################################
### Calculate/Load Results of numerical and analytical model 
##############################################################################

# Ex.run_numerical_model()#write_to_file = False)
Ex.read_num_from_pickle() 
Ex.select_times()

Ex2.read_num_from_pickle(file_results = r'../results/river_numerical_01.p')
Ex2.select_times()

###############################################################################
### Plot Results
###############################################################################

plt.close('all')
plt.figure(1)

### adapt displayed time indeces given different resolutions
i0 = Ex.t_indices[14]+9
i1 = Ex.t_indices[10]+6
i2 = Ex.t_indices[3]+1

i0b = Ex2.t_indices[14]
i1b = Ex2.t_indices[10]
i2b = Ex2.t_indices[3]

plt.plot(Ex.x_num, Ex.head_tx_num[i0,:], '-', color='C0', label="T = {:.4f} (coarse) ".format(Ex.times_num[i0]))
plt.plot(Ex2.x_num[::100], Ex2.head_tx_num[i0b,::100], 'o', color='C0', label="T = {:.4f} (fine)".format(Ex2.times_num[i0b]))

plt.plot(Ex.x_num, Ex.head_tx_num[i1,:], '-', color='C1', label="T = {:.4f} (coarse)".format(Ex.times_num[i1]))
plt.plot(Ex2.x_num[::100], Ex2.head_tx_num[i1b,::100], 'o', color='C1', label="T = {:.4f} (fine)".format(Ex2.times_num[i1b]))

plt.plot(Ex.x_num, Ex.head_tx_num[i2,:], '-', color='C2', label="T = {:.4f} (fine)".format(Ex.times_num[i2]))
plt.plot(Ex2.x_num[::100], Ex2.head_tx_num[i2b,::100], 'o', color='C2', label="T = {:.4f} (coarse)".format(Ex2.times_num[i2b]))

plt.title("Numerical result (spatial distribution)")
plt.xlabel("Distance from coast $x$ [m]")
plt.ylabel("Hydraulic head h(x) [m]")
plt.grid(True)
plt.legend()
plt.savefig('../results/Check_Numerical_Results_Coarsened_River',png = 300)