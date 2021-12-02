import numpy as np
from WavePropagationAquifers import WavePropagationAquifers,damping_coefficient

###############################################################################
### Select Model settings
###############################################################################

BC1 = dict(
    BC_setting = 'tide_wave', 
    # normalize_time = True,
    normalize_time = False,
    normalize_hs = True,
    )
BC2 = dict(
    BC_setting = 'square_wave', 
    normalize_time = False,
    normalize_hs = True
    )
BC3 = dict(
    BC_setting = 'river_levels', 
    normalize_time = True,
    normalize_hs = True,
    )

file_name = '../results/{}_diff_fit_locs.txt'

###############################################################################
### SETTINGS 
flow_setting = 'confined'
BC = BC1
x_range = np.arange(10,6510,10) ### range of piezometric locations 

###############################################################################
### process data on inverse estimation at various locations

### setup instance of model and read in BC input wave
Ex = WavePropagationAquifers(
    flow_setting = flow_setting,
    **BC)
Ex.read_wave(**BC) 

### runs in numerical simualation to produce "actificial measurements"
# Ex.run_numerical_model(**BC) 

### reads in numerical simualation results from file
Ex.read_num_from_pickle()

### lists of estimation results
x_piez_range = []
diff_fit_range = []
eps_diff_range = []
A_piez_range = []

for ix,x in enumerate(x_range):

    ### extracting head observation time series as selected spatial location x
    x_piez, t_piez, h_piez = Ex.prepare_piezometric_data(x=x, write_to_file=False)
    x_piez_range.append(x_piez) # actual x-locations (of numerical data)
    A_piez_range.append(Ex.extract_piez_wave_component()[0])

    ### inverse estimation of diffusivity through fitting to analytical solution
    diff_fit, eps_diff = Ex.fit_data_to_analytical_confined(verbose = False)
    diff_fit_range.append(diff_fit)
    eps_diff_range.append(eps_diff)

### calculating dominant wave component and damping coefficients
A_wave = Ex.extract_dominant_input_wave_component()[0]
dc = damping_coefficient(A_piez_range,A_wave)

### saving all estimates to files
fit_locs = np.vstack((x_piez_range,diff_fit_range,eps_diff_range,dc))
np.savetxt(file_name.format(Ex.task_name),fit_locs.T,delimiter = ',',header = 'x, diffusivity, rel diff, damping coeff')
