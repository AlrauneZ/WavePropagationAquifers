import numpy as np
from WavePropagationAquifers import WavePropagationAquifers

###############################################################################
### Model settings
###############################################################################

# flow_setting = 'leakage'   
# flow_setting = 'confined'   
flow_setting = 'barrier'   

BC1 = dict(
    BC_setting =  'tide_wave',   
    normalize_time = True,
    normalize_hs = True,
    coarsen_x=1.,
    )
BC2 = dict(
    BC_setting = 'square_wave',    
    normalize_time = False,
    normalize_hs = False,
    coarsen_x=1.,
    )
BC3 = dict(
    BC_setting = 'river_levels', 
    normalize_time = True,
    normalize_hs = True,
    Lx = 50000,
    dx_max = 5,
    coarsen_x=1.,
    )

BC = BC1

Ex = WavePropagationAquifers(
    flow_setting = flow_setting,
    **BC)

Ex.read_wave(**BC) 

##############################################################################
### Calculate/Load Results of numerical and analytical model 
##############################################################################

Ex.run_numerical_model(**BC) 
#Ex.read_num_from_pickle() 
