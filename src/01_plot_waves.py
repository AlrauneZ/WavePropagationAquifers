import numpy as np
import matplotlib.pyplot as plt
from WavePropagationAquifers import WavePropagationAquifers

###############################################################################
### Select Model settings
###############################################################################

BC1 = dict(
    task_name =  'tide_wave',    
    cut_input = 1310,
    # normalize_time = False,
    normalize_time = True,
    # normalize_hs = False,
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


# Ex1 = WavePropagationAquifers(**BC1)
# Ex1.read_wave(**BC1)

# Ex2 = WavePropagationAquifers(**BC2)
# Ex2.read_wave(**BC2)

Ex3 = WavePropagationAquifers(**BC3)
Ex3.read_wave(**BC3)

###############################################################################
### Plot Results
###############################################################################

plt.close('all')

# Ex1.plot_wave()
# Ex2.plot_wave()
Ex3.plot_wave()

# Ex1.reconstruct_wave_fft()
# fft = Ex.fft
# plt.plot(Ex1.wave_time,Ex1.wave_reconst,label = 'wave construction fron FFT')
# plt.legend()
