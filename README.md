[![DOI](https://zenodo.org/badge/389735344.svg)](https://zenodo.org/badge/latestdoi/389735344)

# Overview

This project provides a python scripts to simulate the wave propagation from 
an open water body into an 
aquifer. Simulations can be performed via the numerical
solver ModFlow throught the Python API flopy or using a semi-analytical solution.
The code comes along three input waves: a tidal wave, an artifical square wave
and the time series of a river.

## Structure

The project is organized as follows:

- `README.md`   - description of the project
- `LICENSE`     - the default license is MIT
- `results/`    - folder with simulation results and plots
- `data/`       - input data of waves in open water bodies and summarized 
                  numerical simulation results
- `src/`        - folder containing the Python scripts of the project:
  + `WavePropagationAquifers.py` 
                 - central class file containign class "WavePropagationAquifers" 
                   for 1D aquifer model (confined, leaky) with all characterizing
                   properties and routines for input wave handling, fft wave 
                   decomposition and analytical solution calculation, inverse 
                   estimation and running numerical simulation 
  + `flopy_model.py` 
                 - functions to start 1D Modflow models (via flopy) for simulating
                   wave propagation in aquifers with time dependent boundary 
                   condition under various aquifer conditions (confined, leaky, 
                   barrier)
  + `01_run_num_model.py` 
                - Settings of three examples, including the command for starting
                  the numerical simulation with the selected aquifer setting and
                  boundary condition; simulation results are saved as pickle-files
                  in the results directory
                  saved results are used for plot preparation (Fig 3 and S3)
  + `02_impact_observation_loc.py` 
                - script to process numerical simuation data for studying impact
                  of observation location on inverse estimation procedure
                  saved results are used for plot preparation of Fig 5 
  + `03_run_num_resistance.py` 
                - script to prepare inverse estimation results for various values
                  of resistances for the aquifer setting of leakage and barrier
                  and save results for files:
                    -  running numerical simulations for range of resistance values
                    -  postprocess simuation data by selecting head values at
                       specified spatial locations
                    -  perform inverse estimation of diffusivity from 
                       numerical results with confined aquifer solution
                    -  perform inverse estimation of diffusivity from 
                       numerical results with leakage aquifer solution
                  saved results are used for plot preparation Fig 6 & 7
  + `F02_Input_BC_Waves.py` 
                - reproducing Figure 2 of the manuscript: settings of three 
                  example, including plotting of input waves
  + `F03_Solution_Performance.py` 
                - reproducing Figure 3 of the manuscript: comparing analytical 
                  solution to numerical simulation results from pickle-file 
  + `F04_Fit_Confined.py` 
                - reproducing Figure 4 of the manuscript: fitting analytical 
                  solution to numerical simulation results from pickle-file 
  + `F05_Impact_Location.py` 
                - reproducing Figure 5 of the manuscript: values of relative 
                  differences of diffusivity estimate as function of piezometer
                  location (distance x to open water body)
                  using results prepared with 02_impact_observation_loc.py
  + `F06_Impact_Resistance.py` 
                - reproducing Figure 6 of the manuscript: relative difference 
                  for diffusivity estimates when fitting the confined aquifer 
                  solution to complex aquifer settings with leakage and flow 
                  barrier as function of the resistance value
                  using results prepared with 03_run_num_resistance.py
  + `F07_Fit_Leakage.py` 
                - reproducing Figure 7 of the manuscript: relative difference 
                  for diffusivity and resistance factor estimates when fitting 
                  the analytical solution to numerical simulation results for 
                  the leaky aquifer setting as function of the confining layer
                  resistance value
                  using results prepared with 03_run_num_resistance.py
  + `SF01_Input_Wave_Reconstruction.py`
                - reproducing Figure S1 of the supporting information:     
                  comparing input BC time series to FFT reconstructions
  + `SF02_Model_Fitting_Max_Wave`
                - reproducing Figure S2 of the supporting information: fitting 
                  analytical solution to numerical simulation results from 
                  pickle-file (similar to Fig 4) with fit of dominant wave
                  component                          

## Python environment

To make the example reproducible, we provide the following files:
- `requirements.txt` - requirements for [pip](https://pip.pypa.io/en/stable/user_guide/#requirements-files) to install all needed packages


## Contact

You can contact us via <a.zech@uu.nl>.


## License

MIT © 2021
