# Overview

This project provides a python scripts to simulate the wave propagation from 
an open water body into an aquifer. Simulations can be performed via the numerical
solver ModFlow throught the Python API flopy or using a semianalytical solution.
The code comes along three input waves: a tidal wave, an artifical square wave
and the time series of a river.

## Structure

The project is organized as follows:

- `README.md`   - description of the project
- `LICENSE`     - the default license is MIT
- `results/`    - folder with simulation results and plots
- `data/`       - input data of waves in open water bodies
- `src/`        - folder containing the Python scripts of the project:
  + `01_plot_waves.py` 
                 - Settings of three example, including the plotting of
                   the three input waves
  + `02_run_num_model.py` 
                 - running numerical simulation of selected example and saving
                   simulation results in pickle-file (in results)
  + `03_run_ana_model.py` 
                 - running analytical model of selected example and saving
                   simulation results in pickle-file (in results)
  + `04_compare_models.py` 
                 - compare numerical and analytical model results for 
                   selected example 
  + `05_check_num_river.py` 
                 - comparing river example at two temporal resolutions of 
                   input wave                   
  + `WavePropagationAquifers.py` 
                 - class of 1D aquifer model containing routines for numerical
                   simulation, input wave handling, fft decomposition and 
                   analytical solution
  + `flopy_model.py` 
                 - functions of 1D Modflow models (via flopy) for simulating
                   wave propagation in aquifers from time dependent boundary condition
                       
## Python environment

To make the example reproducible, we provide the following files:
- `requirements.txt` - requirements for [pip](https://pip.pypa.io/en/stable/user_guide/#requirements-files) to install all needed packages

## Workflow

After finalizing your work, you should tag the repository with a version like `v1.0`.

Then, a [Zenodo](https://zenodo.org/) release will be created, so you can cite the repository in you publication.

Please keep your `master` branch in line with the latest release.
For further development use the `develop` branch and update `master` with pull-requests.


## Contact

You can contact us via <a.zech@uu.nl>.


## License

MIT © 2021
