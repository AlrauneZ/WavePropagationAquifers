import numpy as np
from scipy.special import erfc
# import copy
import os #, sys, subprocess, shutil
import pickle
# import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')

import flopy_model


class WavePropagationAquifers:
    
    def __init__(self,
                 task_name = 'wave',
                 task_root = '../results' ,
                 hk  = 25,
                 ss = 1E-5,               
                 **settings,
                 ):
    
        self.task_name = task_name
        self.task_root = task_root

        self.hk = hk  # horizontal hydraulic conductivity
        self.ss = ss  # storage term
        
        # self.settings = copy.copy(DEF_settings)
        # self.settings.update(**settings)

        self.wave = None
        self.fft = None

    def read_wave(self,
                  file_wave = None,
                  cut_input = False,
                  normalize_time = True,
                  coarsen_time = 0,
                  normalize_hs = True,
                  **kwargs,
                  ):

        """ read in complex wave data from file and preprocess
        """
        
        if file_wave == None:
            file_wave = '../data/{}.txt'.format(self.task_name)
            # file_wave = '../data/{}.csv'.format(self.task_name)
        if not os.path.isfile(file_wave):
            raise ValueError("wave input file is not valid: " + str(file_wave))

        tide = np.loadtxt(file_wave)#,delimiter=',')
        if cut_input:
            tide = tide[:int(cut_input),:]

        self.dt = (tide[-1,0]-tide[0,0])/len(tide[:,0])
        
        if normalize_time is False:
            self.wave_time = tide[:,0]
            self.wave = tide[:,1]
        else:
            wave_time = tide[:,0]
            wave = tide[:,1]

            if coarsen_time>0:
                self.wave_time = np.arange(wave_time[0],wave_time[-1]+self.dt,coarsen_time*self.dt)
                self.dt *= coarsen_time
            else:
                self.wave_time = np.arange(wave_time[0],wave_time[-1]+self.dt,self.dt)
            self.wave = np.interp(self.wave_time,wave_time,wave)

        if normalize_hs:
            if normalize_time:
                hs0 = np.mean(self.wave)
            else:
                hs0 = np.average(self.wave[:-1], weights=np.diff(self.wave_time))
            self.wave -= hs0

        return self.wave_time,self.wave
    
    def plot_wave(self,
                  save = False,
                  **kwargs,
                  ):
        
        plt.plot(self.wave_time,self.wave,marker = '.',ls = '--',label='input {}'.format(self.task_name))        
        plt.grid(True)
        plt.xlabel("Time $t$")
        plt.ylabel("Water table $h(t)$")
        plt.legend()

        if save:
            plt.savefig('{}.png'.format(self.task_root), dpi=300)
            
    def run_numerical_model(self,
                            dir_sim = None,
                            write_to_file = True,
                            **kwargs,
                            ):

        if self.wave is None:
            raise ValueError('Input wave not give. Read in wave data first.')

        if dir_sim is None:
            dir_sim ='{}/NumModel_{}/'.format(self.task_root, self.task_name) 
        try:
            os.mkdir(dir_sim)
        except:
            pass
        print("Run ModFlow Model \n for boundary condition: {}\n in directory: {} \n".format(self.task_name,dir_sim))

        times_num,x_num,head_tx_num = flopy_model.flopy_model(
            dt=self.dt,
            wave=self.wave,    
            task_name = self.task_name,
            task_root = dir_sim,
            hk = self.hk,
            ss = self.ss,
            ) 
        
        self.times_num = np.array(times_num)
        self.x_num = np.array(x_num)
        self.head_tx_num = np.array(head_tx_num)
        
        if write_to_file:
            self.write_num_to_pickle()

        return self.times_num,self.x_num,self.head_tx_num
        # return self.head_tx_num

    def write_num_to_pickle(self,
                            file_results = None,                            
                            **kwargs,
                            ):

        if file_results is None:
             file_results= r'{}/{}_numerical.p'.format(self.task_root,self.task_name)

        # Write data in binary format
        with open(file_results, 'wb') as filehandle:
            # store the data as binary data stream
            pickle.dump([self.times_num,self.x_num,self.head_tx_num], filehandle)
        print("\nSimulation results saved to binary (pickle) file: \n {}".format(file_results))

    def read_num_from_pickle(self,
                             file_results = None,                            
                             **kwargs,
                             ):
        
        if file_results is None:
            file_results= r'{}/{}_numerical.p'.format(self.task_root,self.task_name)
            
        if not os.path.isfile(file_results):
            raise ValueError('File with stored numerical model results not existant at: \n {}'.format(file_results))

        with open(file_results, 'rb') as filehandle:
            data = pickle.load(filehandle)

        self.x_num = data[1]
        self.times_num = data[0]
        self.head_tx_num = data[2]
        # self.x_num = data[0]
        # self.times_num = data[1]
        # self.head_tx_num = data[2].T
        # read the data as binary data stream
        
        return self.times_num,self.x_num,self.head_tx_num
            

    def decompose_wave_fft(self,
                           amplitude_threshold= 0.0,
                           **kwargs,
                           ):
        #################################################################################
        # Decompose tidal wave
        #################################################################################
        fft = np.fft.fft(self.wave)
        
        nt = len(self.wave_time)
        freqs = np.fft.fftfreq(len(self.wave_time),self.dt)
        
        if (nt % 2) == 0:     # Check if n_obs is even or odd
            middle = nt//2
        else:
            middle = nt//2 + 1
        
        wavelength   = freqs[0:middle] * 2.*np.pi
        phase_shift  = np.arctan2(fft.imag[0:middle],fft.real[0:middle])
        # phase_shift  = np.arctan2(fft.imag[0:middle],fft.real[0:middle]) + np.pi/2
        amplitude    = np.zeros(middle)
        
        for i in range(middle):
            if abs(fft[i])/(nt) > amplitude_threshold:
                if i == 0:
                    coeff = 2
                else:
                    coeff = 1
                amplitude[i] = 1/(nt*coeff/2.) * np.sqrt(fft.real[i]**2 + fft.imag[i]**2)
    
        self.fft = dict(
            wavelength = wavelength,
            phase_shift = phase_shift,
            amplitude = amplitude,
            )
 
        return self.fft
 
    def reconstruct_wave_fft(self,
                             **kwargs,
                             ):
    
        #################################################################################
        # Reconstruct tidal wave and compare with the original wave
        #################################################################################

        if self.fft is None: 
            self.decompose_wave_fft(**kwargs)
       
        nt = len(self.wave_time)
        n_fft = len(self.fft['amplitude'])

        amp = np.tile(self.fft['amplitude'],(nt,1))
        wl = np.tile(self.fft['wavelength'],(nt,1))
        ps = np.tile(self.fft['phase_shift'],(nt,1))
        wt = np.tile(self.wave_time,(n_fft,1)).T
 
        self.wave_reconst = np.sum(amp*np.cos(wl*wt + ps),axis = 1)
 
        return self.wave_reconst

    def select_times(self,
                     nt = 15,
                     wout_version = True,
                     ):


        ### Wout's way (but error prone --> replace by new way after result checking)
        if wout_version:
            index = np.ones(nt)        
            for i in range(nt):
                index[i] = (i) * int(len(self.wave_time)/nt)
            self.t_indices = np.array(index,dtype = int)        

        else:
            ### improved way
            self.t_indices = np.linspace(0,len(self.wave_time)-1,nt,endpoint=True,dtype = int)
            # self.t_indices = np.linspace(0,len(self.wave_time),nt,endpoint=False,dtype = int)

        self.t_select = self.wave_time[self.t_indices]
        self.t_rel  = self.t_select/self.wave_time[-1]
 
        return self.t_indices,self.t_select,self.t_rel

    def run_analytical_model(self,
                                t_ana = 15,
                                x_ana = np.linspace(0,3000,20),   
                                write_to_file = True,
                                as_loop = False,
                                **kwargs,
                                ):

        if isinstance(t_ana,int):
            ### select t_ana time points within the wave_time (at equal distance)
            self.select_times(nt = t_ana)
            t_ana= self.t_select
        elif isinstance(t_ana,list):
            t_ana = np.array(t_ana)
        elif isinstance(t_ana,np.ndarray):
            pass
        else:
            raise ValueError('Input of time series not correct: integer, list or array')

        #################################################################################
        ### Check availability of wave decomposition components of input wave       
        #################################################################################
        if self.fft is None: 
            self.decompose_wave_fft(**kwargs)
        
        #################################################################################
        ### Set up analytical solution
        #################################################################################
        n_fft = len(self.fft['amplitude'])
        
        nx = len(x_ana) # number of spatial steps
        nt = len(t_ana)          
        
        # corr = np.ones(shape=(nx,nt))
        A = self.fft['amplitude']
        w = self.fft['wavelength']
        phi = self.fft['phase_shift']

        ### TODO: optimize calculation avoiding loops using matrix multiplication
        h = np.ones(shape=(n_fft,nt,nx))
        for i in range(0,n_fft):
            for j in range(0,nx):
                for k in range(0,nt):
                    if w[i] == 0:
                        h[i,k,j] = A[i] * erfc(x_ana[j] / (2*np.sqrt(self.hk/self.ss * t_ana[k])))
                    else:
                        h[i,k,j] = A[i] * np.exp(-x_ana[j] * np.sqrt(w[i]*self.ss / (2*self.hk))) * np.cos(w[i]*t_ana[k] + phi[i] - x_ana[j]*np.sqrt(w[i]*self.ss / (2*self.hk)))

        self.x_ana = x_ana
        self.times_ana = t_ana
        self.head_tx_ana = h.sum(axis=0)
 
        if write_to_file:
            self.write_ana_to_pickle()

        return self.times_ana,self.x_ana,self.head_tx_ana

    
    def write_ana_to_pickle(self,
                        file_results = None,                            
                        **kwargs,
                        ):

        if file_results is None:
             file_results= r'{}/{}_analytical.p'.format(self.task_root,self.task_name)

        # Write data in binary format
        with open(file_results, 'wb') as filehandle:
            # store the data as binary data stream
            pickle.dump([self.times_ana,self.x_ana,self.head_tx_ana], filehandle)
        print("\nAnalytical results saved to binary (pickle) file: \n {}".format(file_results))

    def read_ana_from_pickle(self,
                             file_results = None,                            
                             **kwargs,
                             ):
        
        if file_results is None:
            file_results= r'{}/{}_analytical.p'.format(self.task_root,self.task_name)

        if not os.path.isfile(file_results):
            raise ValueError('File with stored analytical model results not existant at: \n {}'.format(file_results))

        # read the data as binary data stream
        with open(file_results, 'rb') as filehandle:
            data = pickle.load(filehandle)
        self.x_ana = data[0]
        self.times_ana = data[1]
        self.head_tx_ana = data[2].T
        
        return self.times_ana,self.x_ana,self.head_tx_ana

def find_nearest_value(value, array):

    """Find match of given 'value' in given 'array'.

        Parameter
        ---------
        
        value       :   float    
                        value to find in array
        array       :   array
                        search for nearest match in this numpy array
        
        Return
        ------
        match       :   float
                        value in 'array' being closest to 'value'
    """


    array=np.array(array)

    index = (np.abs(array-value)).argmin()
    return array[index]

def find_nearest_index(value, array):

    """Find index of nearest match of given 'value' in given 'array'.

        Parameter
        ---------
        
        value       :   float    
                        index of this number will be found in array
        array       :   array
                        search for nearest match in this numpy array
        
        Return
        ------
        index       :   int
                        index of entry in 'array' being closest to 'value'
    """

    array=np.array(array)
    index = (np.abs(array-value)).argmin()

    return index
