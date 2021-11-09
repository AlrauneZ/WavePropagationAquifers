import numpy as np
from scipy.optimize import curve_fit
# from scipy.special import erfc
# import copy
# import lmfit
import os #, sys, subprocess, shutil
import pickle
import matplotlib.pyplot as plt
plt.close('all')

import flopy_model

class WavePropagationAquifers:
    
    def __init__(self,
                 task_name      = None,
                 task_root      = '../results' ,
                 BC_setting     = 'wave',
                 flow_setting   = 'confined',   
                 hk             = 25.,      # horizontal hydraulic conductivity (confined aquifer)
                 ss             = 5e-5,     # storage term (confined aquifer)
                 d_conf         = 10.,      # thickness (confined aquifer)
                 ss_unconfined  = 0.025,    # storage term (unconfined aquifer) = specific yield/thickness
                 c_L            = 100.,     # resistance = 1/leakage coefficient 
                 d_barrier      = 0.1,      # Thickness barrier [m]
                 c_barrier      = 100.,      # resistance barrier
                 **settings,
                 ):

        if task_name is None:
            self.task_name = '{}_{}'.format(BC_setting,flow_setting)
        else:
            self.task_name = task_name

        self.task_root = task_root
        self.flow_setting = flow_setting
        self.BC_setting = BC_setting

        ### confined aquifer parameters
        self.hk = hk  # horizontal hydraulic conductivity
        self.ss = ss  # storage term
        self.d_conf = d_conf # thickness of confined aquifer

        if flow_setting == 'leakage':
            ### parameters for unconfined aquifer and confining layer
            self.ss_unconfined = ss_unconfined  # storage term (unconfined aquifer)
            self.c_L = c_L      # resistance = 1/leakage coefficient= thickness confining/K_confining
        elif flow_setting == 'barrier':
            ### parameters for flow barrier (low-K interface)
            self.d_barrier = d_barrier  # Thickness barrier [m]
            self.c_barrier = c_barrier        # resistance barrier = thickness barrier/K_barrier
       
        self.wave = None
        self.fft = None
        self.head_tx_num = None

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
            file_wave = '../data/{}.txt'.format(self.BC_setting)
            # file_wave = '../data/{}.csv'.format(self.task_name)
        if not os.path.isfile(file_wave):
            raise ValueError("wave input file is not valid: " + str(file_wave))

        tide = np.loadtxt(file_wave)#,delimiter=',')
        if cut_input:
            tide = tide[:int(cut_input),:]

        self.dt = (tide[-1,0]-tide[0,0])/(len(tide[:,0])-1)
        # self.dt = (tide[-1,0]-tide[0,0])/(len(tide[:,0]))
 
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
               
    def run_numerical_model(self,
                            dir_sim = None,
                            write_to_file = True,
                            **kwargs,
                            ):

        if self.wave is None:
            raise ValueError('Input wave not given. Read in wave BC data first.')

        if dir_sim is None:
            dir_sim ='{}/NumModel_{}/'.format(self.task_root, self.task_name) 
        try:
            os.mkdir(dir_sim)
        except:
            pass
        print("Run ModFlow Model \n for boundary condition: {}\n in directory: {} \n".format(self.task_name,dir_sim))
                            
        if self.flow_setting == 'confined':   
            times_num,x_num,head_tx_num = flopy_model.model_confined(
    		    dt=self.dt,
    		    wave=self.wave,    
    		    task_name = self.task_name,
    		    task_root = dir_sim,
    		    hk = self.hk,
    		    ss = self.ss,
                ztop = self.d_conf,
                **kwargs
    		    ) 

        elif self.flow_setting == 'leakage':   
            times_num,x_num,head_tx_num = flopy_model.model_leakage(
                dt=self.dt,
                wave=self.wave,    
                task_name = self.task_name,
                task_root = dir_sim,
                hk_unconfined  = self.hk,
                hk_confining = 1./self.c_L,                
                hk_confined = self.hk,
                ss_confined  = self.ss, 
                ss_unconfined  = self.ss_unconfined, 
                thickness_lay3 = self.d_conf,
                **kwargs
                ) 

        elif self.flow_setting == 'barrier':   
            times_num,x_num,head_tx_num = flopy_model.model_barrier(
    		    dt=self.dt,
    		    wave=self.wave,    
    		    task_name = self.task_name,
    		    task_root = dir_sim,
    		    hk = self.hk,
    		    ss = self.ss,
                ztop = self.d_conf,
                d_barrier = self.d_barrier,    # Thickness sheetpiles [m]
                c_barrier = self.c_barrier,    # resistance barrier [d]
                **kwargs
                ) 
        
        self.t_num = np.array(times_num)-self.dt
        self.x_num = np.array(x_num)
        self.head_tx_num = np.array(head_tx_num)
        
        if write_to_file:
            self.write_num_to_pickle(**kwargs)

        return self.t_num,self.x_num,self.head_tx_num

    def write_num_to_pickle(self,
                            file_results = None,     
                            coarsen_x = False,
                            **kwargs,
                            ):

        if coarsen_x is False:
            x_num,head_tx_num = self.x_num,self.head_tx_num

        elif (isinstance(coarsen_x, float) or isinstance(coarsen_x, int)) and coarsen_x>0:
            ix_max = find_nearest_index(coarsen_x,np.diff(self.x_num))

            x_coarse = np.arange(self.x_num[0],self.x_num[ix_max],coarsen_x)
            x_num = np.concatenate((x_coarse,self.x_num[ix_max:]))


            head_tx_num = np.zeros((len(self.t_num),len(x_num)))
            head_tx_num[:,len(x_coarse):] = self.head_tx_num[:,ix_max:]
            for it in range(len(self.t_num)):
                head_tx_num[it,:len(x_coarse)] = np.interp(x_coarse,self.x_num[:ix_max],self.head_tx_num[it,:ix_max])
            
        elif isinstance(coarsen_x, np.ndarray) and len(coarsen_x)<len(self.x_num):
            x_num = coarsen_x
            head_tx_num = np.zeros((len(self.t_num),len(x_num)))
            for it in range(len(self.t_num)):
                head_tx_num[it,:] = np.interp(x_num,self.x_num,self.head_tx_num[it,:])

        else:
            raise ValueError('Specification of coarse_x not correct: either float or np.array')

        if file_results is None:
             file_results= r'{}/{}_numerical.p'.format(self.task_root,self.task_name)
           
        # Write data in binary format
        with open(file_results, 'wb') as filehandle:
            pickle.dump([self.t_num,x_num,head_tx_num], filehandle)
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

        self.t_num = data[0]
        self.x_num = data[1]
        self.head_tx_num = data[2]

        return self.t_num,self.x_num,self.head_tx_num          

    def decompose_wave_fft(self,
                           amplitude_threshold= 0.0,
                           wave_data = 'input_BC', #'piez'
                           **kwargs,
                           ):
        #################################################################################
        # Decompose tidal wave
        #################################################################################

        if wave_data == 'input_BC':
            if self.wave is None:
                raise ValueError('Input wave not given. Read in wave data first.')           
            wave = self.wave
            wave_time = self.wave_time
            dt = self.dt
        elif wave_data == 'piez':
            if self.h_piez is None:
                raise ValueError('Piezometer data not given. Read/prepare/set data first.')           
            wave_time = self.t_piez 
            wave = self.h_piez 
            dt = (wave_time[-1]-wave_time[0])/(len(wave_time)-1)
        elif wave_data is None:
            raise ValueError('Provide wave data for decomposition.')           
        else:
            wave_time = wave_data[:,0]
            wave = wave_data[:,1]
            dt = (wave_time[-1]-wave_time[0])/(len(wave_time)-1)

        fft = np.fft.fft(wave)
        
        nt = len(wave_time)
        freqs = np.fft.fftfreq(nt,dt)
        
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
    
        fft = dict(
            wavelength = wavelength,
            phase_shift = phase_shift,
            amplitude = amplitude,
            )
        
        if wave_data == 'input_BC':
            self.fft = fft
        elif wave_data == 'piez':
            self.fft_piez = fft
    
        return fft
 
    def reconstruct_wave_fft(self,
                             wave_data = 'input_BC', #'piez'
                             extend_time = False, # if time should be extended at margins, give value between [0,1] of extended time
                             **kwargs,
                             ):
    
        #################################################################################
        # Reconstruct tidal wave and compare with the original wave
        #################################################################################

        if wave_data == 'input_BC':
            if self.wave is None:
                raise ValueError('Input wave not given. Read in wave data first.')           
            wave_time = self.wave_time
        elif wave_data == 'piez':
            if self.h_piez is None:
                raise ValueError('Piezometer data not given. Read/prepare/set data first.')           
            wave_time = self.t_piez 
        elif wave_data is None:
            raise ValueError('Provide wave data for decomposition.')           
        else:
            wave_time = wave_data[:,0]

        fft = self.decompose_wave_fft(wave_data = wave_data,**kwargs)
       
        nt = len(wave_time)
        n_fft = len(fft['amplitude'])

        amp = np.tile(fft['amplitude'],(nt,1))
        wl = np.tile(fft['wavelength'],(nt,1))
        ps = np.tile(fft['phase_shift'],(nt,1))

        if extend_time>0:
            t_ext = extend_time*max(wave_time)
            wave_time = np.linspace(wave_time[0]-t_ext,wave_time[-1]+t_ext,nt)
        
        wave_reconst = np.sum(amp*np.cos(wl*np.tile(wave_time,(n_fft,1)).T + ps),axis = 1)
            
        return wave_time,wave_reconst

    def prepare_piezometric_data(self,
                             x = 400,
                             write_to_file = True,
                             file_piez = None,                            
                             **kwargs,
                             ):

        if self.head_tx_num is None:
            self.read_num_from_pickle(**kwargs)
                    
        index = find_nearest_index(x,self.x_num)        

        self.x_piez = self.x_num[index]
        self.t_piez = self.t_num
        self.h_piez = self.head_tx_num[:,index] 
               
        if write_to_file:
            if file_piez is None:
                file_piez= r'{}/{}_piez.txt'.format(self.task_root,self.task_name)
            # Write data in text format
            np.savetxt(file_piez, np.vstack((self.t_piez, self.h_piez)).T,delimiter=',',header='times , heads, at x = {}'.format(self.x_piez))
            print("\nPrepared simulation results to piezometric data, \nsaved to text file: \n {}".format(file_piez))
        
        return self.x_piez, self.t_piez, self.h_piez

    def read_piezometric_data(self,
                              x = None,
                              file_piez = None,                            
                              **kwargs,
                              ):
        if file_piez is None:
            file_piez= r'{}/{}_piez.txt'.format(self.task_root,self.task_name)
        if not os.path.isfile(file_piez):
            raise ValueError("piezometric input data file is not valid: " + str(file_piez))

        data = np.loadtxt(file_piez,delimiter=',',skiprows =1)
        self.t_piez = data[:,0]
        self.h_piez = data[:,1]

        if x is None:
            ### Extract distance from piezometer data file ()
            fp = open(file_piez)
            line = fp.readline()
            self.x_piez = np.float(line.split('=')[1])
        else:
            self.x_piez = x
        
        return self.x_piez, self.t_piez, self.h_piez

    def set_piezometric_data(self,
                             x,
                             time,
                             head):
        
        self.x_piez = x
        self.t_piez = time
        self.h_piez = head

        return self.x_piez, self.t_piez, self.h_piez

    def relative_sim_times(self):
        
        self.times_rel = self.t_num/(self.t_num[-1]-self.t_num[0])
        
        return self.times_rel

    def select_sim_times(self,
                     times_select = [0.3,0.6,0.9],
                     t_rel = True,
                     ):

        times_select = np.array(times_select)
        
        if t_rel is True:
            self.relative_sim_times()

        self.head_num_select = np.zeros((len(times_select),len(self.x_num)))
        for ix in range(len(self.x_num)):
            self.head_num_select[:,ix] = np.interp(times_select, self.times_rel,self.head_tx_num[:,ix])

        return self.head_num_select

    def run_analytical_model(self,
                                t_ana = [0.3,0.6,0.9],
                                t_rel = True,
                                x_ana = np.linspace(0,3000,20), 
                                diffusivity = None,
                                cS = None,
                                # cT = None,
                                write_to_file = False,
                                **kwargs,
                                ):

        if t_rel is True:
            self.t_ana = np.array(t_ana)* (self.wave_time[-1]-self.wave_time[0])
        else:
            self.t_ana = np.array(t_ana)

        self.x_ana = np.array(x_ana,ndmin = 1)

        if diffusivity is None:
            diffusivity = self.ss/self.hk
        else:
            diffusivity = np.float(diffusivity)
        
        #################################################################################
        ### Check availability of wave decomposition components of input wave       
        #################################################################################
        if self.fft is None: 
            self.decompose_wave_fft(**kwargs)
        
        #################################################################################
        ### Calculate analytical solution
        #################################################################################      
        A = self.fft['amplitude']
        w = self.fft['wavelength']
        phi = self.fft['phase_shift']
        
        if self.flow_setting == 'confined':
            term1 = self.x_ana[np.newaxis,:,np.newaxis]*np.sqrt(0.5*w[np.newaxis,:] * diffusivity)
            term2 = term1
        elif self.flow_setting == 'leakage':
            if cS is None:
                cS = self.c_L*self.ss*self.d_conf
            else:
                cS = np.float(cS)
            # if cT is None:
            #     cT = self.c_L*self.hk*self.d_conf
            # else:
            #     cT = np.float(cT)
            

            cT = cS / diffusivity
            p = np.sqrt(0.5*np.sqrt(1/cT**2 + (w[np.newaxis,:]*diffusivity)**2) + 0.5/cT)
            term1 = self.x_ana[np.newaxis,:,np.newaxis] * p
            term2 = self.x_ana[np.newaxis,:,np.newaxis] * (0.5*w[np.newaxis,:] * diffusivity)/p

            # p2 = np.sqrt(np.sqrt(1./(cS*w[np.newaxis,:])**2 + 1) + 1./(cS*w[np.newaxis,:]))
            # term1 = self.x_ana[np.newaxis,:,np.newaxis] * np.sqrt(0.5*w[np.newaxis,:]*diffusivity) * p2
            # term2 = self.x_ana[np.newaxis,:,np.newaxis] * np.sqrt(0.5*w[np.newaxis,:]*diffusivity) / p2

        else:
            raise ValueError('No analytical solution available for flow setting: {}'.format(self.flow_setting))
            
        h_xt = A[np.newaxis,:] * np.exp(-term1) * np.cos(w[np.newaxis,:]*self.t_ana[:,np.newaxis,np.newaxis] + phi[np.newaxis,:] - term2)        
        self.head_tx_ana = h_xt.sum(axis=2)
       
        if write_to_file:
            self.write_ana_to_pickle()

        return self.t_ana,self.x_ana,self.head_tx_ana

    def write_ana_to_pickle(self,
                        file_results = None,                            
                        **kwargs,
                        ):

        if file_results is None:
             file_results= r'{}/{}_analytical.p'.format(self.task_root,self.task_name)

        # Write data in binary format
        with open(file_results, 'wb') as filehandle:
            # store the data as binary data stream
            pickle.dump([self.t_ana,self.x_ana,self.head_tx_ana], filehandle)
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
        self.t_ana = data[1]
        self.head_tx_ana = data[2].T
        
        return self.t_ana,self.x_ana,self.head_tx_ana
       
    def extract_dominant_input_wave_component(self,
                               max_period = 50,
                               **kwargs):

        if self.fft is None:
            self.decompose_wave_fft(wave_data = 'input_BC',**kwargs)

        ### extract all fft components with periods above max_period
        condition = 2.*np.pi / self.fft['wavelength'] < max_period
        
        ### identify most prominent wave: that with maximum amplitude
        A_cut = np.compress(condition,self.fft['amplitude'])
        i_cut = np.argmax(A_cut)
        index_max = np.compress(condition,np.arange(len(self.fft['amplitude'])))[i_cut]

        self.A_max = self.fft['amplitude'][index_max]
        self.w_max = self.fft['wavelength'][index_max]
        self.phi_max = self.fft['phase_shift'][index_max]

        return self.A_max, self.w_max, self.phi_max

    def extract_piez_wave_component(self,
                                    reshape_time= False,
                                    **kwargs,
                                    ):

        self.extract_dominant_input_wave_component(**kwargs)

        ### redistribute observed values at identical time points as input wave 
        ### to allow proper wave decomposition
        if reshape_time is True:
            h_piez_reshape = np.interp(self.wave_time,self.t_piez,self.h_piez)
            self.h_piez = h_piez_reshape
            self.t_piez = self.wave_time
                             
        self.decompose_wave_fft(wave_data = 'piez', **kwargs)

        index = find_nearest_index(self.w_max,self.fft_piez['wavelength'])
        w_piez_max = self.fft_piez['wavelength'][index]
               
        if relative_difference(self.w_max,w_piez_max)>1:
            print('Warning: dominant wave component of observed data could not be matched well to dominant input wave component')
        
        A_piez_max = self.fft_piez['amplitude'][index]
        phi_piez_max = self.fft_piez['phase_shift'][index]
        self.wave_piez_max = A_piez_max * np.cos( w_piez_max * self.t_piez + phi_piez_max)
   
        return A_piez_max, w_piez_max, phi_piez_max

    def fit_data_to_analytical(self,                
                               **kwargs):

        print('\n################################')
        print("Run Inverse Parameter Estimation \n for boundary condition: {}\n and flow conditions: {} \n".format(self.BC_setting,self.flow_setting))
        
        if self.flow_setting == 'confined':
            results = self.fit_data_to_analytical_confined(**kwargs)
        elif self.flow_setting == 'leakage':
            results = self.fit_data_to_analytical_leakage(**kwargs)
        else:
            raise ValueError('No analytical solution available for flow setting: {}'.format(self.flow_setting))

        return results 
        
    def fit_data_to_analytical_confined(self,
                 p0_diff = 1E-6,
                 bounds_diff = (1E-9, 1E-3),
                 verbose = True,
                 **kwargs,
                 ):
       
        self.extract_dominant_input_wave_component()
        self.extract_piez_wave_component()
        
        def model_confined(t,diffusivity):
            ### simple wave solution
            a = self.x_piez * np.sqrt(self.w_max * 0.5 * diffusivity)           
            h_xt = self.A_max * np.exp(-a) * np.cos(self.w_max*t + self.phi_max - a)            
            return h_xt

        popt, pcov = curve_fit(
            model_confined,
            self.t_piez,
            self.wave_piez_max,
            p0 = p0_diff,
            bounds=bounds_diff,
            )

        self.diff_fit = popt[0]

        ### calculate relative difference of fitted diffusivity to input value
        self.eps_diff = relative_difference(self.ss / self.hk, self.diff_fit)

        ### calculate dominant wave componant of analytical solution for fitted diffusivity value
        self.wave_ana_max = model_confined(self.t_piez,self.diff_fit)

        ### calculate head with analytical solution for fitted diffusivity value
        if self.flow_setting == 'confined':
            self.head_ana_fit = self.run_analytical_model(
                t_ana = self.t_piez,
                t_rel = False,
                x_ana = self.x_piez,
                diffusivity = self.diff_fit,
                )[2]

        if verbose:
            print('Inverse Estimation Results')
            print('-------------------------- \n')
            print("Input Diffusivity S_s/K = {:.2e} d/m^2".format(self.ss/self.hk))
            print("Fitted diffusivity = {:.2e} d/m^2".format(self.diff_fit))
            print("Relative difference = {:.2f}%".format(self.eps_diff))
            
        return self.diff_fit, self.eps_diff #,self.head_ana_fit      

    def fit_data_to_analytical_leakage(self,
                 p0 = [2e-6,0.05],
                 bounds = ([1e-9,1e-3], [1e-3,10]),
                 verbose = True,
                 **kwargs,
                 ):
        
        # if self.flow_setting != 'leakage':
        #     raise ValueError('Fitting routine does not match flow setting! (choose leakage)')

        self.extract_dominant_input_wave_component()
        self.extract_piez_wave_component()

        def model_leakage(t,diffusivity,cS):
            ### leakage model for dominant wave component
            # p2 = np.sqrt(np.sqrt(1./cS**2 + 1) + 1./(cS*self.w_max))
            p2 = np.sqrt(np.sqrt(1./(cS*self.w_max)**2 + 1) + 1./(cS*self.w_max))
            a1 = self.x_piez *np.sqrt(self.w_max * 0.5 * diffusivity)* p2
            a2 = self.x_piez *np.sqrt(self.w_max * 0.5 * diffusivity) / p2
            h_xt = self.A_max * np.exp(-a1) * np.cos(self.w_max*t + self.phi_max -a2)
            return h_xt
                    
        popt, pcov = curve_fit(
            model_leakage,
            self.t_piez,
            self.wave_piez_max,
            p0 = p0,
            bounds=bounds,
            )

        self.diff_fit = popt[0]
        self.cS_fit = popt[1]

        ### calculate relative difference of fitted diffusivity to input value
        self.eps_diff = relative_difference(self.ss / self.hk, self.diff_fit)
        self.eps_cS = relative_difference(self.c_L*self.ss*self.d_conf, self.cS_fit)

        ### calculate dominant wave componant of analytical solution for fitted values
        self.wave_ana_max = model_leakage(self.t_piez,self.diff_fit,self.cS_fit)

        ### calculate head with analytical solution for fitted diffusivity value
        if self.flow_setting == 'confined':
            self.head_ana_fit = self.run_analytical_model(
                t_ana = self.t_piez,
                t_rel = False,
                x_ana = self.x_piez,
                diffusivity = self.diff_fit,
                cS = self.cS_fit,
                # cT = self.cT_fit,
                )[2]

        if verbose:
            print('Inverse Estimation Results')
            print('-------------------------- \n')
            print("Input Diffusivity S_s/K = {:.2e} d/m^2".format(self.ss/self.hk))
            print("Fitted diffusivity = {:.2e}".format(self.diff_fit))
            print("Relative difference = {:.2f}%".format(self.eps_diff))
            
            print("\nInput Factor c*S = {:.2e} d".format(self.c_L*self.ss*self.d_conf))
            print("Fitted factor c*S = {:.2e} d".format(self.cS_fit))
            print("Relative difference = {:.2f}%".format(self.eps_cS))
        
        return self.diff_fit,self.cS_fit, self.eps_diff, self.eps_cS


###############################################################################
### General/Auxiliary functions 
###############################################################################
def save_heads_cL(h_piez,
                  x_piez,
                  t_piez,
                  cL_range,
                  file_results_cL='./heads_cL.txt'):
    
    ### save postprocessed numerical simulation data at x_piez to one file
    data_save = np.zeros((h_piez.shape[0]+1,h_piez.shape[1]+1))
    data_save[1:,1:]    = h_piez
    data_save[0,0]      = x_piez
    data_save[1:,0]     = t_piez
    data_save[0,1:]     = cL_range

    np.savetxt(file_results_cL,data_save,delimiter =',')
    print('Simulated heads at specified location for all cL saved in file\n',file_results_cL)

def read_heads_cL(file_results_cL):
    
    '''
        h_piez,x_piez,t_piez,cL_range = read_head_cL(file_results_cL) 
    '''
    data_save = np.loadtxt(file_results_cL,delimiter =',')
    h_piez = data_save[1:,1:]
    x_piez = data_save[0,0]
    t_piez = data_save[1:,0]
    cL_range = data_save[0,1:]

    return h_piez,x_piez,t_piez,cL_range

def damping_coefficient(A_piez,A_wave):
    
    damp = 1- A_piez/A_wave
    
    return damp

def relative_difference(value_init, value_fit):

    """Calculate relative difference between two values.

        Parameter
        ---------
        value_init   :   float, reference value for error estimation 
        value_fit       :   float, values to identify relative difference to
        
        Return
        ------
        rel_diff       :   float between 0 and 100 relative difference between values in percent
    """
    
    rel_diff = abs((value_fit - value_init)/value_init) * 100    

    return rel_diff

def find_nearest_index(value, array):

    """Find index of nearest match of given 'value' in given 'array'.

        Parameter
        ---------
        value       :   float, index of this number will be found in array
        array       :   array, search for nearest match in this numpy array
        
        Return
        ------
        index       :   int, index of entry in 'array' being closest to 'value'
    """

    array=np.array(array)
    index = (np.abs(array-value)).argmin()

    return index
