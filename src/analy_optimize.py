#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 21:37:43 2021

@author: zech0001
"""

import numpy as np

a = np.arange(10)
b = np.arange(15)
c = np.arange(20)

na = len(a)
nb = len(b)
nc = len(c)

# a1 = np.tile(a,(nb,nc,1))
# a2 = np.moveaxis(a1,-1,0)

# b1 = np.tile(b,(na,nc,1))
# b2 = np.moveaxis(b1,-1,1)

# c1 = np.tile(c,(nb,na,1))
# c2 = np.moveaxis(c1,0,1)

a1 = np.tile(a,(nb,1))


        else:

            # a1 = np.tile(a,(nb,nc,1))
            # a2 = np.moveaxis(a1,-1,0)
            
            # b1 = np.tile(b,(na,nc,1))
            # b2 = np.moveaxis(b1,-1,1)
            
            # c1 = np.tile(c,(nb,na,1))
            # c2 = np.moveaxis(c1,0,1)

            x_term = np.tile(x_ana,(n_fft,1)) 

            amp_tile = np.tile(self.fft['amplitude'],(nt,1))
            wl_tile = np.tile(self.fft['wavelength'],(nt,1))
            ps_tile = np.tile(self.fft['phase_shift'],(nt,1))
            t_tile = np.tile(t_ana,(n_fft,1)).T
            x_tile = np.tile(x_ana,(n_fft,1))                        
            wl_term = np.sqrt(wl_tile * 0.5*self.ss/self.hk)

            corr = A[0] - A[0] * erfc(np.tile(x_ana,(nt,1))/(2*np.sqrt(np.tile(t_ana,(nx,1)).T * self.hk/self.ss)))

            time_curve = amp_tile * np.exp(x_tile * wl_term) * np.cos(t_tile*wl_tile + ps_tile -x_tile * wl_term)

            h = time_curve.sum(axis=0)  - corr
            
# def model(t, x, A, w, phi, headstart, S_K):
#     a = -x * np.sqrt(w * S_K/2 )
#     corr = A[0] - A[0] * erfc(x/(2*np.sqrt(t / S_K)))
    
#     A = np.array([A,])
#     a = np.array([a,])
#     w = np.array([w,])
#     phi = np.array([phi,])
    
#     tA = np.transpose(A)
#     ta = np.transpose(a)
#     tw = np.transpose(w)
#     tphi = np.transpose(phi)
    
#     time_curve = np.repeat((tA * np.exp(ta)).reshape(len(fft_wave)),len(t)).reshape(len(fft_wave),len(t)) * np.sin( t * tw + tphi + ta)
#     time_curve = time_curve.sum(axis=0) + headstart - corr
#     return (time_curve)





