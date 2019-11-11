# -*- coding: utf-8 -*-
"""
Metric values calculation
@copywrite, Ruizhe Lin and Peter Kner, University of Georgia, 2019
"""

import numpy as N
import Utility as U
import scipy.ndimage.filters as filters
from scipy.fftpack import fft2, fftshift
from scipy import signal

def rms(img):
    nx,ny = img.shape
    n = nx*ny
    m = N.mean(img, dtype=N.float64)
    a = (img-m)**2
    r = N.sqrt(N.sum(a)/n)
    return r
    
def sharp(img):
    nx,ny = img.shape
    n = nx*ny
    img = (img-img.min()).astype(N.float64)
    a = N.sum(img*img)
    b = (N.sum(img))**2
    s = (n*a)/(b)
    return s

def std(img):
    p = N.std(img, dtype=N.float64)
    q = N.mean(img, dtype=N.float64)
    s = p/q
    return s

def peak(data):
    neighborhood_size = 6
    threshold = 45
    data_max = filters.maximum_filter(data, neighborhood_size)
    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0
    s = N.sum(N.multiply(data,maxima))
    m = N.mean(data)
    p = s / m
    return p

def snr(img):
    img = N.array(img)
    nx,ny = img.shape
    w = window(img)
    img = img*w    
    wl = 0.515 #0.67
    na = 1.2
    siglp = 4
    sighp = 64
    dp = 1/(nx*.089)
    radius = (na/wl)/dp
    msk = U.discArray((nx,nx),radius)
    lp=gaussianArr(shape=(nx,nx), sigma=siglp, peakVal=1, orig=None, dtype=N.float32)
    hp=1-gaussianArr(shape=(nx,nx), sigma=sighp, peakVal=1, orig=None, dtype=N.float32)
    aft = fftshift(fft2(img))
    aft = aft*msk
    hpC = (hp*N.absolute(aft)).sum()
    lpC = (lp*N.absolute(aft)).sum()
    res = hpC/lpC
    return res

def window(img):
    nx,ny = img.shape
    wx = signal.tukey(nx)
#    wy = signal.tukey(nx)
    winx = N.tile(wx,(ny,1))
    winy = winx.swapaxes(0,1)
    win = winx * winy
    return win

def gaussianArr(shape,sigma,peakVal,orig=None,dtype=N.float32):
    nx,ny = shape
    if orig == None:
        ux = nx/2.
        uy = ny/2.
    else:
        ux,uy = orig
    return N.fromfunction(lambda i,j: N.exp(-((i-ux)**2.+(j-uy)**2.)/(2.*sigma**2.)),(nx,ny),dtype=dtype)*peakVal#/(sigma*sqrt(2*N.pi))
