# -*- coding: utf-8 -*-
"""
metric values calculation
11/06/2018 Ruizhe Lin
"""

import numpy as np
# import Utility as U
import circle as c
import scipy.ndimage.filters as filters
from scipy.fftpack import fft2, fftshift
from scipy import signal

na = 1.2

def rms(img):
    nx,ny = img.shape
    n = nx*ny
    m = np.mean(img, dtype=np.float64)
    a = (img-m)**2
    r = np.sqrt(np.sum(a)/n)
    return r
    
def sharp(img):
    nx,ny = img.shape
    n = nx*ny
    m = img.min()
    img[img<=m] = 0.
    img[img>m] = img[img>m] - m
    a = np.sum(img*img).astype(float)
    b = ((np.sum(img))**2).astype(float)
    s = -(n*a)/(b)
    return s

def std(img):
    p = np.std(img, dtype=np.float64)
    q = np.mean(img, dtype=np.float64)
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
    s = np.sum(np.multiply(data,maxima))
    m = np.mean(data)
    p = s / m
    return p

def peakv(img):
    nx,ny = img.shape
    w = c.circle(64,nx)
    # w = U.discArray((nx,ny),64)
    img = img*w
    maxv = img.max()
    return maxv

def snr(wl,img,lpr,hpr):
    nx,ny = img.shape
    m = img.min()
    img[img<=m] = 0.
    img[img>m] = img[img>m] - m
#    w = U.discArray((nx,ny),256)
#    img = img*w
    na = 1.2
    dp = 1/(nx*0.089)
    radius = (na/wl)/dp
    msk = c.circle(2.0*radius,nx)
    siglp = radius* 0.1
    sighp = radius * 0.4
    lp=gaussianArr(shape=(nx,nx), sigma=siglp, peakVal=1, orig=None, dtype=np.float32)
    hp=1-gaussianArr(shape=(nx,nx), sigma=sighp, peakVal=1, orig=None, dtype=np.float32)
    aft = fftshift(fft2(img))
    hpC = (np.abs(hp*aft*msk)).sum()
    lpC = (np.abs(lp*aft*msk)).sum()
    res = hpC/lpC
    return res
    
# def snr(wl,img,siglp,sighp):
#     img = np.array(img)
#     nx = 256
#     img = img[128:384,128:384]
#     nx,ny = img.shape
#     w = c.circle(100,nx)
#     # w = U.discArray((nx,ny),100)
#     img = img*w
#     na = 1.2
#     dp = 1/(nx*.089)
#     radius = (na/wl)/dp
#     msk = c.circle(2.0*radius,nx)
#     # msk = U.discArray((nx,nx),2.0*radius)# - U.discArray((nx,nx),0.05*radius)
#     lp=gaussianArr(shape=(nx,nx), sigma=siglp, peakVal=1, orig=None, dtype=np.float32)
#     hp=1-gaussianArr(shape=(nx,nx), sigma=sighp, peakVal=1, orig=None, dtype=np.float32)
#     aft = fftshift(fft2(img))
#     aft = aft*msk
#     hpC = (hp*np.abs(aft)).sum()
#     lpC = (lp*np.abs(aft)).sum()
#     res = hpC/lpC
#     return res

def hf(wl,img):
    nx,ny = img.shape
#    m = img.min()
#    img[img<=m] = 0.
#    img[img>m] = img[img>m] - m
    # w = c.circle(256,nx)
    # w = U.discArray((nx,ny),256)
    # img = img*w
    na = 1.2
    dp = 1/(nx*0.089)
    radius = (na/wl)/dp
    msk = c.circle(1.8*radius,nx) - c.circle(0.1*radius,nx)
    # msk = U.discArray((nx,nx),1.8*radius) - U.discArray((nx,nx),0.1*radius)
    lp=gaussianArr(shape=(nx,nx), sigma=0.8*radius, peakVal=1, orig=None, dtype=np.float32)
    hp=1-gaussianArr(shape=(nx,nx), sigma=0.8*radius, peakVal=1, orig=None, dtype=np.float32)
    aft = fftshift(fft2(img))
    aft = aft*(msk*lp*hp)
    res = (np.abs(aft)).sum()
    return res

def gaussianArr(shape,sigma,peakVal,orig=None,dtype=np.float32):
    nx,ny = shape
    if orig == None:
        ux = nx/2.
        uy = ny/2.
    else:
        ux,uy = orig
    return np.fromfunction(lambda i,j: np.exp(-((i-ux)**2.+(j-uy)**2.)/(2.*sigma**2.)),(nx,ny),dtype=dtype)*peakVal#/(sigma*sqrt(2*np.pi))
