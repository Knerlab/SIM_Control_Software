# -*- coding: utf-8 -*-
"""
Created on Mon Apr 09 16:36:51 2012

@author: kner
"""

import numpy as N
import zernike as Z

def discArray(shape=(128,128),radius=64,origin=None,dtype=N.float64):
    nx = shape[0]
    ny = shape[1]
    ox = nx/2
    oy = ny/2
    x = N.linspace(-ox,nx-ox,nx)
    y = N.linspace(-oy,ny-oy,ny)
    X,Y = N.meshgrid(x,y)
    rho = N.sqrt(X**2 + Y**2)
    disc = (rho<radius).astype(dtype)
    if not origin==None:
        s0 = origin[0]-int(nx/2)
        s1 = origin[1]-int(ny/2)
        disc = N.roll(N.roll(disc,s0,0),s1,1)
    return disc
    
def radialArray(shape=(128,128), func=None, origin=None, dtype=N.float64):
    nx = shape[0]
    ny = shape[1]
    ox = nx/2
    oy = ny/2
    x = N.linspace(-ox,nx-ox,nx)
    y = N.linspace(-oy,ny-oy,ny)
    X,Y = N.meshgrid(x,y)
    rho = N.sqrt(X**2 + Y**2)
    rarr = func(rho)
    if not origin==None:
        s0 = origin[0]-nx/2
        s1 = origin[1]-ny/2
        rarr = N.roll(N.roll(rarr,s0,0),s1,1)
    return rarr
    
def shift(arr,shifts=None):
    if shifts == None:
        shifts = N.array(arr.shape)/2
    if len(arr.shape)==len(shifts):
        for m,p in enumerate(shifts):
            arr = N.roll(arr,p,m)
    return arr
    
def buildphiZ(phi,shape,radius):
    nx = shape[0]
    pupil = N.zeros(shape)
    for m,amp in enumerate(phi):
        pupil = pupil + amp*Z.Zm(m,radius,None,nx)
    return pupil