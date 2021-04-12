# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 09:58:54 2020

@author: linruizhe
"""

import circle as cz
import pylab
import numpy as np

def get_gs(nz=64, nx=16, verbose=False):
    ''' gram-schmidt orthogonalization '''
    zarr = np.zeros((nz, nx, nx))
    for ii in range(nz):
        zarr[ii] = cz.zernike_noll(ii + 1, nx)
    phi = np.zeros((nz, nx, nx))
    # orthogonalize
    for ii in range(nz):
        phi[ii] = zarr[ii]
        for jj in range(ii):
            c = (zarr[ii]*phi[jj]).sum()
            phi[ii] = phi[ii] - c*phi[jj]
        phi[ii] = phi[ii]/np.sqrt((phi[ii]**2).sum())
    # look at cross terms
    if (verbose):
        cmat = np.zeros((nz, nz))
        for ii in range(nz):
            for jj in range(nz):
                w1 = phi[ii]
                w2 = phi[jj]
                cmat[ii, jj] = (w1 * w2.conj()).sum() / (w2 * w2.conj()).sum()
        pylab.imshow(cmat, vmin=0, vmax=0.1)
    return phi

def decomwf(wf, zn):
    """
     Decompose the wavefront into Zernike polynomials.
     Args:
        wf (ndarray): The wavefront to be decomposed
        zn (ndarray): The number of zernike modes
     Returns:
        ndarray: The Zernike mode coefficients
     """
    nx, ny = wf.shape
    a = np.zeros(zn)
    zern = get_gs(zn,nx)
    for i in range(zn):
        wz = zern[i]
        a[i] = ( wf * wz.conj() ).sum() / ( wz * wz.conj() ).sum()
    return a

def recomwf(a, d):
    """
     Re-compose the wavefront from Zernike polynomials.
     Args:
        zn (ndarray): The number of zernike modes
        a (ndarray): The Zernike mode coefficients
     Returns:
        ndarray: The wavefront
     """
    zn = a.shape[0]
    wf = np.zeros((d,d))
    zern = get_gs(zn,d)
    for i in range(zn):
        wf += a[i] * zern[i]
    return wf