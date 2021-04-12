# -*- coding: utf-8 -*-
"""
@copywrite, Peter Kner, University of Georgia, 2019
"""

import circle as cz
import pylab
import numpy as np


def getcross(nz=16,nx=64):
    ''' look at cross terms of zernikes '''
    zarr = np.zeros((nz,nx,nx))
    for ii in range(nz):
        zarr[ii] = cz.zernike_noll(ii+1,nx)
    cmat = np.zeros((nz,nz))
    for ii in range(nz):
        for jj in range(nz):
            w1 = zarr[ii]
            w2 = zarr[jj]
            cmat[ii,jj] = ( w1 * w2.conj() ).sum() / ( w2 * w2.conj() ).sum()
    pylab.imshow(cmat, vmin=0, vmax=0.1)
    return True

def get_gs(nz=16, nx=64):
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
    cmat = np.zeros((nz, nz))
    for ii in range(nz):
        for jj in range(nz):
            w1 = phi[ii]
            w2 = phi[jj]
            cmat[ii, jj] = (w1 * w2.conj()).sum() / (w2 * w2.conj()).sum()
    pylab.imshow(cmat, vmin=0, vmax=0.1)
    return phi
