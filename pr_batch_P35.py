"""
@copywrite, Peter Kner, University of Georgia, 2019
"""

import numpy as N
import matplotlib.pyplot as P
import tifffile as T

class Smat(object):

    def __init__(self,stack,radius=None):
        self.wfstack = stack.copy()
        nz,nx,ny = self.wfstack.shape
        if not(nz == 69):
            raise 'There must be 69 images!'
        # set size of array
        if radius==None:
            radius = 15
        self.sz= 2*radius
        # cut out wavefront
#        R0 = self.wfstack[:,(nx/2-radius):(nx/2+radius),(nx/2-radius):(nx/2+radius)]
        R0 = self.wfstack
        msk = (R0[0]!=0.0).astype(N.float32)
        nz, nx, ny = R0.shape
        # remove mean
        for m in range(nz):
            mn = R0[m].sum()/msk.sum()
            R0[m] = msk*(R0[m]-mn)
        R = N.zeros((nx*ny,nz))
        for m in range(nz):
            R[:,m] = R0[m].reshape(nx*ny)
        self.R = R
        self.calcS(21)
        #Y.view(R0)

    def __del__(self):
        pass

    def calcS(self,Ns=21): # SVD decomposition of R. Use largest Ns vectors.
        ''' Ns is number of singular values to retain '''
        u,s,vh = N.linalg.svd(self.R,0,1)
        ut = N.transpose(u)
        self.s = s
        sd = 1./s
        s = N.zeros((69,69))
        for i in range(Ns):
            s[i,i] = sd[i]
        v = N.transpose(vh)
        t = N.dot(v,s)
        self.u = u
        self.v = v
        self.S = N.dot(t,ut)
        return True

    def calcSalpha(self,Ns=21,alpha=2.0): # SVD decomposition of R. Use largest Ns vectors.
        ''' Ns is number of singular values to retain '''
        u,s,vh = N.linalg.svd(self.R,0,1)
        ut = N.transpose(u)
        self.s = s
        sd = 1./N.sqrt(s**2 + alpha**2)
        s = N.zeros((69,69))
        for i in range(Ns):
            s[i,i] = sd[i]
        v = N.transpose(vh)
        t = N.dot(v,s)
        self.u = u
        self.v = v
        self.S = N.dot(t,ut)
        return True

    def view_eigen_vectors(self):
        t = (self.u.reshape(self.sz,self.sz,69)).swapaxes(1,2).swapaxes(0,1)
        T.imshow(t,vmin=None)
        P.plot(self.s)
        return True
