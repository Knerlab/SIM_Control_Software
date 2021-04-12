# -*- coding: utf-8 -*-
"""
Calculates the amplitude and phase in the back pupil plane from a stack of
images of a point source.

Copyright (C) 2008-2015  
Peter Kner
University of California, San Francisco
University of Georgia
kner@engr.uga.edu

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

To cite this program please use
P. Kner, L. Winoto, D. A. Agard, and J. W. Sedat,
"Closed loop adaptive optics for microscopy without a wavefront sensor,"
Proc. SPIE 7570, pp. 6-9 (2010).

center of psf in (x,y) should be at (0,0)
Look at the function test() to see how to run the program

Created on Sat Sep 19 14:09:44 2015
"""

import numpy as np
import pylab
import Utility36 as U
import tifffile as tf
import VecZern_P33 as vz

pi = np.pi
fft2 = np.fft.fft2
ifft2 = np.fft.ifft2
fftshift = np.fft.fftshift

class phase(object):

    def __init__(self,fns=None):
        self.stacks = None
        self.ms = [3,2,1,0,10,4,5,6,7,8,9]
        self.av= [-0.04, -0.06, -0.08, -0.1, 0., 0.02, 0.04, 0.06, 0.08, 0.1, -0.02]
        if fns==None:
            raise Exception("must give list of files!")
        else:
            self.fns = fns
#        self.ps = {'Nx':128,'dx':0.0889,'dz':0.2,'n2':1.5,'na':1.3,'wl':0.515,\
#            'zcenter':3.0,'offset':100}
        self.Nx = 64
        self.dx = 0.089
        self.dz = 0.2
        self.n2 = 1.525
        self.na = 1.2
        self.wl = 0.58
        self.zcenter = 3.0
        self.offset = 100
        self.cps = [250,225,128] # cs,cy and size for cutting out bead

    def __del__(self):
        pass

    def readfiles(self,bgr=900):
        a = []
        for fn in self.fns:
            a.append(tf.imread(fn)-bgr)
        self.stacks = a
        return True

    def getphase(self,st,iters=128,cntr=False):
        if cntr:
            mx = st.max()
            nz,nx,ny = np.unravel_index(st.argmax(), st.shape)
            stack = self.cut(st,int(nx), int(ny),self.cps[2])
        else:
            stack = self.cut(st,self.cps[0],self.cps[1],self.cps[2])
        p = mre(stack)
        p.set_params(self.dx, self.na, self.wl)
        p.dz = self.dz
        p.zcenter = self.zcenter
        p.offset = self.offset
        p.run(iters)
        amp, phase = p.get_amp_phase()
        return amp, phase, p.getStrehl()

    def cut(self,st0,cx,cy,sz):
        b = int(sz/2)
        st1 = st0[:,(cx-b):(cx+b),(cy-b):(cy+b)]
        st2 = np.zeros(st1.shape)
        st2[:,:b,:b] = st1[:,b:,b:]
        st2[:,b:,b:] = st1[:,:b,:b]
        st2[:,b:,:b] = st1[:,:b,b:]
        st2[:,:b,b:] = st1[:,b:,:b]
        return st2

    def run(self,iters=128,cut=True):
        if self.stacks==None:
            self.readfiles()
        self.amp = []
        self.phase = []
        self.zarr = []
        for m,st in enumerate(self.stacks):
            print(m)
            amp, phase, strehl = self.getphase(st,iters,cut)
            self.phase.append(phase)
            self.amp.append(amp)
            #bp1,zar = mny.uw(mny.Cshift(bp),self.ps)
            #self.zarr.append(zar)
        return strehl

    def analyze(self,mcenter):
        wfc = []
        yc = []
        for bp in self.phase:
            bp0 = bp - self.phase[mcenter]
            wfc.append(bp0)
            #yc.append(bp0[:,64])
        #P.plot(yc)
        #self.yc = N.array(yc)
        self.wfc = np.array(wfc)
        return True

    def analyze2(self,mcenter,remttf=False):
        ''' this worked well - do zernike decomposition and rebuild phase to deal with
            phase wrapping and noise -- first 37 modes '''
        nx = self.ps['Nx']
        dx = self.ps['dx']
        wl = self.ps['wl']
        nap = self.ps['na']
        n2 = self.ps['n2']
        dp = 1/(nx*dx)
        radius = (2*nap/wl)/2/dp
        msk = U.discArr((nx, nx), radius)
        self.wfc = []
        ph0 = np.angle(mny.Cshift(self.bps[mcenter]))
        for bp in self.bps:
            ph = np.angle(mny.Cshift(bp))
            zarr = mny.VecZernDecomp(ph-ph0,self.ps,phase=True)
            if remttf:
                zarr[:4]=0.0
            phi = mny.BuildBPP(zarr,self.ps)
            self.wfc.append(phi)
        return True

    def getZarr(self):
        # get mask
        nx = self.ps['Nx']
        dx = self.ps['dx']
        wl = self.ps['wl']
        nap = self.ps['na']
        n2 = self.ps['n2']
        dp = 1/(nx*dx)
        radius = (2*nap/wl)/2/dp
        msk = U.discArr((nx, nx), radius)
        # ###############
        bpp = self.bps[0]
        phase = np.angle(mny.Cshift(bpp))*msk
        zarr = mny.VecZernDecomp(phase,self.ps,phase=True)
        return zarr
        
        
    def showphaseandamp(self,no):
        # get mask
#        nx = self.ps['Nx']
#        dx = self.ps['dx']
#        wl = self.ps['wl']
#        nap = self.ps['na']
#        n2 = self.ps['n2']
#        dp = 1/(nx*dx)
#        radius = (2*nap/wl)/2/dp
#        msk = F.discArr((nx, nx), radius)
#        # ###############
#        bpp = self.bps[no]
#        phase = N.angle(mny.Cshift(bpp))*msk
#        amp = abs(mny.Cshift(bpp))
        amp = self.amp[no]
        phase = self.phase[no]
        tf.imshow(phase,vmin=None)
        tf.imshow(amp,vmin=None)
        return True

class mre(object):
    
    def __init__(self, stack):
        self.stack = stack
        nz, nx, ny = stack.shape
        self.Nx = nx # image size, in pixels
        self.nz = nz # number of slices
        self.dx = 0.089 # pixel size, in microns
        self.dz = 0.20 # distance in microns between slices
        self.n2 = 1.4 # refractive index at sample
        self.na = 1.3 # objective numerical aperture
        self.wl = 0.515 # wavelength, microns
        self.zcenter = 0.0 # which z slice is the focal plane (microns)
        self.offset = 200 # the mean background level (cannot be zero for unsigned int data)
        # further derived parameters
        self.dp = 1.0 / ( self.Nx * self.dx )
        self.radius = (2*self.na/self.wl)/2/self.dp
        # initialize
        self.gamma = self.getgamma() # gamma is the defocus parameter (see updateU and updateA)
    
    def __del__(self):
        pass
    
    def set_params(self,dx=None,na=None,wl=None):
        if not dx==None:
            self.dx = dx
        if not na==None:
            self.na = na
        if not wl==None:
            self.wl = wl
        self.dp = 1.0/(self.Nx*self.dx)
        self.radius = (2*self.na/self.wl)/2/self.dp
        self.gamma = self.getgamma() 
        
    def run(self,ncyc,bpp=None):
        ''' run phase retrieval '''
        nx = self.Nx
        ny = nx
        nz = self.nz
        self.A = np.zeros((nx,ny),dtype=np.complex128) # wavefront
        self.U = np.zeros((nz,nx,ny),dtype=np.complex128) # amplitude of the light
        self.I = np.zeros((nz,nx,ny),dtype=np.float64) # Intensity |u|^2+offset
        if bpp==None: # the process is normally started with a set starting guess
            self.InitializeA(1.0)
        else: # you can provide your own starting guess
            self.A[:] = bpp
        xl = []
        err = []
        for m in range(ncyc): # Here we start the iteration
            self.updateU() # U(x,y,z) is the FT of A with the proper amount of defocus
            self.updateI() # |u|^2+offset
            self.updateA() # go back to A, inverse FT of U with defocus accounted for
            xl.append(self.getL()) # at the moment, I forget the significance of this metric
            err.append(self.geterr()) # difference between I and stack (would be zero ideally with no noise)
            print('%d / %d : %f, %f' % (m,ncyc,xl[m],err[m]))
##            if (m%8==0):
##                print "averaging"
##                A[:] = bppavg(A)
        print('Strehl ratio: %f' % self.getStrehl())
        return True
        
    def InitializeA(self,t):
        ''' set A, the wavefront, to a flat disk of the proper numerical aperture '''
        Nx = self.Nx
        radius = self.radius
        msk = U.discArray((Nx,Nx),radius,origin=None,dtype=np.float64)
        self.A[:,:] = t*fftshift(msk)
        return True
        
    def updateU(self):
        ''' update the field around focus from the wavefront A '''
        dz = self.dz
        nz = self.nz
        zcenter = self.zcenter
        for m in range(nz): # calculate for each z slice
            z = dz*m-zcenter
            self.U[m] = ifft2(self.A*np.exp(1j*z*self.gamma)) # FFT
        return True
        
    def updateI(self):
        mu = self.offset
        self.I[:] = np.abs(self.U)**2+mu
        return True
    
    def updateA(self):
        ''' This routine adjusts U to agree with the measured data in stack '''
        nx = self.Nx
        #dx = self.dx
        dz = self.dz
        nz = self.nz
        #n2 = self.n2
        #nap = self.na
        #wl = self.wl
        zcenter = self.zcenter
        radius = self.radius
        msk = U.discArray((nx,nx),radius,origin=(0,0),dtype=np.float64)
        self.A[:,:] = 0.0
        for m in range(nz):
            z = dz*m-zcenter
            kernel = ((self.stack[m]/self.I[m])*self.U[m]) # here's the key step
            pf = np.exp(-1j*z*self.gamma) # defocus factor
            self.A[:,:] = msk*(self.A[:,:] + (1./nz)*pf*fft2(kernel)) #FFT
        #A[:,:] = abs(A.max())*msk*N.exp(1j*N.angle(A)).astype(N.complex64)
        return True    
    
    def getgamma(self):
         ''' gamma is the defocus parameter '''
         Nx = self.Nx
         radius = self.radius
         msk = U.discArray((Nx,Nx),radius,origin=None,dtype=np.float64)
         sinphim = (self.na/self.n2)
         w = np.zeros((Nx,Nx), np.float64)
         rho = np.fromfunction(lambda i,j: (1./radius)*np.sqrt((i-Nx/2)**2+(j-Nx/2)**2), (Nx,Nx))
         rho = msk*rho
         w = msk*2*pi*self.n2*np.sqrt(1-(sinphim*rho)**2)/self.wl
         return fftshift(w)
         
    def getL(self):
        t = 0.0
        for m in range(self.nz):
            t = t + (self.stack[m]*np.log(self.I[m])).sum()
        return t
        
    def geterr(self):
        norm = self.Nx*self.Nx*self.nz
        t = ((self.I-self.stack)**2).sum()/norm
        return t
        
    def getStrehl(self):
        ''' I think this is more correct '''
        Nx = self.Nx
        radius = self.radius
        msk = U.discArray((Nx,Nx),radius,origin=None,dtype=np.float32)
        Np = msk.sum()
        num = abs(self.A.sum())**2
        den = Np*(abs(self.A)**2).sum()
        return (num/den)
        
    def get_amp_phase(self, wrapped=False):
        nx = self.Nx
        radius = self.radius
        msk = U.discArray((nx,nx),radius,origin=None,dtype=np.float64)
        amp = np.abs(fftshift(self.A))
        phase = msk*np.angle(fftshift(self.A))
        if not wrapped:
            zarr = vz.VecZernDecomp(phase,nx,radius)
            phase = vz.buildphiZ(zarr,amp.shape,rad=radius)
        return (amp,phase)
    
def test():
    half_size = 64
    img = tf.imread('bd_stack_3_948nm_20090413.tif')
    # get image center
    n3 = img.argmax()
    nz, nx, ny = img.shape
    zc = n3/(nx*ny)
    n2 = n3 % (nx*ny)
    xc = n2 / ny
    yc = n2 % ny    
    # cutout image and shift center to (0,0)
    imgout = fftshift(img[:,(xc-half_size):(xc+half_size),(yc-half_size):(yc+half_size)],(1,2))
    # start phase retrieval
    p = mre(imgout)
    p.set_params(dx=0.0948) # set pixel size, microns
    p.dz = 0.2 # set spacing between z slices, microns
    p.zcenter = p.dz*zc
    p.offset = 235
    p.run(32)
    amp, phase = p.get_amp_phase()
    pylab.figure()
    pylab.imshow(amp, interpolation='nearest')
    pylab.figure()
    pylab.imshow(phase, interpolation='nearest')
    return imgout.min()