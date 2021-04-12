# -*- coding: utf-8 -*-
"""
An example of running the 3D structured illumination microscopy image reconstruction codes
@copywrite, Ruizhe Lin and Peter Kner, University of Georgia, 2019
"""

import numpy as N
from scipy.signal import fftconvolve as corr
from skimage.filters import threshold_otsu
from scipy.ndimage import center_of_mass as com
import circle as c
fft2 = N.fft.fft2
ifft2 = N.fft.ifft2
fftshift = N.fft.fftshift
pi = N.pi

class Wavefront_Reconstruction():
    def __init__(self):
        # set up inital parameters to determine size of the scene composition
#        self.radius = 8 # 1/2 the total number of lenslets in linear direction
        self.diameter = 16
        self.x_center_base = 1186
        self.y_center_base = 903
        self.x_center_offset = 1186
        self.y_center_offset = 903
        self.px_spacing = 26 # spacing between each lenslet
        self.hsp = 12 # size of subimage is 2*hsp
        self.calfactor = (.0065/4.1)*(150) # pixel size * focalLength * pitch
        # set up seccorr center
        section = N.ones((2*self.hsp,2*self.hsp))
        sectioncorr = corr(1.0*section, 1.0*section[::-1,::-1], mode='full')
        self.CorrCenter = N.unravel_index(sectioncorr.argmax(), sectioncorr.shape)
                
    def GetAberration2img(self,baseimg,offsetimg,method=1):
        # initialize Arrays
        self.nx = self.diameter
        self.ny = self.diameter
        self.im = N.zeros((2, 2*self.hsp*self.diameter, 2*self.hsp*self.diameter))
        self.gradxy = N.zeros((2, self.diameter, self.diameter))
        self.gradx = N.zeros((self.ny,self.nx))
        self.grady = N.zeros((self.ny,self.nx))
        # self.findcenter(baseimg, offsetimg)
#        self.hudgins_prep()
        if (method==0):
            print('Correlation')
            # baseimg = self.suback(baseimg)
            # offsetimg = self.suback(offsetimg)
            gradx,grady = self.GetGradientsCorr(baseimg,offsetimg)
            self.gradxy[0] = self.gradx
            self.gradxy[1] = self.grady
        if (method==1):
            print('CenterOfMass')
            baseimg = baseimg / baseimg.max()
            offsetimg = offsetimg / offsetimg.max()
            # baseimg = self.suback(baseimg)
            # offsetimg = self.suback(offsetimg)
            self.gradx, self.grady = self.GetGradientsCom(baseimg,offsetimg)
            self.gradxy[0] = self.gradx
            self.gradxy[1] = self.grady
        self.extx, self.exty = self.hudgins_extend_mask0(self.gradx, self.grady)
        self.phi = self.recon_hudgins(self.extx, self.exty)
        self.phi = self.phi*self.calfactor
        self.phicorr = self.RemoveGlobalWaffle(self.phi)
        self.msk = c.circle(self.diameter/2.,self.diameter)
        self.phicorr = self.phicorr*self.msk
        return self.phicorr
    
    def findcenter(self,base, offset):
        secbs = base[self.y_center_base-self.hsp : self.y_center_base+self.hsp, self.x_center_base-self.hsp : self.x_center_base+self.hsp]
        secof = offset[self.y_center_offset-self.hsp : self.y_center_offset+self.hsp, self.x_center_offset-self.hsp : self.x_center_offset+self.hsp]
        ind_bas = com(secbs)
        ind_off = com(secof)
        self.x_center_base = self.x_center_base-self.hsp+round(ind_bas[1])
        self.y_center_base = self.y_center_base-self.hsp+round(ind_bas[0])
        self.x_center_offset = self.x_center_offset-self.hsp+round(ind_off[1])
        self.y_center_offset = self.y_center_offset-self.hsp+round(ind_off[0])

    def GetGradientsCorr(self,baseimg,offsetimg):
        ''' Determines Gradients by Correlating each section with its base reference section'''
        for ii in range(self.nx):
            for jj in range(self.ny):
                gradx,grady = self.findDotsCorrelateoff(baseimg,offsetimg,jj,ii)
                self.gradx[jj,ii] = gradx# + 0.5
                self.grady[jj,ii] = grady# + 0.5    
        return self.gradx,self.grady
    
    def GetGradientsCom(self,baseimg,offsetimg):
        ''' Determines Gradients by Correlating each section with its base reference section'''
        for ii in range(self.nx):
            for jj in range(self.ny):
                gradx,grady = self.findDotsCOM(baseimg,offsetimg,jj,ii)
                self.gradx[jj,ii] = gradx# + 0.5
                self.grady[jj,ii] = grady# + 0.5    
        return self.gradx,self.grady

    def findDotsCorrelateoff(self,baseimg,offsetimg,iy,ix):
        '''finds new spots using each section correlated with the center'''
        hsp = self.hsp
        r = int(N.floor(self.diameter/2.))
        bot_base = self.y_center_base - r*self.px_spacing
        left_base = self.x_center_base - r*self.px_spacing
        vert_base = int(bot_base + iy*self.px_spacing)
        horiz_base = int(left_base + ix*self.px_spacing)
        bot_offset = self.y_center_offset - r*self.px_spacing
        left_offset = self.x_center_offset - r*self.px_spacing
        vert_offset = int(bot_offset + iy*self.px_spacing)
        horiz_offset = int(left_offset + ix*self.px_spacing)
        sec = offsetimg[(vert_offset-hsp):(vert_offset+hsp),(horiz_offset-hsp):(horiz_offset+hsp)]
        secbase = baseimg[(vert_base-hsp):(vert_base+hsp),(horiz_base-hsp):(horiz_base+hsp)]
        # indsec = N.where(sec == N.amax(sec))
        # indsecbase = N.where(sec == N.amax(secbase))        
        sec = self.suback(sec)
        secbase = self.suback(secbase)
        seccorr = corr(1.0*secbase, 1.0*sec[::-1,::-1], mode='full')
#        secbasecorr = corr(1.0*secbase, 1.0*secbase[::-1,::-1], mode='full')
        px,py = self.Parabolicfit(seccorr)
#        self.CorrCenter = N.unravel_index(secbasecorr.argmax(), secbasecorr.shape)
        gradx = self.CorrCenter[1] - px
        grady = self.CorrCenter[0] - py
        self.im[0, iy*2*self.hsp : iy*2*self.hsp+2*self.hsp, ix*2*self.hsp : ix*2*self.hsp+2*self.hsp] = secbase
        self.im[1, iy*2*self.hsp : iy*2*self.hsp+2*self.hsp, ix*2*self.hsp : ix*2*self.hsp+2*self.hsp] = sec
        return gradx,grady
    
    def findDotsCOM(self,baseimg,offsetimg,iy,ix):
        ''' Find Dot location using center of Mass'''
        hsp = self.hsp
        r = int(N.floor(self.diameter/2.))
        bot_base = self.y_center_base - r*self.px_spacing
        left_base = self.x_center_base - r*self.px_spacing
        vert_base = int(bot_base + iy*self.px_spacing)
        horiz_base = int(left_base + ix*self.px_spacing)
        bot_offset = self.y_center_offset - r*self.px_spacing
        left_offset = self.x_center_offset - r*self.px_spacing
        vert_offset = int(bot_offset + iy*self.px_spacing)
        horiz_offset = int(left_offset + ix*self.px_spacing)
        sec = offsetimg[(vert_offset-hsp):(vert_offset+hsp),(horiz_offset-hsp):(horiz_offset+hsp)]
        secbase = baseimg[(vert_base-hsp):(vert_base+hsp),(horiz_base-hsp):(horiz_base+hsp)]
        py,px = com(sec)
        sy,sx = com(secbase)
        gradx = px - sx
        grady = py - sy
        self.im[0, iy*2*self.hsp : iy*2*self.hsp+2*self.hsp, ix*2*self.hsp : ix*2*self.hsp+2*self.hsp] = secbase
        self.im[1, iy*2*self.hsp : iy*2*self.hsp+2*self.hsp, ix*2*self.hsp : ix*2*self.hsp+2*self.hsp] = sec
        return gradx,grady

    def RemoveGlobalWaffle(self,phi):
        wmode = N.zeros((self.diameter,self.diameter))
        constant_num = 0
        constant_den = 0
        # a waffle-mode vector of +-1 for a given pixel of the Wavefront
        for x in range(self.diameter):
            for y in range(self.diameter):
                if (x+y)/2 - N.round((x+y)/2)==0:
                    wmode[y,x] = 1
                else:
                    wmode[y,x] = -1
        wmode = wmode
        self.temp = wmode
        for i in range(self.diameter):
            for k in range(self.diameter):
                temp = phi[i,k]*wmode[i,k]
                temp2 = wmode[i,k]*wmode[i,k]
                constant_num = constant_num+temp
                constant_den = constant_den+temp2
        constant = constant_num/constant_den
        #print constant
        phicorr = (phi - constant*wmode)
        return phicorr
        
    def Parabolicfit(self,sec):
        try:
            MaxIntLoc = N.unravel_index(sec.argmax(), sec.shape)
            secsmall = sec[(MaxIntLoc[0]-1):(MaxIntLoc[0]+2),(MaxIntLoc[1] -1):(MaxIntLoc[1] + 2)]
            gradx =  MaxIntLoc[1] + 0.5*(1.0*secsmall[1,0] - 1.0*secsmall[1,2])/(1.0*secsmall[1,0] + 1.0*secsmall[1,2] - 2.0*secsmall[1,1])
            grady =  MaxIntLoc[0] + 0.5*(1.0*secsmall[0,1] - 1.0*secsmall[2,1])/(1.0*secsmall[0,1] + 1.0*secsmall[2,1] - 2.0*secsmall[1,1])
        except: #IndexError
            gradx =  self.CorrCenter[1]
            grady =  self.CorrCenter[0]
        gradx = gradx
        grady = grady
        return gradx,grady
        
    def recon_hudgins(self, gradx, grady):
        """ wavefront reconstruction from gradients
            Hudgins Geometry, Poyneer 2002 """
        sx = fft2(gradx)
        sy = fft2(grady)
        nx, ny = gradx.shape
        k, l = N.meshgrid(N.arange(ny), N.arange(nx))
        numx = (N.exp(-2j * pi * k / nx) - 1)
        numy = (N.exp(-2j * pi * l / ny) - 1)
        den = 4 * (N.sin(pi * k / nx) ** 2 + N.sin(pi * l / ny) ** 2)
        sw = (numx * sx + numy * sy) / den
        sw[0, 0] = 0.0
        phi = (ifft2(sw)).real
        return phi
    
    def hudgins_extend_mask0(self, gradx, grady):
        """ extension technique Poyneer 2002 """
        nx, ny = gradx.shape
        if nx % 2 == 0:  # even
            mx = nx / 2
        else:  # odd
            mx = (nx + 1) / 2
        if ny % 2 == 0:  # even
            my = ny / 2
        else:  # odd
            my = (ny + 1) / 2
        for jj in range(int(nx)):
            for ii in range(int(my), int(ny)):
                if grady[jj, ii] == 0.0:
                    grady[jj, ii] = grady[jj, ii-1]
            for ii in range(int(my), -1, -1):
                if grady[jj, ii] == 0.0:
                    grady[jj, ii] = grady[jj, ii + 1]
        for jj in range(int(ny)):
            for ii in range(int(mx), int(nx)):
                if gradx[ii, jj] == 0.0:
                    gradx[ii, jj] = gradx[ii-1, jj]
            for ii in range(int(mx), -1, -1):
                if gradx[ii, jj] == 0.0:
                    gradx[ii, jj] = gradx[ii+1, jj]
        gradxe = gradx.copy()
        gradye = grady.copy()
        gradxe[:, ny - 1] = -1.0 * gradx[:, :(ny - 1)].sum(1)
        gradye[nx - 1, :] = -1.0 * grady[:(nx - 1), :].sum(0)
        return gradxe, gradye
    
    def suback(self,img):
        thresh  = threshold_otsu(img)
        binary = img > thresh
        imgsub = (img- thresh) * binary 
        return imgsub
        
