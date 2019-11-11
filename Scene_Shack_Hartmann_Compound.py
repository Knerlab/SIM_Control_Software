# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:04:19 2016

@author: keelan
"""
import scipy.ndimage.measurements
com = scipy.ndimage.measurements.center_of_mass
#import AdaptiveThreshold_Opencv as opencv
import tifffile as tf
import numpy as N
from pylab import imshow
#from scipy.ndimage import correlate as corr
from scipy.signal import fftconvolve as corr
from scipy.signal import convolve as convolve
from scipy.ndimage.measurements import center_of_mass as sccom
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import detrend

from PyQt4 import QtGui
from PyQt4 import QtCore

fft2 = N.fft.fft2
ifft2 = N.fft.ifft2
fftshift = N.fft.fftshift
pi = N.pi

class Wavefront_Sensor(QtCore.QObject):
    
    def __init__(self):
        # set up inital parameters to determine size of the scene composition
        #self.adapt = opencv.adaptivethreshold()
        self.radius = 17 # 1/2 the total number of lenslets in linear direction
        self.x_center = 1022
        self.y_center = 586
        self.px_spacing = 32.85 # spacing between each lenslet
        self.barrier = 256.0 # maximum intensity value accepted by the camera
        self.threshold = 100.0 # minimum intensity value accepted by the camera
        self.nx = 2*self.radius # number of lenslets x-direction
        self.ny = 2*self.radius # number of lenslets y-direction
        self.hsp = 8 # size of subimage is 2*hsp
        self.calfactor = (.00586/6.7)*(150) # pixel size * focalLength * pitch
        self.NaN_Mask = self.mask_nan()
        self.Zero_Mask = self.mask_zero()
        # set up seccorr center
        section = N.ones((2*self.hsp,2*self.hsp))
        sectioncorr = corr(1.0*section, 1.0*section[::-1,::-1], mode='full')
        self.CorrCenter = N.unravel_index(sectioncorr.argmax(), sectioncorr.shape)
        # Set up T and W filters for waffle removal
        # T
        self.T_Filter = N.zeros((3,3))
        self.T_Filter[0,1] = self.T_Filter[1,0] = self.T_Filter[1,2] = self.T_Filter[2,1] = 0.25
        self.T_Filter[1,1] = 1
        self.T_Filter = 0.5*self.T_Filter
        # W
        self.W_Filter = N.ones((3,3))
        self.W_Filter[0,0] = self.W_Filter[0,2] = self.W_Filter[2,0] = self.W_Filter[2,2] = 0.25
        self.W_Filter[0,1] = self.W_Filter[1,0] = self.W_Filter[1,2] = self.W_Filter[2,1] = 0.5
        self.W_Filter = 0.25*self.W_Filter
        self.Test_Waffle = N.zeros((3,3))
        self.Test_Waffle2 = N.zeros((3,3))
        self.Test_Waffle[0,1] = self.Test_Waffle[1,0] = self.Test_Waffle[1,2] = self.Test_Waffle[2,1] = 1.0
        self.Test_Waffle2[0,0] = self.Test_Waffle2[0,2] = self.Test_Waffle2[1,1] = self.Test_Waffle2[2,0] = self.Test_Waffle2[2,2] = 1.0
        # initialize Arrays
        self.gradx = N.zeros((self.ny,self.nx))
        self.grady = N.zeros((self.ny,self.nx))
        self.gradx_guide = N.zeros((self.ny,self.nx))
        self.grady_guide = N.zeros((self.ny,self.nx))
        self.gradx_Ref = N.zeros((self.ny,self.nx))
        self.grady_Ref = N.zeros((self.ny,self.nx))
        self.hudgins_prep()
        
    def GetAberration2img(self,baseimg,offsetimg):
        baseimg = self.high_pass_filter(baseimg)
        offsetimg = self.high_pass_filter(offsetimg)
        baseimg = self.low_pass_filter(baseimg)
        offsetimg = self.low_pass_filter(offsetimg)
        gradx,grady = self.GetGradientsCorr(baseimg,offsetimg)
        #gradx = gradx*self.Zero_Mask
        #grady = grady*self.Zero_Mask
        extx,exty = self.hudgins_extend_mask0(gradx,grady)
        phi = self.recon_hudgins(extx,exty)
        phi = phi*self.calfactor
        phi = self.RemoveGlobalWaffle(phi)
        #phi = phi*self.Zero_Mask
        #phi = self.RemoveLocalWaffle(phi)
        return phi
        
    def GetAberration1img(self,image):
        image = self.high_pass_filter(image)
        image = self.low_pass_filter(image)
        gradx,grady = self.GetGradientsRef(image)
        #gradx = self.mask*gradx
        #grady = self.mask*grady
        median_x = self.Median_reduction(gradx)
        median_y = self.Median_reduction(grady)
        gradx = 1.0*gradx - median_x
        grady = 1.0*grady - median_y
        gradx = gradx*self.Zero_Mask
        grady = grady*self.Zero_Mask
        extx,exty = self.hudgins_extend_mask0(gradx,grady)
        phi = self.recon_hudgins(extx,exty)
        phi = phi*self.calfactor
        phi = self.RemoveGlobalWaffle(phi)
        #phi = phi*self.NaN_Mask
        phi = phi*self.Zero_Mask
        return phi
        
    def GetGradientsCorr(self,baseimg,offsetimg):
        ''' Determines Gradients by Correlating each section
            with its base reference section'''
        for ii in range(self.nx):
            for jj in range(self.ny):
                gradx,grady = self.findDotsCorrelateoff(baseimg,offsetimg,jj,ii)
                self.gradx[jj,ii] = gradx# + 0.5
                self.grady[jj,ii] = grady# + 0.5       
        return self.gradx,self.grady
        
    def GetGradientsGuide(self,baseimg,offsetimg):
        ''' Determines Gradients using center of mass'''
        for ii in range(self.nx):
            for jj in range(self.ny):
                gradx,grady = self.findDotsCOM(baseimg,offsetimg,jj,ii)
                self.gradx_guide[jj,ii] = gradx# + 0.5
                self.grady_guide[jj,ii] = grady# + 0.5       
        return self.gradx_guide,self.grady_guide
        
    def GetGradientsRef(self,image):
        ''' Determines Gradients from a single image by Correlating
            Each section with a reference section'''
        for ii in range(self.nx):
            for jj in range(self.ny):
                gradx,grady = self.findDotsRef(image,jj,ii)
                self.gradx_Ref[jj,ii] = gradx
                self.grady_Ref[jj,ii] = grady
        #self.gradx_Ref = self.gradx_Ref*self.mask_grads
        #self.grady_Ref = self.grady_Ref*self.mask_grads
        #mask[mask == 0.0] = 'nan'
        #self.gradx_Ref[self.gradx_Ref == 0.0] = 'nan'
        #self.grady_Ref[self.grady_Ref == 0.0] = 'nan'
        #self.gradx_Ref = (self.Median_reduction(self.gradx_Ref))#*(-1.0)# + 0.5
        #self.grady_Ref = (self.Median_reduction(self.grady_Ref))#*(-1.0)# + 0.5
        #self.gradx_Ref[N.isnan(self.gradx_Ref)] = 0.0
        #self.grady_Ref[N.isnan(self.grady_Ref)] = 0.0 
        #self.gradx_Ref[self.gradx_Ref == 'nan'] = 0.0
        #self.grady_Ref[self.grady_Ref == 'nan'] = 0.0
        return self.gradx_Ref,self.grady_Ref
        
    def Median_reduction(self,gradient):
        # Reduces the values of all gradients by the values of the median gradient
        gradient = gradient*self.NaN_Mask
        median = N.nanmedian(gradient)
        return median
        
    def setZernike(self,mode, amp):
        phiin = amp*zernike.Zm(mode,rad=11, orig=None,Nx=256) # was 28
        # set mirror with new dm shape
        return phiin
        
    def ZernikeLoop(self,amp=0.5,radius=11):
        phi_array = N.zeros((25,2*radius,2*radius))
        for k in range(25):
            phi = self.setZernike(k,amp)
            phi_array[k] = phi
        return phi_array
        
#    def LinearReg_reduction(self,gradient):
        # Redcues the values off all gradients by the central value of the linear fit
        
        
    def findDot(self,image,iy,ix): 
        hsp = self.hsp
        bot = self.y_center - (self.radius)*self.px_spacing
        left = self.x_center - (self.radius)*self.px_spacing
        vert = int(bot + iy*self.px_spacing)
        horiz = int(left + ix*self.px_spacing)
        sec = image[(vert-hsp):(vert+hsp),(horiz-hsp):(horiz+hsp)]
        return sec
        
    def findDots(self,image):
        dotarray = N.zeros((self.nx*self.ny,2*self.hsp,2*self.hsp))
        count = 0
        for ii in range(self.nx):
            for jj in range(self.ny):
                sec = self.findDot(image,jj,ii)
                dotarray[count] = sec
                #print jj,ii
                count = count + 1
        return dotarray
        
    def findDotsCorrelateoff(self,baseimg,offsetimg,iy,ix):
        '''finds new spots using each section correlated with the center'''
        hsp = self.hsp
        bot = self.y_center - (self.radius)*self.px_spacing
        left = self.x_center - (self.radius)*self.px_spacing
        vert = int(bot + iy*self.px_spacing)
        horiz = int(left + ix*self.px_spacing)
        sec = offsetimg[(vert-hsp):(vert+hsp),(horiz-hsp):(horiz+hsp)]
        secbase = baseimg[(vert-hsp):(vert+hsp),(horiz-hsp):(horiz+hsp)]
        seccorr = corr(1.0*secbase, 1.0*sec[::-1,::-1], mode='full')
        #seccorrnorm = seccorr*(1.0/N.sum(seccorr))
        try:
        # Corrects for sub-images leaving the sub-field of view by setting the local gradient to zero
            px,py = self.Parabolicfit(seccorr)
        except: #IndexError
            px,py = self.CorrCenter[1],self.CorrCenter[0]
            #print ix,iy, "Scene Left SubImage"
        gradx = self.CorrCenter[1] - px
        grady = self.CorrCenter[0] - py
        return gradx,grady
        
    def findDotsCOM(self,baseimg,offsetimg,iy,ix):
        ''' Find Dot location using center of Mass'''
        hsp = self.hsp
        bot = self.y_center - (self.radius)*self.px_spacing
        left = self.x_center - (self.radius)*self.px_spacing
        vert = int(bot + iy*self.px_spacing)
        horiz = int(left + ix*self.px_spacing)
        sec = offsetimg[(vert-hsp):(vert+hsp),(horiz-hsp):(horiz+hsp)]
        secbase = baseimg[(vert-hsp):(vert+hsp),(horiz-hsp):(horiz+hsp)]
        py,px = sccom(sec)
        sy,sx = sccom(secbase)
        gradx = px - sx
        grady = py - sx
        return (gradx,grady)
        
    def findDotsRef(self,image,iy,ix):
        ''' find dot location using correlations and parabolic fit'''
        #image = self.high_pass_filter(image)
        hsp = self.hsp
        bot = self.y_center - (self.radius)*self.px_spacing
        left = self.x_center - (self.radius)*self.px_spacing
        vert = int(bot + iy*self.px_spacing)
        horiz = int(left + ix*self.px_spacing)
        vertmiddle = int(self.y_center)
        horizmiddle = int(self.x_center)
        sec = image[(vert-hsp):(vert+hsp),(horiz-hsp):(horiz+hsp)]
        #sec = self.Threshold_local(self,sec)
        #sec = self.high_pass_filter(sec)
        secmiddle = image[(vertmiddle-hsp):(vertmiddle+hsp),(horizmiddle-hsp):(horizmiddle+hsp)]
        #secmiddle = self.Threshold_local(self,secmiddle)
        #secmiddle = self.high_pass_filter(secmiddle)
        seccorr = corr(1.0*secmiddle, 1.0*sec[::-1,::-1], mode='full')
        px,py = self.Parabolicfit(seccorr)
        gradx = px
        grady = py
        return gradx,grady
        
    def FindDotsMaxInt(self,image):
        image = self.high_pass_filter(image)
        intensityarray = N.zeros((self.nx,self.nx))
        count = 0
        for ii in range(self.nx):
            for jj in range(self.ny):
                sec = self.findDot(image,jj,ii)
                maximum = sec.max()
                intensityarray[jj,ii] = maximum
                #print jj,ii
                count = count + 1
        return intensityarray
        
    #Removes Waffle mode see "Fast wave-front reconstruction in large adaptive optics systems....Poyneer,Gavel,Brase
#    def get_mask(self):
#        self.wv_mask = N.zeros((2*self.radius,2*self.radius))
#        for i in range(2*self.radius):
#            for k in range(2*self.radius):
#                if self.placement_mask[i,k] >=0:
#                    self.wv_mask[i,k] = 1
    
    def RemoveGlobalWaffle(self,phi):
        wmode = N.zeros((2*self.radius,2*self.radius)) # was originally (50,50)
        constant_num = 0
        constant_den = 0
        # a waffle-mode vector of +-1 for a given pixel of the Wavefront
        for x in range(2*self.radius):
            for y in range(2*self.radius):
                if (x+y)/2 - N.round((x+y)/2)==0:
                    wmode[y,x] = 1
                else:
                    wmode[y,x] = -1
        wmode = wmode
        self.temp = wmode
        for i in range(2*self.radius):
            for k in range(2*self.radius):
                temp = phi[i,k]*wmode[i,k]
                temp2 = wmode[i,k]*wmode[i,k]
                constant_num = constant_num+temp
                constant_den = constant_den+temp2
        constant = constant_num/constant_den
        #print constant
        phicorr = (phi - constant*wmode)
        return phicorr
        
    def RemoveLocalWaffle(self,phi):
        # initialize DM_map corrected for local waffle
        phi_corrected = N.zeros((2*self.radius,2*self.radius))
        for ii in range(2*self.radius-2):
            for jj in range(2*self.radius-2):
                img3x3 = phi[((1+jj)-1):((1+jj)+2),((1+ii)-1):((1+ii)+2)]
                phi_corrected[jj+1,ii+1] = N.sum(img3x3*self.T_Filter)
        phi_corrected[:,0] = phi[:,0]
        phi_corrected[:,2*self.radius-1] = phi[:,2*self.radius-1]
        phi_corrected[0,] = phi[0,]
        phi_corrected[2*self.radius-1,] = phi[2*self.radius-1,]
        return phi_corrected   
        
    def high_pass_filter(self,data):
        y,x = data.shape
        data[data<self.threshold] = 0.0
        return data
        
    def low_pass_filter(self,data):
        y,x = data.shape
        data[data>self.barrier] = 0.0
        return data
        
    def Threshold_local(self,section,threshold=0.9):
        maximum = section.max()
        section[section<threshold*maximum] = 0.0
        return section
        
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
        
    def Get_Parabolic_WF(self):
        actsimage = N.zeros((19,19))
        centerx = 9
        centery = 9
        for ii in range(19):
            for jj in range(19):
                radius = N.sqrt((centerx-ii)**2 + ((centery-jj)**2))
                actsimage[jj,ii] = 0.2 - radius*0.01
        actsimage = actsimage - actsimage.mean()
        return actsimage
        
    def SVD_to_CMD(self,SVD):
        '''Creates a 1-d vector of size 100 for the using the center 10x10 actuators '''
        cmd = 0.5*N.ones(140)
        count = 0
        for k in range(10):
            for m in range(10):
                x = m
                y = 12*(k+1)
                cmd[(x+y)-1] = SVD[count]
                count = count + 1
        return cmd
        
    def compositeimgoff(self,baseimg,threshold=0.75):
        images = []
        for i in range(self.nx):
            row = []
            for j in range(self.nx):
                dot = self.findDot(baseimg,i,j)
                dot = self.low_pass_filter(dot)
                dot = self.Threshold_local(dot,threshold)
                row.append(dot)
            images.append(row)
            
        flat = []
        for image_row in images:
            for flatrow_index in xrange(len(images[0][0])):
                flatrow = []
                for image in image_row:
                    for pixel in image[flatrow_index]:
                        flatrow.append(pixel)
                flat.append(flatrow)
        flatarray = N.asarray(flat)
        return flatarray
        
    def mask_zero(self):
        # Mask
        a = self.radius #14  # center variable of mask (a,b)
        b = self.radius #14  # center variable of mask (a,b)
        Ly = 2*self.radius #28 # y dimension of image
        Lx = 2*self.radius #28 # x dimension of image
        r1 = self.radius -1 #13 # radius of mask
        y, x = N.ogrid[-a:Ly-a, -b:Lx-b] #creates a meshgrid centered at (a,b)
        mask = x**2 + y**2 <= r1**2
        mask[1,a] = mask[a,1] = mask[a,Lx-1] = mask[Ly-1,a] = 0.0
        mask = mask.astype('float')
        #mask[mask == 0.0] = 'nan'
        return mask
    
    def mask_nan(self):
        a = self.radius  # center variable of mask (a,b)
        b = self.radius  # center variable of mask (a,b)
        Ly = 2*self.radius # y dimension of image
        Lx = 2*self.radius # x dimension of image
        r1 = self.radius -1 # radius of mask
        y, x = N.ogrid[-a:Ly-a, -b:Lx-b] #creates a meshgrid centered at (a,b)
        mask = x**2 + y**2 <= r1**2
        mask[1,a] = mask[a,1] = mask[a,Lx-1] = mask[Ly-1,a] = 0.0
        mask = mask.astype('float')
        mask[mask == 0.0] = 'nan'
        return mask
        
    def recon_hudgins(self,gradx,grady):
        ''' wavefront reconstruction from gradients
            Hudgins Geometry, Poyneer 2002 '''
        sx = fft2(fftshift(gradx))
        sy = fft2(fftshift(grady))
        numx = self.numx
        numy = self.numy
        den = self.den
        den[0,0] = 1.e-6
        self.sw = (numx*sx + numy*sy)/den
        self.sw[0,0] = 0.0
        #self.sw[nx/2,nx/2] = 0.0
        phi = fftshift(ifft2(self.sw)).real
        return phi

    def hudgins_extend_masknan(self,gradx,grady):
        ''' extension tecchnique Poyneer 2002 '''
        nx = self.nx
        ny = self.ny
        if (nx % 2 == 0): #even
            mx = nx/2
        else: #odd
            mx = (nx+1)/2        
        for jj in range(int(nx)):
            for ii in range(int(mx),int(nx)):
                if grady[jj,ii] == ('nan'):
                    grady[jj,ii] = grady[jj,ii-1]
                if gradx[ii,jj] == ('nan'):
                    gradx[ii,jj] = gradx[ii-1,jj]
            for ii in range(int(mx),int(-1),int(-1)):
                if grady[jj,ii] == ('nan'):
                    grady[jj,ii] = grady[jj,ii+1]
                if gradx[ii,jj] == ('nan'):
                    gradx[ii,jj] = gradx[ii+1,jj]
        gradx[:,nx-1] = -gradx[:,:(nx-1)].sum(1)
        grady[nx-1,:] = -grady[:(nx-1),:].sum(0)
        return gradx,grady
    
    def hudgins_extend_mask0(self,gradx,grady):
        ''' extension tecchnique Poyneer 2002 '''
        nx = self.nx
        ny = self.ny
        if (nx % 2 == 0): #even
            mx = nx/2
        else: #odd
            mx = (nx+1)/2        
        for jj in range(int(nx)):
            for ii in range(int(mx),int(nx)):
                if grady[jj,ii] == (0.0):
                    grady[jj,ii] = grady[jj,ii-1]
                if gradx[ii,jj] == (0.0):
                    gradx[ii,jj] = gradx[ii-1,jj]
            for ii in range(int(mx),int(-1),int(-1)):
                if grady[jj,ii] == (0.0):
                    grady[jj,ii] = grady[jj,ii+1]
                if gradx[ii,jj] == (0.0):
                    gradx[ii,jj] = gradx[ii+1,jj]
        gradx[:,nx-1] = -gradx[:,:(nx-1)].sum(1)
        grady[nx-1,:] = -grady[:(nx-1),:].sum(0)
        return gradx,grady
        
    def hudgins_prep(self):
        nx = self.nx
        ny = self.ny
        k,l = N.meshgrid(N.arange(ny),N.arange(nx))
        self.numx = (N.exp(-2j*pi*k/nx)-1)
        self.numy = (N.exp(-2j*pi*l/nx)-1)
        self.den = 4*(N.sin(pi*k/nx)**2 + N.sin(pi*l/nx)**2)
        return True
        