# -*- coding: utf-8 -*-
"""
@copywrite, Ruizhe Lin and Peter Kner, University of Georgia, 2019
"""

import os, sys
temppath = 'C:/Users/Public/Documents/python_code/DMtestSH/temp'
join = lambda fn: os.path.join(temppath,fn)

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import uic
import numpy as np
import Scene_Shack_Hartmann_Compound as shwfs
from numpy.random import randn
from matplotlib.colorbar import Colorbar as colorbar
from scipy.signal import fftconvolve as corr
import tifffile as tf
import PySpin
from asdk import DM
import time
import pr_batch_P33_0 as prb
import zernike

infl_fcn_file = r'C:/Users/Public/Documents/python_code/SLM_P27/influ_func_08062018.tif'
infl_fcn = tf.imread(infl_fcn_file)
Sm = prb.Smat(infl_fcn)

class Form(QDialog):
    
    def __init__(self):
        super(Form, self).__init__()
        uic.loadUi('Gui_FLIR.ui', self)    
        self.connect(self.QuitButton, SIGNAL("clicked()"),self.shutdown)
        self.connect(self.RunButton, SIGNAL("clicked()"), self.run)
#        self.connect(self.StopButton, SIGNAL("clicked()"), self.stop)
        self.connect(self.SaveButton, SIGNAL("clicked()"), self.Saveimg)
        self.connect(self.SaveWFButton, SIGNAL("clicked()"), self.savedmwf)
        self.connect(self.SaveZMWFButton, SIGNAL("clicked()"), self.savezmwf)
        self.connect(self.AutoScanButton, SIGNAL("clicked()"), self.autoscandm)
        self.connect(self.exposurebox, SIGNAL("returnPressed()"), self.exposureupdate)
        self.connect(self.thresholdbox, SIGNAL("returnPressed()"), self.thresholdupdate)
        self.connect(self.PushDMButton, SIGNAL("clicked()"), self.pushdm)
        self.connect(self.ScanZernikeModeButton, SIGNAL("clicked()"), self.scanzernikemode)
        self.connect(self.LoadDMButton, SIGNAL("clicked()"), self.loaddmfile)
        self.connect(self.InitiatePushButton, SIGNAL("clicked()"), self.set_baseimage)
        self.connect(self.Infl_func_pushButton, SIGNAL("clicked()"), self.influencefunction)
        self.ProgressBar.setRange(0,100)
        self.ProgressBar.setValue(0)
        self.index = 0
        #DM
        self.serialName = 'BAX228'
        self.dm = DM( self.serialName )
        self.nbAct = int( self.dm.Get('NBOfActuator') )
        self.values = [0.] * self.nbAct
        self.dm.Send( self.values )
        #Initial camera
        self.system = PySpin.System.GetInstance()
        self.cam_list = self.system.GetCameras()
        self.cam = self.cam_list[0]
        self.cam.Init()
        self.nodemap = self.cam.GetNodeMap()
        self.node_acquisition_mode = PySpin.CEnumerationPtr(self.nodemap.GetNode('AcquisitionMode'))
        self.node_acquisition_mode_continuous = self.node_acquisition_mode.GetEntryByName('Continuous')
        self.acquisition_mode_continuous = self.node_acquisition_mode_continuous.GetValue()
        self.node_acquisition_mode.SetIntValue(self.acquisition_mode_continuous)
        self.cam.ExposureAuto.SetValue(PySpin.ExposureAuto_Off)
        self.exposure_time_to_set = 200.0
        self.cam.ExposureTime.SetValue(self.exposure_time_to_set)
        self.cam.BeginAcquisition()
        self.image_result = self.cam.GetNextImage()
        self.img = self.image_result.GetNDArray()
        self.image = self.img.copy()
        self.baseimage = self.img.copy()
        self.image_result.Release()
        self.cam.EndAcquisition()
        self.ah = self.CCDImage.axes.imshow(self.image)
        # setup the wavefront
        self.WFS = shwfs.Wavefront_Sensor()
        self.phi = self.WFS.GetAberration2img(self.baseimage,self.baseimage)
        self.phi = self.phi*self.WFS.Zero_Mask
        self.wah = self.WFimage.axes.imshow(self.phi,interpolation='nearest',vmin=-0.8,vmax=0.8)
        self.WFimage.figure.colorbar(self.wah)
        
    def set_baseimage(self,):
        self.cam.BeginAcquisition()
        self.image_result = self.cam.GetNextImage()
        self.baseimage = self.image_result.GetNDArray().copy()
        self.image_result.Release()
        self.cam.EndAcquisition()
        
    def run(self):
        self.cam.BeginAcquisition()
        self.image_result = self.cam.GetNextImage()
        self.image = self.image_result.GetNDArray().copy() 
        self.image_result.Release()
        self.cam.EndAcquisition()
        self.imageupdate()
        
    def imageupdate(self):
        # update ccd image
        self.ah.set_data(self.image)
        self.CCDImage.figure.canvas.draw()
        # measure wavefront
        self.phi = self.WFS.GetAberration2img(self.baseimage,self.image)
        self.phi = self.phi*self.WFS.Zero_Mask
        self.wah.set_data(self.phi)
        self.wah.set_clim(self.phi.min(),self.phi.max())
        self.WFimage.figure.canvas.draw()
        self.ProgressBar.setValue(self.index)
        self.index = (self.index + 1) % 100
        
    def shutdown(self):
        x = self.dm.Reset()
        if (x==0):
            print('DM Reset')
        else:
            print('DM cannot reset')
        del self.cam
        self.cam_list.Clear()
        self.system.ReleaseInstance()
        self.close()

    def pushdm(self):
        push = self.PushdoubleSpinBox.value()
        num = self.ActspinBox.value()
        self.values[num] = push
        self.dm.Send( self.values )
        self.values[num] = 0.
        return True
        
    def influencefunction(self):
        n = self.nbAct
        push = self.PushdoubleSpinBox.value()
        for i in range(n):
            self.values[i] = push
            self.dm.Send( self.values )
            self.values[i] = 0.
            self.run()
            fn = 'act_%d'%i + '_v%f.tif'%push
            tf.imsave(fn,self.phi.astype(np.float32),photometric='minisblack')       
        return True
        
    def autoscandm(self):
        for i in range( self.nbAct ):
            data = []
            for p in range (19):
                push = -0.9+p*0.1
                self.values[i] = push
                self.dm.Send( self.values )
                self.cam.BeginAcquisition()
                self.image_result = self.cam.GetNextImage()
                self.image = self.image_result.GetNDArray().copy() 
                self.image_result.Release()
                self.cam.EndAcquisition()
                self.phi = self.WFS.GetAberration2img(self.baseimage,self.image)
                self.phi = self.phi*self.WFS.Zero_Mask
                data[p] = self.phi
                self.values[i] = 0.   
                time.sleep(0.1)
            fn = 'act_%d.tif'%(i)
            tf.imsave(join(fn),data.astype(np.float32),photometric='minisblack')
        return True
        
    def scanzernikemode(self):
        n = self.ZernikeModeBox.value()
        h = self.PushSpinBox.value()
        phiin = h*zernike.Zm(n,rad=17, orig=None,Nx=33)
        dmarr = 0.1*np.dot(Sm.S,phiin.reshape(33*33))
        dmfile = [0.] * 69
        for i in range(69):
            dmfile[i] = dmarr[i]
        if (all(i <= 1.0 for i in dmfile)):
            self.dm.Send(dmfile)
        else:
            raise Exception(' Error: push value greater than 1.0 ')
            
    def loaddmfile(self):
        fns = self.dmfilebox.text()
        fn = r'%s' %fns
        with open(fn, "r") as file:
            cmd = eval(file.readline())
        if (all(i <= 1.0 for i in cmd)):
            self.dm.Send(cmd)
        else:
            raise Exception(' Error: push value greater than 1.0 ')
        
    def savedmwf(self):
        m = self.ActspinBox.value()
        v = self.PushdoubleSpinBox.value()
        fn = 'act_%d'%m + '_v%f.tif'%v
        tf.imsave(fn,self.phi.astype(np.float32),photometric='minisblack')
    
    def savezmwf(self):
        m = self.ZernikeModeBox.value()
        v = self.PushSpinBox.value()
        fn = 'zernike_mode%d'%m + '_v%f.tif'%v
        tf.imsave(fn,self.phi.astype(np.float32),photometric='minisblack')
        
    def Saveimg(self):
        tf.imsave('camera_image.tif',self.img.astype(np.float32),photometric='minisblack')
        tf.imsave('wfs_image.tif',self.phi.astype(np.float32),photometric='minisblack')
   
    def updatespacing(self):
        text = self.Spacinglinebox.text()
        self.px_spacing = float(text)
        
    def updateXcenter(self):
        xcenter = self.xCenterbox.text()
        self.x_center = float(xcenter)

    def updateYcenter(self):
        ycenter = self.yCenterbox.text()
        self.y_center = float(ycenter)
    
    def exposureupdate(self):
        exposure = self.exposurebox.text()
        exposure = int(exposure)
        self.cam.ExposureTime.SetValue(exposure)
    
    def thresholdupdate(self):
        threshold = self.thresholdbox.text()
        self.theshold = float(threshold)
        
    def high_pass_filter(self,data):
        y,x = data.shape
        data[data<self.threshold] = 0.0
        return data
    
if __name__=="__main__":
    app = QApplication(sys.argv)
    form = Form()
    form.show()
    app.exec_()
   