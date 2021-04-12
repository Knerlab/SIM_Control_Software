# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 20:53:10 2020

@author: linruizhe
"""

import sys, time, os
sys.path.append(r'C:/Users/Public/Documents/python_code/SIM_AO_p36')

import numpy as np
import tifffile as tf
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import uic
import microscope_p35 as mm
import metric_tools_p36 as mt
import ao_tools_p36 as ao
import wfzm
import ShackHartmannWavefrontReconstruction as shwfr

class Form(QDialog):

    def __init__(self):
        super(Form, self).__init__()
        uic.loadUi('Microscope_GUI_p36.ui', self)
        self.progressBar.setRange(0,100)
        self.index = 0
        self.progressBar.setValue(0)
        self.quitButton.clicked.connect(self.quit)
        self.startButton.clicked.connect(self.video)
        self.stopButton.clicked.connect(self.stopvideo)
        self.leftButton.clicked.connect(self.priorleft)
        self.rightButton.clicked.connect(self.priorright)
        self.upButton.clicked.connect(self.priorup)   
        self.downButton.clicked.connect(self.priordown)
        # Main camera display
        self.data = np.zeros((512,512,4),dtype=np.uint8)
        self.data[:,:,3] = 0
        self.scene = QGraphicsScene()
        self.view.setScene(self.scene)
        self.display_image(self.data)
        # Microscope
        self.om = mm.scope(self)
        self.om.update.connect(self.updateimage)
        self.om.update2.connect(self.updateimage2)
        self.temperatureLCD.display(self.om.ccd.GetTemperature())
        self.SetNewCoordButton.clicked.connect(self.set_coord)
        self.ResetCoordButton.clicked.connect(self.reset_coord)
        # Piezo
        self.piezoSpinBox.valueChanged.connect(self.piezomove)
        temppos = self.om.zst.getPosition()
        self.piezoSpinBox.setValue(temppos)
        # Stack Z range
        self.StartSpinBox.setValue(temppos-2.5)
        self.StopSpinBox.setValue(temppos+2.5)
        self.StepSpinBox.setValue(0.2)
        # attenuator
        self.attButton.clicked.connect(self.setattenuat)
        # Take Image Stack
        self.simzstackButton.clicked.connect(self.take_z_stack)
        self.sectioning_zstackButton.clicked.connect(self.take_sec_z_stack)
        self.stackButton.clicked.connect(self.takestack)
        self.TimeStackButton.clicked.connect(self.take_t_stack)
        self.TimeLapseButton.clicked.connect(self.taketimelapse)
        # Save
        self.saveButton.clicked.connect(self.savestack)
        # Polarizer
        self.polarizerSpinBox.valueChanged.connect(self.polarizermove)
        # Laser power
        self.PowerSpinBox_R.valueChanged.connect(self.rlaserpower)
        self.PowerSpinBox_Y.valueChanged.connect(self.ylaserpower)
        self.PowerSpinBox_UV.valueChanged.connect(self.uvlaserpower)
        # Filter Wheel
        # self.FilterButton.clicked.connect(self.setfilter)
        # Initiate DM
        if (all(i <= 1.0 for i in ao.cmd_best)):
            self.om.dm.Send(ao.cmd_best)
        else:
            raise Exception(' Error: push value greater than 1.0 ')
        # DM
        self.setnewDMButton.clicked.connect(self.setnewdm)
        self.DMPushButton.clicked.connect(self.pushdm)
        self.Infl_func_pushButton.clicked.connect(self.influencefunction)
        self.DMScanButton.clicked.connect(self.actuator_amp_scan)
        self.ZMPushButton.clicked.connect(self.setZM)
        self.LoadDMButton.clicked.connect(self.loaddmfile)
        # Adaptive Optics
        self.AOButton_SNR.clicked.connect(self.ao_optimize_snr)
        self.AOButton_MAX.clicked.connect(self.ao_optimize_max)
        self.AOButton_HPF.clicked.connect(self.ao_optimize_hpf)
        self.SaveAOButton.clicked.connect(self.write_dmfile)
        self.zeroDMButton.clicked.connect(self.cleardm)
        self.resetDMButton.clicked.connect(self.reloadflatfile)
        # SLM
        self.SelectOrderButton.clicked.connect(self.selectSLM)
        self.ActivateSLMButton.clicked.connect(self.activateSLM)
        self.DeactivateSLMButton.clicked.connect(self.deactivateSLM)
        self.AOButton_SI.clicked.connect(self.ao_optimize_si)
        # SHWFS
        self.WFS = shwfr.Wavefront_Reconstruction()
        self.base = np.zeros((2048,2048))
        self.offset = np.zeros((2048,2048))
        self.wf = np.zeros((self.WFS.diameter, self.WFS.diameter))
        # SHWFS camera display
        self.md = 0 # 0:correlation, 1:center of mass
        self.SHimg = np.zeros((2048,2048,4),dtype=np.uint8)
        self.SHimg[:,:,3] = 0
        self.scene_sh = QGraphicsScene()
        self.view_sh.setScene(self.scene_sh)
        self.display_shimage(self.SHimg)
        # WF display
        self.WFimg = np.zeros((self.WFS.diameter, self.WFS.diameter, 4), dtype=np.uint8)
        self.WFimg[:,:,3] = 0
        self.scene_wf = QGraphicsScene()
        self.view_wf.setScene(self.scene_wf)
        self.display_wfimage(self.WFimg)
        self.ReconWFPushButton.clicked.connect(self.wf_recon)
        self.ParameterSetButton.clicked.connect(self.setparameters)
        self.BaseSetButton.clicked.connect(self.load_baseimage)
        self.WFSInitiatePushButton.clicked.connect(self.set_baseimage)     
        self.WFSRunPushButton.clicked.connect(self.WFrun)
        self.WFSSavePushButton.clicked.connect(self.SaveWFimg)
        self.CalculateWFPushButton.clicked.connect(self.calculate_wf)
        self.LoadWFButton.clicked.connect(self.setnewwf)

    def display_image(self, img):
        self.scene.clear()
        w, h, z = img.shape
        self.qimg = QImage(img,w,h,QImage.Format_RGB32)
        pixMap = QPixmap()
        pixMap.convertFromImage(self.qimg)
        self.scene.addPixmap(pixMap)
        self.view.fitInView(QRectF(0, 0, w, h), Qt.KeepAspectRatio)
        self.scene.update()
        
    def display_shimage(self, img):
        self.scene_sh.clear()
        w, h, z = img.shape
        self.qimg_sh = QImage(img,w,h,QImage.Format_RGB32)
        pixMap = QPixmap()
        pixMap.convertFromImage(self.qimg_sh)
        self.scene_sh.addPixmap(pixMap)
        self.view_sh.fitInView(QRectF(0, 0, w, h), Qt.KeepAspectRatio)
        self.scene_sh.update()
        
    def display_wfimage(self, img):
        self.scene_wf.clear()
        w, h, z = img.shape
        self.qimg_wf = QImage(img,w,h,QImage.Format_RGB32)
        pixMap = QPixmap()
        pixMap.convertFromImage(self.qimg_wf)
        self.scene_wf.addPixmap(pixMap)
        self.view_wf.fitInView(QRectF(0, 0, w, h), Qt.KeepAspectRatio)
        self.scene_wf.update()
        
    def updateimage(self):
        self.progressBar.setValue(self.index)
        self.index = (self.index + 1) % 100
        temp = self.om.data
        nx, ny = temp.shape
        self.data = np.zeros((nx, ny, 4),dtype=np.uint8)
        self.data[:,:,3] = 0
        temp = temp - temp.min()
        temp = (255 * (temp / (temp.max()))).astype(np.uint8)
        self.data[:,:,0] = temp
        self.data[:,:,1] = temp
        self.data[:,:,2] = temp
        self.display_image(self.data)
        self.intensityBar.setValue(temp.max())
        self.MaxLCD.display(self.om.data.max())
        self.MinLCD.display(self.om.data.min())
        QApplication.processEvents()
        
    def updateimage2(self):
        self.progressBar.setValue(self.index)
        self.index = (self.index + 1) % 100
        QApplication.processEvents()
        
    def quit(self):
        self.close()
        
    def run(self):
        self.show()
        self.raise_()
        self.exec_()

    def savestack(self):
        text = self.fileEdit.text()
        self.om.saveTifA(slideName=text)

    def takestack(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setUvlaser = self.L405Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        setstart = self.StartSpinBox.value()
        setstop = self.StopSpinBox.value()
        setstep = self.StepSpinBox.value()
        led = self.setled()
        self.om.l = self.PatternSpinBox.value()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,UVlaser=setUvlaser,LED12=led)
        self.om.StackExt(setstart,setstop,setstep)

    def taketimelapse(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setUvlaser = self.L405Button.isChecked()
        led = self.setled()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        time = self.TimeSpinBox.value()
        pol = self.polarizerSpinBox.value()
        self.om.l = self.PatternSpinBox.value()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,UVlaser=setUvlaser,LED12=led)
        self.om.TimeLapseExt(time, pol)

    def take_z_stack(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setUvlaser = self.L405Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        setstart = self.StartSpinBox.value()
        setstop = self.StopSpinBox.value()
        setstep = self.StepSpinBox.value()
        led = self.setled()
        self.om.l = self.PatternSpinBox.value()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,UVlaser=setUvlaser,LED12=led)
        self.om.Stack_Patterns(setstart,setstop,setstep)
        
    def take_sec_z_stack(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setUvlaser = self.L405Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        setstart = self.StartSpinBox.value()
        setstop = self.StopSpinBox.value()
        setstep = self.StepSpinBox.value()
        led = self.setled()
        self.om.l = self.PatternSpinBox.value()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,UVlaser=setUvlaser,LED12=led)
        self.om.Stack_Sectioning(setstart,setstop,setstep)
        
    def take_t_stack(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setUvlaser = self.L405Button.isChecked()
        led = self.setled()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        time = self.TimeSpinBox.value()
        pol = self.polarizerSpinBox.value()
        angle = self.angle_spinBox.value()
        self.om.l = self.PatternSpinBox.value()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,UVlaser=setUvlaser,LED12=led)
        self.om.Image_Patterns(angle, time, pol)

    def param_update(self):
        if (not self.om.is_videomode):
            self.temperatureLCD.display(self.om.ccd.GetTemperature())
        self.exposureLCD.display(self.om.ccd.GetExposureTime())
        # self.exposureLCD_2.display(self.cam.)
        self.piezoLCD.display(self.om.zst.getPosition())
        if (self.om.laser647):
            self.redlaserpowerLCD.display(self.om.ll647.GetPower())
        if (self.om.laser561):
            self.greenlaserpowerLCD.display(self.om.ll561.GetPower())
        if (self.om.laser405):
            self.UVlaserpowerLCD.display(self.om.ll405.GetPower())
        self.om.prior.getPosition()

    def video(self):
        self.param_update()
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setUvlaser = self.L405Button.isChecked()
        led = self.setled()
        setFTM = self.FTMButton.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        self.om.normal_mode = self.Mode_selectButton.isChecked()
        self.om.norecordimg = self.RecordingImgButton.isChecked()
        self.om.l = self.PatternSpinBox.value()
        self.om.k = self.angle_spinBox.value()
        self.om.open_Acq(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,UVlaser=setUvlaser,
                         LED12=led,FTM=setFTM,conv=False,ccd=True,trig=1)
        self.om.start()

    def stopvideo(self):
        self.om.is_videomode = False
        time.sleep(0.01)
        self.param_update()
        
    def set_coord(self):
        xcoord = self.coordXSpinBox.value()
        ycoord = self.coordYSpinBox.value()
        nxy = self.NXYSpinBox.value()
        f = self.om.setCoord2(xcoord,ycoord,nxy)
        print(f)
        print('New Coordinate Set')
        return True
        
    def reset_coord(self):
        f = self.om.reSetCoord(mode=0)
        print(f)
        print('Coordinate Reset')

    def selectSLM(self):
        no = self.SLMSpinBox.value()
        mm.qx.selecteorder(no)
        return True
        
    def activateSLM(self):
        mm.qx.activate()
        return True
    
    def deactivateSLM(self):
        mm.qx.deactivate()
        return True
    
    def rlaserpower(self):
        power = self.PowerSpinBox_R.value()
        self.om.ll647.SetPowerLevel(power)
        self.param_update()
    
    def ylaserpower(self):
        power = self.PowerSpinBox_Y.value()
        self.om.ll561.SetPowerLevel(power)
        self.param_update()
        
    def uvlaserpower(self):
        power = self.PowerSpinBox_UV.value()
        self.om.ll405.SetPowerLevel(power)
        self.param_update()
        
    def setattenuat(self):
        pov = self.attenuatorSpinBox.value()
        self.om.att.setattenuation(pov)
        
    def setled(self):
        if (self.LED1Button.isChecked()):
            led = 1
        elif (self.LED2Button.isChecked()):
            led = 2
        else:
            led = False
        return (led)
    
    def piezomove(self):
        newpos = self.piezoSpinBox.value()
        self.om.zst.setPositionf(newpos)
        self.param_update()
        
    def priorleft(self):
        self.om.prior.MoveRelX(-1)
        self.param_update()
        
    def priorright(self):
        self.om.prior.MoveRelX(1)
        self.param_update()
        
    def priorup(self):
        self.om.prior.MoveRelY(-1)
        self.param_update()
        
    def priordown(self):
        self.om.prior.MoveRelY(1)
        self.param_update()
        
    def polarizermove(self):
        newpol = self.polarizerSpinBox.value()
        print('moving polarizer %f, %f' % (newpol, self.om.pol.getCurPos()))
        self.om.pol.MoveAbs(newpol)
        self.param_update()
        
    # def setfilter(self):
    #     n = self.filterspinBox.value()
    #     print('set filter %i' % n)
    #     self.om.fw.setPositionf(n)
        
    def setZM(self):
        mode = self.ZMSpinBox.value()
        amplitude = self.ZM_doubleSpinBox.value()        
        ao.ao_zernike_set(self.om.dm, mode, amplitude, ao.Sm.S)
        return True
    
    def reloadflatfile(self):
        with open(ao.FLAT_FILE, "r") as file:
            cmd = eval(file.readline())
        ao.cmd_best = cmd
        if (all(i <= 1.0 for i in cmd)):
            self.om.dm.Send(cmd)
            print('DM reset to flatfile')
        else:
            raise Exception(' Error: push value greater than 1.0 ')
        return True
        
    def write_dmfile(self):
        ao.writedmfile(self.om.path, self.cmd, self.res)
        return True
        
    def loaddmfile(self):
        fn = self.DMfileEdit.text()
        fns = r'%s' %fn
        with open(fns, "r") as file:
            cmd = eval(file.readline())
        if (all(i <= 1.0 for i in cmd)):
            self.om.dm.Send(cmd)
            print('DM load successfully')
        else:
            raise Exception(' Error: push value greater than 1.0 ')
        return True
        
    def cleardm(self):
        ao.dmnull(self.om.dm)
        return True
         
    def ao_optimize_snr(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        zcenter = self.piezoSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,LED12=led)
        amp_start = self.ao_amp_StartSpinBox.value()
        amp_stop = self.ao_amp_StopSpinBox.value()
        amp_step = self.ao_amp_StepSpinBox.value()
        mode_start = self.ao_mode_StartSpinBox.value()
        mode_stop = self.ao_mode_StopSpinBox.value()
        ampzr = np.arange(amp_start,amp_stop,amp_step)
        modes = np.arange(mode_start,mode_stop,1)
        metric = mt.snr
        lpr = self.LPFSpinBox.value()
        hpr = self.HPFSpinBox.value()
        self.res, self.cmd = ao.ao_optimize_snr(self.om, self.om.dm, self.om.path, zcenter, metric, lpr, hpr, modes, ampzr, ao.Sm.S, gain=setemgain, Blaser=setBlaser, Rlaser=setRlaser, Ylaser=setYlaser, LED=led)
        return True
        
    def ao_optimize_hpf(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        zcenter = self.piezoSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,LED12=led)
        amp_start = self.ao_amp_StartSpinBox.value()
        amp_stop = self.ao_amp_StopSpinBox.value()
        amp_step = self.ao_amp_StepSpinBox.value()
        mode_start = self.ao_mode_StartSpinBox.value()
        mode_stop = self.ao_mode_StopSpinBox.value()
        ampzr = np.arange(amp_start,amp_stop,amp_step)
        modes = np.arange(mode_start,mode_stop,1)
        metric = mt.hf
        self.res, self.cmd = ao.ao_optimize_hf(self.om, self.om.dm, self.om.path, zcenter, metric, modes, ampzr, ao.Sm.S, gain=setemgain, Blaser=setBlaser, Rlaser=setRlaser, Ylaser=setYlaser, LED=led)
        return True

    def ao_optimize_max(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        zcenter = self.piezoSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,LED12=led)
        amp_start = self.ao_amp_StartSpinBox.value()
        amp_stop = self.ao_amp_StopSpinBox.value()
        amp_step = self.ao_amp_StepSpinBox.value()
        mode_start = self.ao_mode_StartSpinBox.value()
        mode_stop = self.ao_mode_StopSpinBox.value()
        ampzr = np.arange(amp_start,amp_stop,amp_step)
        modes = np.arange(mode_start,mode_stop,1)
        metric = mt.peakv
        self.res, self.cmd = ao.ao_optimize_max(self.om, self.om.dm, self.om.path, zcenter, metric, modes, ampzr, ao.Sm.S, gain=setemgain, Blaser=setBlaser, Rlaser=setRlaser, Ylaser=setYlaser, LED=led)
        return True

    def ao_optimize_si(self):
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        zcenter = self.piezoSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,LED12=led)
        amp_start = self.ao_amp_StartSpinBox.value()
        amp_stop = self.ao_amp_StopSpinBox.value()
        amp_step = self.ao_amp_StepSpinBox.value()
        mode_start = self.ao_mode_StartSpinBox.value()
        mode_stop = self.ao_mode_StopSpinBox.value()
        ampzr = np.arange(amp_start,amp_stop,amp_step)
        modes = np.arange(mode_start,mode_stop,1)
        metric = mt.snr
        slm = mm.qx
        self.res, self.cmd = ao.ao_optimize_si(self.om, self.om.dm, self.om.path, zcenter, metric, slm, modes, ampzr, ao.Sm.S, gain=setemgain, Blaser=setBlaser, Rlaser=setRlaser, Ylaser=setYlaser, LED=led)
        return True
        
    def setparameters(self):
        self.WFS.hsp = self.radius_SpinBox.value()
        self.WFS.diameter = self.diameter_SpinBox.value()
        self.WFS.x_center_base = self.Xcenter_spinBox.value()
        self.WFS.y_center_base = self.Ycenter_spinBox.value()
        self.WFS.x_center_offset = self.Xcenter_spinBox_2.value()
        self.WFS.y_center_offset = self.Ycenter_spinBox_2.value()
        self.WFS.px_spacing = self.spacing_doubleSpinBox.value()
        print('SHWF Reconstruction Parameters Set')
        
    def set_baseimage(self):
        exposure_pco = self.exposureSpinBox_2.value()
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setUvlaser = self.L405Button.isChecked()
        led = self.setled()
        self.om.setCMOS_ext(exposure_pco,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,UVlaser=setUvlaser,LED12=led)
        self.om.cam.record(number_of_images=self.om.number_of_images, mode='ring buffer')
        mm.daq.CCDTrig_run(self.om.handleA,self.om.handleB)
        ca = self.om.cam.rec.get_status()
        image, meta = self.om.cam.image(-1)
        self.om.cam.stop()
        self.base = np.rot90(np.flip(image,0),k=1)
        print('SHWFS setup')
        
    def load_baseimage(self):
        fn = self.baseimageEdit.text()
        fns = r'%s' %fn
        self.base = tf.imread(fns)
        print('Reference Set')
 
    def WFrun(self):
        self.md = self.mdspinBox.value()
        exposure_pco = self.exposureSpinBox_2.value()
        setBlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setYlaser = self.L561Button.isChecked()
        setUvlaser = self.L405Button.isChecked()
        led = self.setled()
        self.om.setCMOS_ext(exposure_pco,Blaser=setBlaser,Rlaser=setRlaser,Ylaser=setYlaser,UVlaser=setUvlaser,LED12=led)
        self.om.cam.record(number_of_images=self.om.number_of_images, mode='ring buffer')
        mm.daq.CCDTrig_run(self.om.handleA,self.om.handleB)
        ca = self.om.cam.rec.get_status()
        image, meta = self.om.cam.image(-1)
        self.om.cam.stop()
        self.offset = np.rot90(np.flip(image,0),k=1)
        # measure wavefront
        self.wf = self.WFS.GetAberration2img(self.base, self.offset, self.md)
        # update image
        temp = self.WFS.im[1]
        xh, yh = temp.shape
        self.SHimg = np.zeros((xh, yh, 4),dtype=np.uint8)
        self.SHimg[:,:,3] = 0
        temp = temp - temp.min()
        temp = (255 * (temp / (temp.max()))).astype(np.uint8)
        self.SHimg[:,:,0] = temp
        self.SHimg[:,:,1] = temp
        self.SHimg[:,:,2] = temp
        self.display_shimage(self.SHimg)
        temp = self.wf
        temp = temp - temp.min()
        temp = (255 * (temp / (temp.max()))).astype(np.uint8)
        self.WFimg = np.zeros((self.WFS.diameter, self.WFS.diameter, 4), dtype=np.uint8)
        self.WFimg[:,:,3] = 0
        self.WFimg[:,:,0] = temp
        self.WFimg[:,:,1] = temp
        self.WFimg[:,:,2] = temp
        self.display_wfimage(self.WFimg)
        self.MaxLCD_2.display(self.offset.max())
        self.MinLCD_2.display(self.offset.min())
        self.MaxLCD_3.display(self.wf.max())
        self.MinLCD_3.display(self.wf.min())
        self.StdLCD_3.display(self.wf.std())
        
    def calculate_wf(self):
        self.md = self.mdspinBox.value()
        self.wf = self.WFS.GetAberration2img(self.base, self.offset, method=self.md)
        temp = self.wf
        temp = temp - temp.min()
        temp = (255 * (temp / (temp.max()))).astype(np.uint8)
        self.WFimg = np.zeros((self.WFS.diameter, self.WFS.diameter, 4), dtype=np.uint8)
        self.WFimg[:,:,3] = 0
        self.WFimg[:,:,0] = temp
        self.WFimg[:,:,1] = temp
        self.WFimg[:,:,2] = temp
        self.display_wfimage(self.WFimg)
        self.MaxLCD_3.display(self.wf.max())
        self.MinLCD_3.display(self.wf.min())
        self.StdLCD_3.display(self.wf.std())
        
    def wf_recon(self):
        zn = 64
        d = self.WFS.diameter
        a = np.zeros(zn)
        a = wfzm.decomwf(self.wf, zn)
        if self.ExcludeButton.isChecked():
            a[:4] = 0.
        nwf = wfzm.recomwf(a, d)
        temp = nwf
        temp = temp - temp.min()
        temp = (255 * (temp / (temp.max()))).astype(np.uint8)
        self.WFimg = np.zeros((d, d, 4), dtype=np.uint8)
        self.WFimg[:,:,3] = 0
        self.WFimg[:,:,0] = temp
        self.WFimg[:,:,1] = temp
        self.WFimg[:,:,2] = temp
        self.display_wfimage(self.WFimg)
        self.MaxLCD_3.display(nwf.max())
        self.MinLCD_3.display(nwf.min())
        self.StdLCD_3.display(nwf.std())
        
    def SaveWFimg(self):
        t = time.localtime()
        x = np.array([1e4,1e2,1])
        t1 = int((t[0:3]*x).sum())
        t2 = int((t[3:6]*x).sum())
        tns = "%s_%s_" %(t1,t2)
        pns = self.om.path + '/' + 'Wavefront/'
        try:
            os.mkdir(pns)
        except:
            print('Directory already exists')
        fn1 = os.path.join(pns,tns + 'camera_image.tif')
        fn2 = os.path.join(pns,tns + 'wf_image.tif')
        fn3 = os.path.join(pns,tns + 'sh_image_stack.tif')
        tf.imsave(fn1, self.offset.astype(np.float32),photometric='minisblack')
        tf.imsave(fn3, self.WFS.im.astype(np.float32),photometric='minisblack')
        tf.imsave(fn2, self.wf.astype(np.float32),photometric='minisblack')
        
    def setnewdm(self):
        exclude4 = True
        zn = 64
        if self.ExcludeButton.isChecked():
           exclude4 = True
        self.res, self.cmd = ao.set_corr_wf(self.om, self.om.dm, wfzm, self.wf, zn, ao.Sm.S, exclude4)
        
    def setnewwf(self):
        exclude4 = True
        fn = self.DMfileEdit.text()
        fns = r'%s' %fn
        zn = 64
        if self.ExcludeButton.isChecked():
           exclude4 = True
        self.res, self.cmd = ao.set_new_wf(self.om, self.om.dm, wfzm, fns, zn, ao.Sm.S, exclude4)

    def pushdm(self):
        push = self.Act_doubleSpinBox.value()
        num = self.ActSpinBox.value()
        values = ao.cmd_best
        values[num] = values[num] + push
        if (all(i <= 1.0 for i in values)):
            self.om.dm.Send( values )
            print('Push DM successfully')
        else:
            raise Exception(' Error: push value greater than 1.0 ')
        return True

    def influencefunction(self):
        t = time.localtime()
        x = np.array([1e4,1e2,1])
        t1 = int((t[0:3]*x).sum())
        t2 = int((t[3:6]*x).sum())
        tns = "%s_%s_" %(t1,t2)
        pns = self.om.path + '/' + tns + '_Influencefunction/'
        try:
            os.mkdir(pns)
        except:
            print('Directory already exists')
        push = self.Act_doubleSpinBox.value()
        for i in range(69):
            values = ao.cmd_zero
            self.om.dm.Send( values )
            time.sleep(1.)
            self.set_baseimage()
            values[i] = values[i] + push
            self.om.dm.Send( values )
            time.sleep(1.)
            self.WFrun()
            fn = 'act_%d'%i + '_v%.3f.tif'%push
            fns = os.path.join(pns,fn)
            tf.imsave(fns,self.wf.astype(np.float32),photometric='minisblack')
            fn = 'act_%d'%i + '_gradxy' +  '_v%.3f.tif'%push
            fns = os.path.join(pns,fn)
            tf.imsave(fns,self.WFS.gradxy.astype(np.float32),photometric='minisblack')
        return True
        
    def actuator_amp_scan(self):
        t = time.localtime()
        x = np.array([1e4,1e2,1])
        t1 = int((t[0:3]*x).sum())
        t2 = int((t[3:6]*x).sum())
        tns = "%s_%s_" %(t1,t2)
        pns = self.om.path + '/' + 'Wavefront/'
        try:
            os.mkdir(pns)
        except:
            print('Directory already exists')
        num = self.ActSpinBox.value()
        ran = self.Range_doubleSpinBox.value()
        step = self.Step_doubleSpinBox.value()
        for i in range(int(2*ran/step + 1)):
            values = self.om.dmfile
            if (all(i <= 1.0 for i in values)):
                self.om.dm.Send( values )
            else:
                raise Exception(' Error: push value greater than 1.0 ')
            push = -ran + step*i
            values[i] = values[i] + push
            if (all(i <= 1.0 for i in values)):
                self.om.dm.Send( values )
                print('Push DM successfully')
            else:
                raise Exception(' Error: push value greater than 1.0 ')
            self.WFrun()
            fn1 = 'act_%d'%num + '_v%.3f_wf.tif'%push
            fns1 = os.path.join(pns,tns+fn1)
            fn2 = 'act_%d'%num + '_v%.3f_shimg.tif'%push
            fns2 = os.path.join(pns,tns+fn2)
            tf.imsave(fns1,self.wf.astype(np.float32),photometric='minisblack')
            tf.imsave(fns2 ,self.offset.astype(np.float32),photometric='minisblack')
        self.om.dm.Send( self.om.dmfile )
        print('Done')

if __name__ == "__main__":
    app = QApplication(sys.argv)
    form = Form()
    form.show()
    app.exec_()