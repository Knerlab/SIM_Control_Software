# -*- coding: utf-8 -*-
"""
GUI for microscope. Can be run from IPython interpreter but interpreter must
be set to display images using qt in separate window

@copywrite, Ruizhe Lin and Peter Kner, University of Georgia, 2019
"""

from PyQt4 import uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import pylab

import sys, time, os
import numpy as np

sys.path.append(r'C:/Users/Public/Documents/python_code')
#from pyAndor import pyAndor_P35_a
import micro_P27_vslm1 as mm
import metric_tools_p27 as mt
#import dm_tools_p27 as dt
from thread import start_new_thread
import ao_tools_p33_0 as ao
import qxga_exec_p27 as qx

class Form(QDialog):

    def __init__(self, parent=None):
        super(Form, self).__init__(parent)
        uic.loadUi('VideoMode_slm_sa1.ui', self)
        self.progressBar.setRange(0,100)
        self.index = 0
        self.progressBar.setValue(0)
        self.connect(self.quitButton, SIGNAL("clicked()"), self.quit)
        self.connect(self.startButton, SIGNAL("clicked()"), self.video)
        self.connect(self.stopButton, SIGNAL("clicked()"), self.stopvideo)
        self.connect(self.leftButton, SIGNAL("clicked()"), self.priorleft)
        self.connect(self.rightButton, SIGNAL("clicked()"), self.priorright)   
        self.connect(self.upButton, SIGNAL("clicked()"), self.priorup)   
        self.connect(self.downButton, SIGNAL("clicked()"), self.priordown)
        self.data = np.zeros((512,512)) #np.random.rand(64,64)
        self.imgh = self.mainplot.axes.imshow(self.data,interpolation='nearest',vmin = 1000,vmax = 2000 )
        ## microscope ##
        self.om = mm.scope(self)
        #self.om.fw.setPositionf(1)
        self.connect(self.om,SIGNAL('update'),self.updateimage)
        self.temperatureLCD.display(self.om.ccd.GetTemperature())
        ## spin box
        self.connect(self.piezoSpinBox, SIGNAL("valueChanged(double)"), self.piezomove)
        temppos = self.om.zst.getPosition()
        self.piezoSpinBox.setValue(temppos)
        ## Taking Stacks
        self.StartSpinBox.setValue(temppos-2.0)
        self.StopSpinBox.setValue(temppos+2.0)
        self.StepSpinBox.setValue(0.2)
        #take image stack
        self.connect(self.simzstackButton, SIGNAL("clicked()"), self.take_z_stack)
        self.connect(self.sectioning_zstackButton, SIGNAL("clicked()"), self.take_sec_z_stack)
        self.connect(self.stackButton, SIGNAL("clicked()"), self.takestack)
        ## Save Stack
        self.connect(self.saveButton, SIGNAL("clicked()"), self.savestack)
        #time-stack
        self.connect(self.TimeStackButton, SIGNAL("clicked()"), self.take_t_stack)
        # set DM
#        dmflatfile = 'C:/Users/Public/Documents/python_code/AdaptiveOptics/flatfile_20180131.mro'
#        dmflat,t = self.om.dm.read_command_file(dmflatfile)
#        self.om.dm.set_mirror(dmflat)
        # polarizer
        self.connect(self.polarizerSpinBox, SIGNAL("valueChanged(double)"), self.polarizermove)
        #rlaser power
        self.connect(self.PowerSpinBox, SIGNAL("valueChanged(double)"), self.rlaserpower)
        #filter
        self.connect(self.FilterButton, SIGNAL("clicked()"), self.setfilter)
        #setdm
        self.connect(self.pushDMButton, SIGNAL("clicked()"), self.setzernike)
        self.connect(self.resetDMButton, SIGNAL("clicked()"), self.reloadflatfile)
#        self.connect(self.loadDMButton, SIGNAL("clicked()"), self.loaddmfile)
        self.connect(self.MetricCalculationButton, SIGNAL("clicked()"), self.update_metric_values)
        self.connect(self.ZernikeScanButton, SIGNAL("clicked()"), self.ao_manual)
        self.connect(self.SelectOrderButton, SIGNAL("clicked()"), self.selectSLM)
        self.connect(self.ActivateSLMButton, SIGNAL("clicked()"), self.activateSLM)
        self.connect(self.DeactivateSLMButton, SIGNAL("clicked()"), self.deactivateSLM)
        self.connect(self.AOButton, SIGNAL("clicked()"), self.ao_optimize)
        self.connect(self.SaveAOButton, SIGNAL("clicked()"), self.write_dmfile)
        
    def quit(self):
        #self.thread1.run()
        self.close()
        
    def run(self):
        self.show()
        self.raise_()
        self.exec_()

    def savestack(self):
        text = self.fileEdit.text()
        self.om.saveTifA(slideName=text)
    
    def takestack(self):
        setlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        setstart = self.StartSpinBox.value()
        setstop = self.StopSpinBox.value()
        setstep = self.StepSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,laser=setlaser,Rlaser=setRlaser,LED12=led)
        self.om.StackExt(setstart,setstop,setstep)

    def take_z_stack(self):
        setlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        setstart = self.StartSpinBox.value()
        setstop = self.StopSpinBox.value()
        setstep = self.StepSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,laser=setlaser,Rlaser=setRlaser,LED12=led)
        self.om.Stack_Patterns(setstart,setstop,setstep)
        
    def take_sec_z_stack(self):
        setlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        setstart = self.StartSpinBox.value()
        setstop = self.StopSpinBox.value()
        setstep = self.StepSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,laser=setlaser,Rlaser=setRlaser,LED12=led)
        self.om.Stack_Sectioning(setstart,setstop,setstep)
        
    def take_t_stack(self):
        setlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        led = self.setled()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        time = self.TimeSpinBox.value()
        pol = self.polarizerSpinBox.value()
        angle = self.angle_spinBox.value()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,laser=setlaser,Rlaser=setRlaser,LED12=led)
        self.om.Image_Patterns(angle, time, pol)
    
    def param_update(self):
        if (not self.om.is_videomode):
            self.temperatureLCD.display(self.om.ccd.GetTemperature())
        self.exposureLCD.display(self.om.ccd.GetExposureTime())
        self.piezoLCD.display(self.om.zst.getPosition())
        self.redlaserpowerLCD.display(self.om.ll647.GetPower())
        self.om.prior.getPosition()
    
    def video(self):
        if self.cmapButton.isChecked():
            self.imgh.set_cmap(pylab.cm.bone)
        else:
            self.imgh.set_cmap(pylab.cm.jet)
        self.param_update()
        setlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        led = self.setled()
        setFTM = self.FTMButton.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        self.om.normal_mode = self.Mode_selectButton.isChecked()
        self.om.open_Acq(exposure=setexposure,emgain=setemgain,laser=setlaser,Rlaser=setRlaser,
                         cntLaser=False,LED12=led,FTM=setFTM,conv=False,ccd=True,
                         trig=1,regen=False)
        self.om.start()
        
    def stopvideo(self):
        self.om.is_videomode = False
        time.sleep(0.01)
        self.param_update()
        
    def updateimage(self):
        self.progressBar.setValue(self.index)
        self.index = (self.index + 1) % 100
        self.data = self.om.data
        self.imgh.set_data(self.data)
        self.imgh.set_clim(self.data.min(),self.data.max())
        self.mainplot.figure.canvas.draw()
        self.intensityBar.setValue(self.data.max())
        QApplication.processEvents()
        
    def update_metric_values(self):
        img = self.data
        #calculate
#        self.peak = mt.peak(img)
#        self.std = mt.std(img)
        self.snr = mt.snr(img)
#        self.rms = mt.rms(img)
        self.sharp = mt.sharp(img)
        #display
#        self.peak_display.display(self.peak)
#        self.std_display.display(self.std)
        self.SNRLCD.display(self.snr)
#        self.rms_display.display(self.rms)
        self.SHARPLCD.display(self.sharp)
        return True
        
    def selectSLM(self):
        no = self.SLMSpinBox.value()
        qx.selecteorder(no)
        return True
        
    def activateSLM(self):
        qx.activate()
        return True
    
    def deactivateSLM(self):
        qx.deactivate()
        return True
    
    def rlaserpower(self):
        power = self.PowerSpinBox.value()
        self.om.ll647.SetPowerLevel(power)
        self.param_update()
        
    def setled(self):
        if (self.LED1Button.isChecked()):
            led = 1
        elif (self.LED2Button.isChecked()):
            led = 2
        else:
            led = False
        return (led)
    
    def piezomove(self):
        #print("moving ")
        newpos = self.piezoSpinBox.value()
        self.om.zst.setPositionf(newpos)
        self.param_update()
        
    def priorleft(self):
        self.om.prior.MoveRelX(-5)
        self.param_update()
        
    def priorright(self):
        self.om.prior.MoveRelX(5)
        self.param_update()
        
    def priorup(self):
        self.om.prior.MoveRelY(-5)
        self.param_update()
        
    def priordown(self):
        self.om.prior.MoveRelY(5)
        self.param_update()
        
    def polarizermove(self):
        newpol = self.polarizerSpinBox.value()
        print('moving polarizer %f, %f' % (newpol, self.om.pol.getCurPos()))
        self.om.pol.MoveAbs(newpol)
        self.param_update()
        
    def setfilter(self):
        n = self.filterspinBox.value()
        print('set filter %i' % n)
        self.om.fw.setPositionf(n)
        
    def setzernike(self):
        mode = np.arange(4,11,1)
        amplitude = np.zeros((7))
        amplitude[0] = self.z4_doubleSpinBox.value()
        amplitude[1] = self.z5_doubleSpinBox.value()
        amplitude[2] = self.z6_doubleSpinBox.value()
        amplitude[3] = self.z7_doubleSpinBox.value()
        amplitude[4] = self.z8_doubleSpinBox.value()
        amplitude[5] = self.z9_doubleSpinBox.value()
        amplitude[6] = self.z10_doubleSpinBox.value()
        self.dmfile = ao.setZernike(self.om, self.om.dm, mode, amplitude, ao.Sm.S)
        return True
    
    def reloadflatfile(self):
        with open(ao.FLAT_FILE, "r") as file:
            cmd = eval(file.readline())
        if (all(i <= 1.0 for i in cmd)):
            self.om.dm.Send(cmd)
            print('DM reset to flatfile')
        else:
            raise Exception(' Error: push value greater than 1.0 ')
        return True
      
    def ao_manual(self):
        ps = self.om.zst.getPosition()
        setlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,laser=setlaser,Rlaser=setRlaser,LED12=led)
        t=time.localtime()
        dpth=str(int(t[0]*1e4+t[1]*1e2+t[2]))+'_'+str(int(t[3]*1e4+t[4]*1e2+t[4]))
        newpath=os.path.join(self.om.path,dpth)
        try:
            os.mkdir(newpath)
        except:
            print('Directory already exists')
            
        start = self.StartSpinBox_z.value()
        stop = self.StopSpinBox_z.value()
        step = self.StepSpinBox_z.value()
        amp = np.arange(start,stop,step)

        if self.Mode_selectButton_z4.isChecked():
            mode = 4
            for counter, v in enumerate(amp):
                print('zernike_mode_#%d'%mode,counter)
                ao.ao_zernike_scan(self.om, self.om.dm, mode, v, ao.Sm.S)
                self.om.zst.setPositionf(ps)
                self.om.singleSnapExt(verbose=False)
                self.om.saveTifA('widefield_image_with_zerinkemode#%d_amp%f'%(mode,v))
        if self.Mode_selectButton_z5.isChecked():
            mode = 5
            for counter, v in enumerate(amp):
                print('zernike_mode_#%d'%mode,counter)
                ao.ao_zernike_scan(self.om, self.om.dm, mode, v, ao.Sm.S)
                self.om.zst.setPositionf(ps)
                self.om.singleSnapExt(verbose=False)
                self.om.saveTifA('widefield_image_with_zerinkemode#%d_amp%f'%(mode,v))
        if self.Mode_selectButton_z6.isChecked():
            mode = 6
            for counter, v in enumerate(amp):
                print('zernike_mode_#%d'%mode,counter)
                ao.ao_zernike_scan(self.om, self.om.dm, mode, v, ao.Sm.S)
                self.om.zst.setPositionf(ps)
                self.om.singleSnapExt(verbose=False)
                self.om.saveTifA('widefield_image_with_zerinkemode#%d_amp%f'%(mode,v))
        if self.Mode_selectButton_z7.isChecked():
            mode = 7
            for counter, v in enumerate(amp):
                print('zernike_mode_#%d'%mode,counter)
                ao.ao_zernike_scan(self.om, self.om.dm, mode, v, ao.Sm.S)
                self.om.zst.setPositionf(ps)
                self.om.singleSnapExt(verbose=False)
                self.om.saveTifA('widefield_image_with_zerinkemode#%d_amp%f'%(mode,v))
        
        self.om.dm.Send(self.om.dmfile)
        print('Finish scan')
        return True
        
    def ao_optimize(self):
        setlaser = self.L488Button.isChecked()
        setRlaser = self.L647Button.isChecked()
        setemgain = self.emgainSpinBox.value()
        setexposure = self.exposureSpinBox.value()
        led = self.setled()
        self.om.setCCD_EM_ext(exposure=setexposure,emgain=setemgain,laser=setlaser,Rlaser=setRlaser,LED12=led)
        amp_start = self.ao_amp_StartSpinBox.value()
        amp_stop = self.ao_amp_StopSpinBox.value()
        amp_step = self.ao_amp_StepSpinBox.value()
        mode_start = self.ao_mode_StartSpinBox.value()
        mode_stop = self.ao_mode_StopSpinBox.value()
        ampzr = np.arange(amp_start,amp_stop,amp_step)
        modes = np.arange(mode_start,mode_stop,1)
        if self.SharpMetricButton.isChecked():
            metric = mt.sharp
        if self.SNRMetricButton.isChecked():
            metric = mt.snr
        self.res, self.cmd = ao.ao_optimize(self.om, self.om.dm, metric, modes, ampzr, ao.Sm.S, gain=setemgain, blaser=setlaser, rlaser=setRlaser, LED=led)
        return True
        
    def write_dmfile(self):
        ao.writedmfile(self.cmd, self.res)
        return True
        
#    def savedmfile(self):
#        dt.writedmfile(self.dmfile)
#        return True

if __name__ == "__main__":
    app = QApplication(sys.argv)
    form = Form()
    form.show()
    app.exec_()
