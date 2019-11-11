# -*- coding: utf-8 -*-
"""
Latest version of microscope code.
New 64-bit computer, python 2.7

@copywrite, Ruizhe Lin and Peter Kner, University of Georgia, 2019
"""

datapath = r'C:/Users/Public/Documents/data'
basepath = r'C:/Users/Public/Documents/python_code'
import os,sys,time,random
try:
    sys.path.index(basepath)
except:
    sys.path.append(basepath)
    sys.path.append(basepath+'\\mirao')

import numpy as N
import tifffile as T
from DAQ27 import daq
from pyAndor27 import pyAndor_P271_slm
import pol_rotator27.motorized_pol_v2 as mp2
from prior27 import zstageP27 as zstage
from prior27 import priorxyP27
import filterWheel_v2 as filterWheel
import Coherent_P33
from asdk import DM
import qxga_exec_p27 as qx
#
from getpass import getuser
from PyQt4.QtCore import QThread, SIGNAL
#
FREQ = 10000 # for daq timing signals
FLAT_FILE = r'C:/Users/Public/Documents/python_code/SIM_P27/20190408140547_flatfile.txt'
TEMPERATURE_SETPOINT = -50
COOLER_MODE = 1 # 1 for stays on on shutdown, 0 for turns off.

class scope(QThread):

    def __init__(self,parent=None):
        super(scope, self).__init__()
        self.delay = 0.02 # delay in loops (StackExt)
        self.handleA = None
        self.handleB = None
        self.QUIT = False
        self.stackparams = {'Date/Time':0,'X':0,'Y':0,'Z1':0,'Z2':0,'Zstep':0,'Exposure(s)':0,'CCD Temperature':0,'Pixel Size(nm)':89,'CCD setting':'','User':''}
        # camera
        self.ccd = pyAndor_P271_slm.ccd()
        self.ccd.CoolerON()
        self.ccd.SetCoolerMode(COOLER_MODE)
        self.ccd.SetTemperature(TEMPERATURE_SETPOINT)
        self.data = N.zeros(self.ccd.image_size, dtype=N.uint16)
        # piezo
        self.zst = zstage.zstage()
        self.zst.setPositionf(50)
        self.zpos = 0
        # Filter wheel
        self.fw = filterWheel.filter_wheel()
        # deformable mirror
        self.serialName = 'BAX228'
        self.dm = DM( self.serialName )
        self.nbAct = int( self.dm.Get('NBOfActuator') )
#        self.values = [0.] * self.nbAct
#        self.dm.Send( self.values )
        with open(FLAT_FILE, "r") as file:
            self.dmfile = eval(file.readline())
        if (all(i <= 1.0 for i in self.dmfile)):
            self.dm.Send(self.dmfile)
        else:
            raise Exception(' Error: push value greater than 1.0 ')
        #motorized polarizer
        print('Initializing Polarizer')
        self.pol = mp2.A8MRU()
        self.pol.MotorOn()
        self.pol.setSpeed(1000,1)
        self.pol_array = N.array([50.,20.,-10.]) #3D
        print(self.pol.getState())
        # priorxy
        self.prior = priorxyP27.prior()
        self.xpos = 0
        self.ypos = 0
        self.prior.limits=[-10000, 10000, -10000, 10000]
        # qxga
        qx.initiate()
        qx.open_usb_port()
        #Coherent 647 laser
        try:
            self.ll647 = Coherent_P33.obis('COM16')
            self.ll647.SetLaserOn()
            self.ll647.SetDigitalModMode()
        except:
            self.ll647 = None
            print('647nm Laser not on') #raise Exception('647nm Laser not on')
        #set save path        
        t=time.localtime()
        self.is_videomode = False
        dpth=str(int(t[0]*1e4+t[1]*1e2+t[2]))+'_'+getuser()
        self.path=os.path.join(datapath,dpth)
        try:
            os.mkdir(self.path)
        except:
            print('Directory already exists')
        # defaults
        self.coord=(200,400,200,400)
        self.Name='default'
        ###
        self.Name = 'default'
        self.coordVal = self.ccd.GetImageCoordinates()
        self.old_Cval = self.ccd.GetImageCoordinates()
        self.normal_mode = False
          
    def __del__(self):
        if (self.handleA != None):
            daq.CCDTrig_close(self.handleA,self.handleB)
        x = self.dm.Reset()
        if (x==0):
            print('DM Reset')
        else:
            print('DM cannot reset')
        qx.close()
#        self.dm.close()
#    
    def get_img(self):
        daq.CCDTrig_run(self.handleA,self.handleB)
        self.ccd.WaitForNewData()
        self.data[:,:] = self.ccd.images
        
    def open_Acq(self,exposure=0.1,emgain=200,laser=False,Rlaser=False,cntLaser=False,LED12=1,FTM=False,conv=False,ccd=True,trig=1,regen=False):
        # setup camera
        self.ccd.SetTriggerMode(trig) #0:internal, 1:External
        self.ccd.SetShutterMode(2) #auto
        self.ccd.SetReadMode(4) #image
        self.ccd.SetADChannel(0) # 14-bit channel
        self.ccd.SetOutputAmplifier(int(conv))
        self.ccd.SetHSSpeed(0)
        self.ccd.SetEMCCDGain(emgain)
        if FTM: # Frame transfer mode
            self.ccd.SetExposureTime(0)
            self.ccd.SetFrameTransferMode(1)
            self.ccd.SetAcquisitionMode(7) # run until abort
        else:
            self.ccd.SetExposureTime(exposure)
            self.ccd.SetFrameTransferMode(0)
            self.ccd.SetAcquisitionMode(5) # run until abort
        q,w,e = self.ccd.GetAcquisitionTimings()
        # setup triggers
        if not self.handleA == None:
            daq.CCDTrig_close(self.handleA,self.handleB)
        if cntLaser:
            digsig = self.getTiming3(FREQ,1e3*exposure,1e3*exposure,1,laser,Rlaser,LED=LED12,CCD=ccd)
        elif regen:
            digsig = self.getTiming2(FREQ,1e3*exposure,1e3*exposure,1,laser,Rlaser,LED=LED12,CCD=ccd)
        else:
            digsig = self.getTiming(FREQ,1e3*exposure,1e3*exposure,1,laser,Rlaser,LED=LED12,CCD=ccd)
        b = daq.CCDTrig_open(FREQ,digsig)
        self.handleA = b[0]
        self.handleB = b[1]
        xs = self.ccd.image_size[0]
        ys = self.ccd.image_size[1]
        self.data = N.zeros((xs,ys), dtype=N.uint16)
        self.ccd.SetShutterMode(1)
        self.ccd.Acquire()
        #self.slm.SLM_on()
        print(q,w,e)
        return (q,w,e)
    
    def close_Acq(self):        
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
        return True
        
    def run(self):
        self.is_videomode = True
        if (self.normal_mode==True):
            while (self.is_videomode):
                self.get_img()
                self.emit(SIGNAL('update'))
        else:
            while (self.is_videomode):
                #self.slm.show_next_patt()
#                psz = qx.getordernum()
                psz =15
                phs = int(psz/3)
                for m in range(3):
                    self.pol.MoveAbs(self.pol_array[m])
                    time.sleep(0.02)
                    for n in range(phs):
                        qx.selecteorder(m*phs+n)
                        qx.activate()
                        self.get_img()
                        self.emit(SIGNAL('update'))
                        qx.deactivate()
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
        #self.slm.SLM_off()
        
############## External ######################################################

    def StackExt(self,start,stop,step=0.2,verbose=True):
        init_loc=self.zst.getPosition()
        no = int((stop-start)/step)+1
        pos = start
        xs = self.ccd.image_size[0]
        ys = self.ccd.image_size[1]
        self.data = N.zeros((no,xs,ys), N.uint16)
        self.ccd.SetShutterMode(1)
        q = self.ccd.Acquire()
        time.sleep(0.2) # was 0.05
        for p in range(no):
            self.zst.setPositionf(pos)
            daq.CCDTrig_run(self.handleA,self.handleB)
            q = self.ccd.WaitForNewData()
            print(p,q)
            self.data[p] = self.ccd.images
            pos += step
            time.sleep(self.delay)
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
        if verbose:
            T.imshow(self.data, vmin=self.data.min(), vmax=self.data.max())
        cur_pos = self.prior.getPosition()
        self.stackTags(cur_pos[0],cur_pos[1],start,stop,step,function='Z-Stack')
        self.zst.setPositionf(init_loc)
        return True
        
    def Stack_Patterns(self,start,stop,step=0.2, verbose=True):
        rots = self.pol_array
        no = int((stop-start)/step)+1
        pos = start
        xs = self.ccd.image_size[0]
        ys = self.ccd.image_size[1]
#        psz = qx.getordernum()
        psz = 15
        phs = int(psz/3)
        self.data = N.zeros((psz*no,xs,ys), dtype=N.uint16)
#        self.slm.SLM_on()
        self.ccd.SetShutterMode(1)
        q = self.ccd.Acquire()
        time.sleep(0.1) # was 0.2,  changed 20141114
        for p in range(no):
            self.zst.setPositionf(pos)
            for w in range(3):
                self.pol.MoveAbs(rots[w])
                time.sleep(0.4)
                for m in range(phs):
                    #self.pr.setVoltage(rots[m])
                    qx.selecteorder(phs*w+m)
                    qx.activate()
                    time.sleep(0.02)
                    #print self.pr.getVoltage()
                    #self.dmd.set_image(patt)
                    #self.slm.show_next_patt()
                    #self.slm.show_patt(m)
                    daq.CCDTrig_run(self.handleA,self.handleB)
                    q = self.ccd.WaitForNewData()
                    print (p,q)
                    self.data[psz*p + 5*w + m] = self.ccd.images
                    qx.deactivate()
#                    time.sleep(0.02)
            pos += step
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
#        self.slm.SLM_off()
        #self.pr.setVoltage(0.0)
        if verbose:
            T.imshow(self.data, vmin=self.data.min(), vmax=self.data.max())
        cur_pos = self.zst.getPosition()
        self.stackTags(cur_pos,start,stop,step,function='Z-Stack patterns')
        return True
        
    def Stack_Sectioning(self,start,stop,step=0.2, verbose=True):
        no = int((stop-start)/step)+1
        pos = start
        xs = self.ccd.image_size[0]
        ys = self.ccd.image_size[1]
#        psz = qx.getordernum()
        psz = 3
        self.data = N.zeros((psz*no,xs,ys), dtype=N.uint16)
        self.ccd.SetShutterMode(1)
        q = self.ccd.Acquire()
        self.pol.MoveAbs(0)
        time.sleep(0.4)
        for p in range(no):
            self.zst.setPositionf(pos)
            for m in range(psz):
                qx.selecteorder(15+m)
                qx.activate()
                time.sleep(0.02)
                daq.CCDTrig_run(self.handleA,self.handleB)
                q = self.ccd.WaitForNewData()
                print (p,q)
                self.data[psz*p + m] = self.ccd.images
                qx.deactivate()
                time.sleep(0.02)
            pos += step
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
        if verbose:
            T.imshow(self.data, vmin=self.data.min(), vmax=self.data.max())
        cur_pos = self.zst.getPosition()
        self.stackTags(cur_pos,start,stop,step,function='Z-Stack patterns')
        return True

    def Image_Patterns(self, angle=0, no=200, pol=0, verbose=True):
        pos = self.zst.getPosition()
        xs = self.ccd.image_size[0]
        ys = self.ccd.image_size[1]
#        psz = self.slm.ni
        psz = 5
        self.data = N.zeros((psz*no,xs,ys), dtype=N.uint16)
#        self.slm.SLM_on()
        self.ccd.SetShutterMode(1)
        q = self.ccd.Acquire()
        time.sleep(0.01) # was 0.2,  changed 20141114
        self.zst.setPositionf(pos)
        for p in range(no):  
            for m in range(angle*5,psz+angle*5):
                #self.pr.setVoltage(rots[m])
                self.pol.MoveAbs(pol)
                qx.selecteorder(m)
                qx.activate()
                #time.sleep(.50)
                #print self.pr.getVoltage()
                #self.dmd.set_image(patt)
                #self.slm.show_next_patt()
#                self.slm.show_patt(m)
                daq.CCDTrig_run(self.handleA,self.handleB)
                q = self.ccd.WaitForNewData()
                print (p,q)
                self.data[psz*p+m%5] = self.ccd.images
                qx.deactivate()
                time.sleep(self.delay)
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
#        self.slm.SLM_off()
        #self.pr.setVoltage(0.0)
        if verbose:
            T.imshow(self.data, vmin=self.data.min(), vmax=self.data.max())
        cur_pos = self.prior.getPosition()
        self.stackTags(cur_pos[0],cur_pos[1],function='Z-Stack patterns')
        return True
        
    def si_2d_pattern(self,start,stop,step=0.2, verbose=True):
        rots = self.pol_array
        no = int((stop-start)/step)+1
        pos = start
        xs = self.ccd.image_size[0]
        ys = self.ccd.image_size[1]
#        psz = qx.getordernum()
        psz = 9
        phs = int(psz/3)
        self.data = N.zeros((psz*no,xs,ys), dtype=N.uint16)
#        self.slm.SLM_on()
        self.ccd.SetShutterMode(1)
        q = self.ccd.Acquire()
        time.sleep(0.1) # was 0.2,  changed 20141114
        for p in range(no):
            self.zst.setPositionf(pos)
            for w in range(3):
                self.pol.MoveAbs(rots[w])
                time.sleep(0.8)
                for m in range(phs):
                    #self.pr.setVoltage(rots[m])
                    qx.selecteorder(19+phs*w+m)
                    qx.activate()
                    time.sleep(0.02)
                    #print self.pr.getVoltage()
                    #self.dmd.set_image(patt)
                    #self.slm.show_next_patt()
                    #self.slm.show_patt(m)
                    daq.CCDTrig_run(self.handleA,self.handleB)
                    q = self.ccd.WaitForNewData()
                    print (p,q)
                    self.data[psz*p + 5*w + m] = self.ccd.images
                    qx.deactivate()
#                    time.sleep(0.02)
            pos += step
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
#        self.slm.SLM_off()
        #self.pr.setVoltage(0.0)
        if verbose:
            T.imshow(self.data, vmin=self.data.min(), vmax=self.data.max())
        cur_pos = self.zst.getPosition()
        self.stackTags(cur_pos,start,stop,step,function='Z-Stack patterns')
        return True        
        
    def setlaserpower(self,power):
        self.ll647.SetPowerLevel(power)
        out = self.ll647.GetPower()
        return out
    
    def singleSnapExt(self,verbose=True):
        xs = self.ccd.image_size[0]
        ys = self.ccd.image_size[1]
        self.data = N.zeros((xs,ys), dtype=N.uint16)
        self.ccd.SetShutterMode(1)
        self.ccd.Acquire()
        time.sleep(0.2) # was 0.05
        daq.CCDTrig_run(self.handleA,self.handleB)
        self.ccd.WaitForNewData()
        self.data[:,:] = self.ccd.images
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
        if verbose:
            T.imshow(self.data, vmin=self.data.min(), vmax=self.data.max())
        cur_pos = self.prior.getPosition()
        self.stackTags(cur_pos[0],cur_pos[1],self.zst.getPosition(),0,0,function='Single snap')
        return True
        
    def setCCD_Conv_ext(self,exposure=0.100,laser=False,Rlaser=False,LED12=1):
        self.is_videMode = False
        # setup camera
        self.ccd.SetTriggerMode(1) #0: internal
        self.ccd.SetShutterMode(2) #0: auto
        self.ccd.SetFrameTransferMode(0)
        self.ccd.SetAcquisitionMode(5) # run until abort
        self.ccd.SetReadMode(4) #image
        self.ccd.SetADChannel(0) # 14-bit channel
        # this sets camera to 1.8us vs speed, and 3MHz readout which doesn't have a lot of pattern
        # to get rid of patterns entirely in conventional mode, go to 1Mhz readout
        self.ccd.SetVSSpeed(3)
        self.ccd.SetHSSpeed(0)
        self.ccd.SetOutputAmplifier(1)
        self.ccd.SetExposureTime(exposure)
        q,w,e = self.ccd.GetAcquisitionTimings()
        # setup triggers
        if not self.handleA == None:
            daq.CCDTrig_close(self.handleA,self.handleB)
        digsig = self.getTiming(FREQ,1e3*exposure,1e3*exposure,1,laser,Rlaser,LED=LED12)
        b = daq.CCDTrig_open(FREQ,digsig)
        self.handleA = b[0]
        self.handleB = b[1]
        return (q,w,e)

    def setCCD_EM_ext(self,exposure=0.1,emgain=200,laser=False,Rlaser=False,LED12=1):
        self.is_videMode = False
        # setup camera
        self.ccd.SetTriggerMode(1) #0:internal, 1:External
        self.ccd.SetShutterMode(2) #0: auto
        self.ccd.SetFrameTransferMode(0)
        self.ccd.SetAcquisitionMode(5) # run until abort
        self.ccd.SetReadMode(4) #image
        self.ccd.SetADChannel(0) # 14-bit channel
        self.ccd.SetOutputAmplifier(0)
        self.ccd.SetHSSpeed(0)
        self.ccd.SetEMCCDGain(emgain)
        self.ccd.SetExposureTime(exposure)
        q,w,e = self.ccd.GetAcquisitionTimings()
        # setup triggers
        if not self.handleA.value == 0:
            daq.CCDTrig_close(self.handleA,self.handleB)
        digsig = self.getTiming(FREQ,1e3*exposure,1e3*exposure,1,laser,Rlaser,LED=LED12)
        b = daq.CCDTrig_open(FREQ,digsig)
        self.handleA = b[0]
        self.handleB = b[1]
        return (q,w,e)
        
    def getTiming(self,freq,exposure,pulse,delay,laser=False,Rlaser=False,LED=1,CCD=True):
        ''' frequency is in Hertz, exposure and delay are in ms
            bit 1: CCD Trigger
            bit 2: LED pulse
            bit 3: Laser shutter
            bit 4: Red Laser
            bit 5: LED 1/2'''
        # laser shutter is active
        count = int(1.1*(exposure+delay)*0.001*freq)
        td = int(delay*0.001*freq)
        texp = int(exposure*0.001*freq)
        tpulse = int(pulse*0.001*freq)
        trigpulse = int(0.004*freq)
        p = range(count)
        laser_arr = N.zeros(count)
        led_arr = N.zeros(count)
        Rlaser_arr = N.zeros(count)
        if LED==1:
            LED12_arr = N.zeros(count)
        else:
            LED12_arr = N.ones(count)
        if laser:
            laser_arr = N.array([(i>td) & (i<(td+texp)) for i in p])
        elif Rlaser:
            Rlaser_arr = N.array([(i>td) & (i<(td+texp)) for i in p])
        else:
            led_arr = N.array([(i>td) & (i<(td+tpulse)) for i in p])
        if (CCD):
            ccdtrig = N.array([(i<=(trigpulse)) for i in p])
        else:
            ccdtrig = N.array([(0) for i in p])
        #out = N.add(N.add(1.*ccdtrig,2.*lha),4.*shutt2).astype(N.int)
##        out = (ccdtrig+2*led+4*laser).astype(N.int)
        out = (ccdtrig+2*led_arr+4*laser_arr+8*Rlaser_arr+16*LED12_arr).astype(N.int)
        return out
    
    def getTiming2(self,freq,exposure,pulse,delay,laser=False,Rlaser=False,LED=2,CCD=True):
        ''' Regen with every pulse with LED           
            exposure time is for laser, followed by pulse time that is for LED.
            '''
        # laser shutter is active
        count = int(1.1*(exposure+pulse+2*delay)*0.001*freq)
        td = int(delay*0.001*freq)
        texp = int(exposure*0.001*freq)
        tpulse = int(pulse*0.001*freq)
        trigpulse = int(0.004*freq)
        p = range(count)
        laser_arr = N.zeros(count)
        led_arr = N.zeros(count)
        Rlaser_arr = N.zeros(count)
        if LED==1:
            LED12_arr = N.zeros(count)
        else:
            LED12_arr = N.ones(count)
        if laser:
            laser_arr = N.array([(i>td) & (i<(td+texp)) for i in p])
        elif Rlaser:
            Rlaser_arr = N.array([(i>td) & (i<(td+texp)) for i in p])
        else:
            led_arr = N.array([(i>(2*td+texp)) & (i<(2*td+texp+tpulse)) for i in p])
        if (CCD):
            ccdtrig = N.array([(i<=(trigpulse)) for i in p])
        else:
            ccdtrig = N.array([(0) for i in p])
        out = (ccdtrig+2*led_arr+4*laser_arr+8*Rlaser_arr+16*LED12_arr).astype(N.int)
        return out

    def getTiming3(self,freq,exposure,pulse,delay,laser=False,Rlaser=False,LED=1,CCD=True):
##    def getTiming(self,freq,exposure,pulse,delay,laser=False):
        ''' frequency is in Hertz, exposure and delay are in ms
            bit 1: CCD Trigger
            bit 2: LED pulse
            bit 3: Laser shutter Contineous
            bit 4: Red Laser contineous
            bit 5: LED 1/2'''
        # laser shutter is active
        count = int(1.1*(exposure+delay)*0.001*freq)
        td = int(delay*0.001*freq)
        texp = int(exposure*0.001*freq)
        tpulse = int(pulse*0.001*freq)
        trigpulse = int(0.004*freq)
        p = range(count)
        led_arr = N.zeros(count)
        ccdtrig = N.zeros(count)
        laser_arr = N.zeros(count)
        Rlaser_arr = N.zeros(count)
        LED12_arr = N.ones(count)
        if LED==1:
            LED12_arr = N.zeros(count)
        if laser:
            laser_arr = N.ones(count)
        if Rlaser:
            Rlaser_arr = N.ones(count)
        if not (laser | Rlaser):
#            led_arr = N.array([(i>td) & (i<(td+tpulse)) for i in p])
            led_arr = N.ones(count)
        if (CCD):
            ccdtrig = N.array([(i<=(trigpulse)) for i in p])
        out = (ccdtrig+2*led_arr+4*laser_arr+8*Rlaser_arr+16*LED12_arr).astype(N.int)
        return out
        
####### Internal Trigger ###################################################
      
    def setCCD_Conv_Int(self,exposure=0.100):
        self.ccd.SetAcquisitionMode(1) # single exposure
        self.ccd.SetOutputAmplifier(1)
        self.ccd.SetTriggerMode(0) #0: internal
        self.ccd.SetShutterMode(2) #0: auto
        self.ccd.SetExposureTime(exposure)
        q = self.ccd.GetAcquisitionTimings()
        return q
        
    def StackInt(self,start,stop,step=0.2):
        no = int((stop-start)/step)+1
        pos = start
        xs = self.ccd.image_size[0]
        ys = self.ccd.image_size[1]
        self.data = N.zeros((no,xs,ys), dtype=N.uint16)
        self.ccd.SetShutterMode(1)
        for p in range(no):
            self.zst.setPositionf(pos)
            q = self.ccd.Acquire()
            q = self.ccd.WaitForNewData()
            q = self.ccd.AbortAcquisition()
            self.data[p] = self.ccd.images
            pos += step
        self.ccd.AbortAcquisition()
        self.ccd.SetShutterMode(2)
        cur_pos = self.prior.getPosition()
        self.stackTags(cur_pos[0],cur_pos[1],start,stop,step,function='Z-Stack')
        T.imshow(self.data, vmin=self.data.min(), vmax=self.data.max())
        return True
        
###########################################################################
        
    def saveTifA(self,slideName='',comments='',Upload=False):
        t=time.localtime()
        x = N.array([1e4,1e2,1])
        t1 = int((t[0:3]*x).sum())
        t2 = int((t[3:6]*x).sum())
        if slideName=='':
            slideName=self.Name
        else:
            self.Name=slideName
        self.stackparams['Comments']=comments
        self.stackparams['Slide Name']=slideName
        fn = "%s-%s_%s_%s" %(t1,t2,slideName,comments)
        fn1 = os.path.join(self.path,fn+'.tif')
        fn2 = os.path.join(self.path,fn+'_ps.txt')
        T.imsave(fn1,self.data)
        self._SaveText(fn2)
        return fn
        
    def saveTiff(self,comments='',fn=None,Upload=False):
        if fn==None:
            return None
        t = fn.partition('.')
        fn1 = t[0]+'.tif'
        fn2 = t[0]+'_ps.txt'
        T.imsave(fn1,self.data)
        if comments:
            self.stackparams['Comments']=comments
        self._SaveText(fn2)
        return fn
        
    def _SaveText(self,fn=None):
        if fn==None:
            return False
        s = []
        for parts in self.stackparams:
            s.append('%s : %s \n' % (parts, self.stackparams[parts]))
        s.sort()
        fid = open(fn,'w')
        fid.writelines(s)
        fid.close()
        return True

    def stackTags(self,xx,yy,z1,z2,zs,function='',ps=89):
        '''Date/Time,X,Y,Z1,Z2,Zstep,Exposure(s),CCD Temperature,
        Pixel Size(nm),CCD setting',User
        '''
        self.stackparams.clear()
        self.stackparams['00 function']=function
        self.stackparams['01 Date/Time']=time.asctime()
        self.stackparams['05 X']=xx
        self.stackparams['06 Y']=yy
        self.stackparams['07 Z1']=z1
        self.stackparams['08 Z2']=z2
        self.stackparams['09 Zstep']=zs
        self.stackparams['02 Exposure(s)']=self.ccd.GetExposureTime()
        self.stackparams['03 CCD Temperature']=self.ccd.GetTemperature()
        self.stackparams['10 User']=getuser()
        self.stackparams['11 Coordinates']=self.ccd.GetImageCoordinates()
        self.stackparams['12 Pixel size']=ps
        ccdsett=self.ccd.SettingsString()
        i=13
        for item in ccdsett.splitlines():
            try:
                csi,csv=item.split(':')
                ncsi=str(i)+' '+csi
                self.stackparams[ncsi]=csv
            except:
                csi=0
            i+=1
        return True
            
    def setCoord(self,xs,xe,ys,ye):
        '''val=(xs,xe,ys,ye)'''
        oxs,oxe,oys,oye = self.coordVal
        self.old_Cval = (oxs,oxe,oys,oye)
        self.coordVal=(xs,xe,ys,ye)
        q=self.ccd.SetImageCoordinates(self.coordVal)
        return q
        
    def setCoord2(self,xs,ys,nn):
        '''Gets one corner (xs,ys) to make (nn x nn) image size'''
        oxs,oxe,oys,oye = self.coordVal
        self.old_Cval = (oxs,oxe,oys,oye)
        self.coordVal=(xs,xs+nn-1,ys,ys+nn-1)
        q=self.ccd.SetImageCoordinates(self.coordVal)
        return q
        
    def reSetCoord(self,mode = 0):
        '''mode 0: Full CCD:(1,512,1,512)
           mode 1: return to previous coordinates: old_Cval        
        '''
        if mode==0:
            xs,xe,ys,ye = (1,512,1,512)
        elif mode==1:
            xs,xe,ys,ye = self.old_Cval
        return self.setCoord(xs,xe,ys,ye)
        
