# -*- coding: utf-8 -*-
"""
AndorCCD cameare control
@copywrite, Peter Kner, University of Georgia, 2019
"""

import numpy as N
import ctypes as ct
#import fftfuncs as F
atm64 = ct.windll.LoadLibrary(r'C:/Users/Public/Documents/python_code/SIM_AO_p36/pyAndor_p35/atmcd64d.dll')
from pyAndor_p35 import pyAndor_codes
import time

class ccd(object):
    "a class for dealing with the camera"

    andorfilelocation = ct.create_string_buffer(b'C:/Program Files (x86)/Andor iXon')
    temperature = 19
    cooler = False         # On
    coolermode = 0 # 0: return to ambient on exit, 1: stay cool
    exposure_time = 0.1
    acquisition_mode = 1
    acquisition_mode_dict = {1: 'Single Scan', 2:'Accumulate', 3:'Kinetics', 4:'Fast Kinetics', 5:'Run until Abort', 6:'Frame Transfer',
                            7:'Run until Abort FT', 9:'Time Delayed Integration'}
    read_mode = 4         #Image
    read_mode_dict = {0:'Full Vertical Binning', 1:'Multi-Track', 2:'Random-Track', 3:'Single-Track', 4:'Image'}
    trigger_mode = 0
    trigger_mode_dict = {0:'Internal', 1:'External', 6:'External Start'}
    fast_ext_trig = 0 # 0: disabled
    ad_channel = 0        # 14-bit -- no 16-bit on this camera
    ad_channel_dict = {0:'14-bit',1:'16-bit'}
    output_amplifier = 1  # Conventional
    output_amplifier_dict = {0:'EM',1:'Conventional'}
    emccdgain = 0
    preampgain = 1 # preamp gain was 1x prior to 4/2/2008
    preampgain_dict = {0:'1x',1:'2.4x',2:'4.6x'}
    shutter_type = 1
    shutter_type_dict = {0:'Shutter Active Low', 1:'Shutter Active High'}
    shutter_mode = 0
    shutter_mode_dict = {0:'Auto', 1:'Open', 2:'Closed'}
    shutter_closing_time = 10
    shutter_opening_time = 10
    VS_speed = 5          # 0.6us
    VS_speed_dict = {0:'0.4us', 1:'0.6us',2:'1us',3:'1.8us',4:'3.4us',5:'6.6us',6:'13us'}
    HS_speed = 0          # 1MHz for conv., 10MHz for EM
    VSvoltage = 0
    VSvoltage_dict = {0:'Normal',1:'+1',2:'+2',3:'+3',4:'+4'}
    set_frame_transfer_mode = 0
    image_coor = (1,512,1,512)
    image_size = (512,512)
    image_bin = (1,1)
    kinetic_length = 1
    kinetic_time = 1
    accumulations_number = 1
    acc_time = 1
    images = N.zeros((1,512,512), dtype = N.uint16)
    images = N.ascontiguousarray(images,dtype=N.uint16)
    image_stack = []
    zpos = 0

    def __init__(self):
        print('Initializing Andor camera.')
        self.andor_dict = dict([(v, k)    for (k, v) in pyAndor_codes.__dict__.items() if k.startswith('DRV_')])
        q = atm64.Initialize(self.andorfilelocation)
        print(self.andor_dict[q])
        print(q)
        self.CoolerON()
        self.SetCoolerMode(self.coolermode)
        self.SetDMAParameters(1,0.0)
        self.SetImageCoordinates(self.image_coor)
        self.SetAcquisitionMode(self.acquisition_mode)#
        self.SetReadMode(self.read_mode)
        self.SetExposureTime(self.exposure_time)
        self.SetTriggerMode(self.trigger_mode)
        self.SetFastExtTrigger(self.fast_ext_trig)
        self.SetShutterMode(self.shutter_mode)
        self.SetADChannel(self.ad_channel)
        self.SetOutputAmplifier(self.output_amplifier)
        self.SetHSSpeed(self.HS_speed)
        self.SetVSSpeed(self.VS_speed)
        self.SetEMCCDGain(self.emccdgain)
        self.SetPreAmpGain(self.preampgain)
        self.SetFrameTransferMode(self.set_frame_transfer_mode)

    def __del__(self):#
        t = self.GetTemperature()
        print('CCD temperature is ',t)
        if (t<0 and self.coolermode==0):
            self.SetTemperature(15)
            while (t<0.0):
                t = self.GetTemperature()
                print('CCD temperature is ',t)
                time.sleep(1)
        print('Shutting down CCD camera')
        atm64.ShutDown()

    def Acquire(self):#
        q = atm64.StartAcquisition()
        #while not (self.andor_dict[atm32.GetStatus()[1]] == 'DRV_IDLE'):
        #    wx.Usleep(10)
        #print self.andor_dict[q]
        return self.andor_dict[q]

    def WaitForFinish(self):
        # add timeout?
        q = self.andor_dict[self.GetStatus()]
        while not (q == 'DRV_IDLE'):
            # wx.Usleep(10)
            q = self.andor_dict[self.GetStatus()]
        return q

    def WaitForNewData_old(self):#
        ''' with 5s timeout '''
        i = 0
        q = self.GetNewImage()
        while (q=='DRV_NO_NEW_DATA'):
            time.sleep(0.01) # was 0.01 4/12/2014
            q = self.GetNewImage()
            i += 1
            if (i > 500): # same
                print('Aborting Acquisition')
                self.AbortAcquisition()
                self.SetShutterMode(2)
                test = 'Camera frozen waiting for data' 
                raise Exception(test)
                break
        return q
    
    def WaitForNewData(self):
        ''' try built in wait function 04/23/2014 '''
        q = atm64.WaitForAcquisitionTimeOut(ct.c_int(1000))
        q = self.GetNewImage()        
        return q

    def GetStatus(self):#
        stat = ct.pointer(ct.c_int(0))
        q = atm64.GetStatus(stat)
        print(self.andor_dict[q])
        return stat.contents.value

    def AbortAcquisition(self):#
        q = atm64.AbortAcquisition()
        return self.andor_dict[q]

    def GetAcquisitionTimings(self):#
        exposure   = ct.pointer(ct.c_float(0))
        accumulate = ct.pointer(ct.c_float(0))
        kinect     = ct.pointer(ct.c_float(0))
        q = atm64.GetAcquisitionTimings(exposure,accumulate,kinect)
        #print 'Exposure time is ',q[1]
        #print 'Accumulate cycle time is ',q[2]
        #print 'Kinetic cycle time is ',q[3]
        return (exposure.contents.value,accumulate.contents.value,kinect.contents.value)

    def GetNewImage(self):#
        image = N.zeros(self.image_size,dtype=N.uint16)
#        imager = N.ascontiguousarray(self.images)
#        q = atm32.GetNewData16(imager.ctypes.data_as(ct.POINTER(ct.c_uint8)),ct.c_ulong(imager.size))
        q = atm64.GetNewData16(image.ctypes.data_as(ct.POINTER(ct.c_uint8)),ct.c_ulong(image.size))
#        self.images = imager.byteswap()  # it was imager[0]
        self.images[:,:]=image
        return self.andor_dict[q]

    def GetImages(self):
        fImg = ct.pointer(ct.c_long(0))
        lImg = ct.pointer(ct.c_long(0))
        q = atm64.GetNumberNewImages(fImg,lImg)
        nni=(fImg.contents.value,lImg.contents.value)
        print('number of images ', nni)
        fVal = ct.pointer(ct.c_long(0))
        lVal = ct.pointer(ct.c_long(0))
        imager = N.ascontiguousarray(self.images)
        q = atm64.GetImages16(ct.c_long(nni[1]),ct.c_long(nni[2]),imager.ctypes.data_as(ct.POINTER(ct.c_uint8)),ct.c_ulong(imager.size),fVal,lVal)
        print(self.andor_dict[q])
        return fVal.contents.value,lVal.contents.value

    def InitializeImages(self):#
        p = (self.kinetic_length,self.image_size[0],self.image_size[1])
        # self.images = F.zeroArrU(p)
        return p

    def InitializeStack(self,layers):#
        p = (layers,self.image_size[0],self.image_size[1])
        self.image_stack = []
        # self.image_stack = F.zeroArrU(p)
        self.zpos = 0
        return p

    def AddImagetoStack(self):#
        self.image_stack[self.zpos,:,:] = self.images[0,:,:]
        self.zpos += 1

    def SetCoolerMode(self,n):#
        ''' 0: temperature returns to ambient on shut down
            1: temperature is maintained on shut down '''
        self.coolermode = n
        q = atm64.SetCoolerMode(ct.c_int(n))
        return self.andor_dict[q]

    def GetCoolerMode(self):#
        return self.coolermode

    def CoolerON(self):#
        self.cooler = True
        q = atm64.CoolerON()
        return self.andor_dict[q]

    def CoolerOFF(self):#
        self.cooler = False
        q = atm64.CoolerOFF()
        return self.andor_dict[q]

    def CoolerStatus(self):#
        return self.cooler

    def SetTemperature(self,T):#
        q = atm64.SetTemperature(ct.c_int(T))
        return self.andor_dict[q]

    def GetTemperature(self):#
        temp = ct.c_int(0)
        q = atm64.GetTemperature(ct.byref(temp))
        print(self.andor_dict[q])
        return temp.value

    def GetTemperatureStatus(self):#
        temp = ct.pointer(ct.c_int(0))
        q = atm64.GetTemperature(temp)
        print(self.andor_dict[q])
        val = (self.andor_dict[q] == 'DRV_TEMP_STABILIZED')
        return val

    def SetFanHigh(self):#
        q = atm64.SetFanMode(ct.c_int(0))
        return self.andor_dict[q]

    def SetFanLow(self):#
        q = atm64.SetFanMode(ct.c_int(1))
        return self.andor_dict[q]

    def SetFanOff(self):#
        q = atm64.SetFanMode(ct.c_int(2))
        return self.andor_dict[q]

    def SetExposureTime(self,val):#
        self.exposure_time = val
        q = atm64.SetExposureTime(ct.c_float(val))
        return self.andor_dict[q]

    def GetExposureTime(self):#
        return self.exposure_time

    def SetAcquisitionMode(self,mode):#
        if (mode in self.acquisition_mode_dict):
            self.acquisition_mode = mode
            atm64.SetAcquisitionMode(ct.c_int(mode))
            return self.acquisition_mode_dict[mode]
        else:
            return 'Error'

    def GetAcquisitionMode(self):#
        return self.acquisition_mode

    def SetReadMode(self,mode):#
        if (mode in self.read_mode_dict):
            self.read_mode = mode
            atm64.SetReadMode(ct.c_int(mode))
            return self.read_mode_dict[mode]
        else:
            return 'Error'

    def GetReadMode(self):#
        return self.read_mode

    def SetTriggerMode(self,mode):#
        if (mode in self.trigger_mode_dict):
            self.trigger_mode = mode
            atm64.SetTriggerMode(ct.c_int(mode))
            return self.trigger_mode_dict[mode]
        else:
            return 'Error'

    def GetTriggerMode(self):#
        return self.trigger_mode

    def SetFastExtTrigger(self,val):#
        val = (val>0)
        self.fast_ext_trig = val
        q = atm64.SetFastExtTrigger(ct.c_int(val))
        return self.andor_dict[q]

    def GetFastExtTrigger(self):
        return self.fast_ext_trig

    def SetADChannel(self,ch):
        if (ch in self.ad_channel_dict):
            self.ad_channel = ch
            atm64.SetADChannel(ct.c_int(ch))
            return self.ad_channel_dict[ch]
        else:
            return 'Error'

    def GetADChannel(self):#
        return self.ad_channel

    def SetOutputAmplifier(self,val):#
        if (val in self.output_amplifier_dict):
            self.output_amplifier = val
            atm64.SetOutputAmplifier(ct.c_int(val))
            return self.output_amplifier_dict[val]
        else:
            return 'Error'

    def GetOutputAmplifier(self):#
        return self.output_amplifier

    def SetEMCCDGain(self,val):#
        val = int(val)
        if ((val>=0) and (val<=4095)):
            self.emccdgain = val
            q = atm64.SetEMCCDGain(ct.c_int(val))
            return self.andor_dict[q]
        else:
            return 'Error'

    def GetEMCCDGain(self):#
        return self.emccdgain

    def SetPreAmpGain(self,val):#
            val = int(val)
            if (val in self.preampgain_dict):
                self.preampgain = val
                q = atm64.SetPreAmpGain(ct.c_int(val))
                return self.andor_dict[q]
            else:
                return 'Error'

    def GetPreAmpGain(self):#
        return self.preampgain

    def SetShutterType(self,val):#
        val = int(val)
        if (val in self.shutter_type_dict):
            self.shutter_type = val
            q = atm64.SetShutter(ct.c_int(self.shutter_type),ct.c_int(self.shutter_mode),ct.c_int(self.shutter_closing_time),ct.c_int(self.shutter_opening_time))
            return self.andor_dict[q]
        else:
            return 'Error'

    def SetShutterMode(self,val):#
        val = int(val)
        if (val in self.shutter_mode_dict):
            self.shutter_mode = val
            q = atm64.SetShutter(ct.c_int(self.shutter_type),ct.c_int(self.shutter_mode),ct.c_int(self.shutter_closing_time),ct.c_int(self.shutter_opening_time))
            return self.andor_dict[q]
        else:
            return 'Error'

    def SetShutterOpeningTime(self,val):#
        if (val>0):
            self.shutter_opening_time = val
            q = atm64.SetShutter(ct.c_int(self.shutter_type),ct.c_int(self.shutter_mode),ct.c_int(self.shutter_closing_time),ct.c_int(self.shutter_opening_time))
            return self.andor_dict[q]
        else:
            return 'Error'

    def SetShutterClosingTime(self,val):#
        if (val>0):
            self.shutter_opening_time = val
            q = atm64.SetShutter(ct.c_int(self.shutter_type),ct.c_int(self.shutter_mode),ct.c_int(self.shutter_closing_time),ct.c_int(self.shutter_opening_time))
            return self.andor_dict[q]
        else:
            return 'Error'

    def GetShutterState(self):#
        q = [self.shutter_type_dict[self.shutter_type],self.shutter_mode_dict[self.shutter_mode],self.shutter_closing_time,self.shutter_opening_time]
        return q

    def SetVSAmplitude(self,val):#
        val = int(val)
        if (val in self.VSvoltage_dict):
            self.VSvoltage = val
            q = atm64.SetVSAmplitude(ct.c_int(val))
            return self.andor_dict[q]
        else:
            return 'Error'

    def GetVSAmplitude(self):#
        return self.VSvoltage

    def SetVSSpeed(self,val):#
        val = int(val)
        if (val in self.VS_speed_dict):
            self.VS_speed = val
            q = atm64.SetVSSpeed(ct.c_int(val))
            return self.andor_dict[q]
        else:
            return 'Error'

    def GetVSSpeed(self):#
        return self.VS_speed

    def SetHSSpeed(self,val):#
        val = int(val)
        self.HS_speed = val
        q = atm64.SetHSSpeed(ct.c_int(self.output_amplifier),ct.c_int(val))
        return self.andor_dict[q]

    def GetHSSpeed(self):#
        return self.HS_speed

    def SetFrameTransferMode(self,val):#
        val = (val>0)
        self.set_frame_transfer_mode = val
        q = atm64.SetFrameTransferMode(ct.c_int(val))
        return self.andor_dict[q]

    def GetFrameTransferMode(self):#
        return self.set_frame_transfer_mode

    def SetImageCoordinates(self,val):#
        " val = (xs,xe,ys,ye)"
        self.image_coor = val
        xs,xe,ys,ye = val
        hbin,vbin = self.image_bin
        xsize = int((xe-xs+1)/hbin)
        ysize = int((ye-ys+1)/vbin)
        self.image_size = (xsize,ysize)
        self.images = N.zeros((1,xsize,ysize), dtype=N.uint16)
        chbin = ct.c_int(hbin)
        cvbin = ct.c_int(vbin)
        cxs = ct.c_int(xs)
        cxe = ct.c_int(xe)
        cys = ct.c_int(ys)
        cye = ct.c_int(ye)
        q = atm64.SetImage(chbin,cvbin,cxs,cxe,cys,cye)
        return self.andor_dict[q]

    def GetImageCoordinates(self):#
        return(self.image_coor)

    def SetBinning(self,hbin,vbin):#
        self.image_bin = (hbin,vbin)
        xs,xe,ys,ye = self.image_coor
        xsize = int((xe-xs+1)/hbin)
        ysize = int((ye-ys+1)/vbin)
        self.image_size = (xsize,ysize)
        # self.images = F.zeroArrU((1,xsize,ysize))
        chbin = ct.c_int(hbin)
        cvbin = ct.c_int(vbin)
        cxs = ct.c_int(xs)
        cxe = ct.c_int(xe)
        cys = ct.c_int(ys)
        cye = ct.c_int(ye)
        q = atm64.SetImage(chbin,cvbin,cxs,cxe,cys,cye)
        return self.andor_dict[q]

    def GetBinning(self):#
        return(self.image_bin)

    def GetImageSize(self):#
        return(self.image_size)

    def GetImageSizeYX(self):#
        x,y = self.image_size
        return((y,x))

    def SetKineticLength(self,val):#
        val = int(val)
        self.kinetic_length = val
        q = atm64.SetNumberKinetics(ct.c_int(val))
        return self.andor_dict[q]

    def GetKineticLength(self):#
        return(self.kinetic_length)

    def SetKineticTime(self,val):#
        self.kinetic_time = val
        q = atm64.SetKineticCycleTime(ct.c_int(val))
        return self.andor_dict[q]

    def GetKineticTime(self):#
        return(self.kinetic_time)

    def SetAccumulations(self,val):#
        val = int(val)
        self.accumulations_number = val
        q = atm64.SetNumberAccumulations(ct.c_int(val))
        return self.andor_dict[q]

    def GetAccumulations(self):#
        return(self.accumulations_number)

    def SetAccTime(self,val):#
        self.acc_time = val
        q = atm64.SetAccumulationCycleTime(ct.c_float(val))
        return self.andor_dict[q]

    def GetAccTime(self):#
        return(self.acc_time)

    def SetDMAParameters(self,MaxImagesPerDMA=1,TimePerDMA=0.0):#
        q = atm64.SetDMAParameters(ct.c_int(MaxImagesPerDMA),ct.c_float(TimePerDMA))
        return self.andor_dict[q]
        
    def ManualShutdown(self):#
        t = self.GetTemperature()
        print('CCD temperature is ',t)
        if (t<0):
            self.SetTemperature(20)            
            while (t<15):
                t = self.GetTemperature()
                print('CCD temperature is ',t)
                time.sleep(1)
        print('Shutting down CCD camera')
        atm64.ShutDown()

    def SettingsString(self):#
        q = (self.emccdgain,self.ad_channel_dict[self.ad_channel],\
        self.output_amplifier_dict[self.output_amplifier],\
        self.preampgain_dict[self.preampgain],self.VS_speed_dict[self.VS_speed],\
        self.VSvoltage_dict[self.VSvoltage])
#        string = '''CCD Settings: EMCCDGain %d\nAD Channel %s\nOutput Amplifier %s\npreamp gain %s\nVertical Shift Speed %s\nVertical Shift Voltage %s\n''' % q
        string = '''\nEMCCDGain: %d\nAD Channel: %s\nOutput Amplifier: %s\npreamp gain: %s\nVertical Shift Speed: %s\nVertical Shift Voltage: %s\n''' % q
        return string

    def GetReadoutSpeed(self):#
        v = ct.pointer(ct.c_float(0))
        q = atm64.GetHSSpeed(ct.c_int(0),ct.c_int(self.output_amplifier),ct.c_int(self.HS_speed),)
        print(self.andor_dict[q])
        return v.contents.value

        
