#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      bent83
#
# Created:     08/10/2013
# Copyright:   (c) bent83 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import sys
import time
import ctypes
fl = 'C:\\Program FIles (x86)\\MicroSMC\\'
sys.path.append(fl)
from pol_rotator_p35.Mot_pol_params import *

byref = ctypes.byref

mot = ctypes.cdll.LoadLibrary('USMCDLL.dll')
#runfile(r'Z:\Mot_pol_params.py', wdir=r'Z:')

class A8MRU(object):

    def __init__(self):
        self.device = 0
        self.devs = USMC_Devices()
        out = mot.USMC_Init(byref(self.devs))
        self.mode,self.params,self.stparams = self.initializeSettings()
        self.setMode()
        self.setParams()
        self.state = USMC_State()
        self.degperstep = 0.075
        self.speed = ctypes.c_float(100)

    def __del__(self):
        self.MoveAbsWait(0)
        self.MotorOff()
        out = self.close()
        print(out)

    def setSpeed(self,vel=100,step=4):
        ''' step = 1,2,4,8
            vel in degrees / sec '''
        if vel > 350:
            step = 1
        if step in (1,2,4,8):
            self.stparams.SDivisor = step
            speedsteps = vel/(step*self.degperstep)
            self.speed = ctypes.c_float(speedsteps)
        out = self.setParams()
        return out

    def MotorOn(self):
        self.mode.ResetD = 0
        out = self.setMode()
        return out

    def MotorOff(self):
        self.mode.ResetD = 1
        out = self.setMode()
        return out

    def MoveAbs(self,pos):
        ''' move to absolute position
            pos in degrees '''
        index = int(pos/self.degperstep)
        out = mot.USMC_Start(self.device,index,byref(self.speed),byref(self.stparams))
        return out

    def MoveAbsWait(self,pos):
        starttime = time.time()
        index = int(pos/self.degperstep)
        out = mot.USMC_Start(self.device,index,byref(self.speed),byref(self.stparams))
        while (abs(pos-self.getCurPos())>2):
            time.sleep(0.01)
        endtime = time.time()
        print (endtime-starttime)
        return self.getCurPos()

    def MoveRel(self,jog=1):
        self.getState()
        index = self.state.CurPos + int(jog/self.degperstep)
        out = mot.USMC_Start(self.device,index,byref(self.speed),byref(self.stparams))
        return out

    def init(self):
        devs = USMC_Devices()
        out1 = mot.USMC_Init(byref(devs))
        return (devs,out1)

    def close(self):
        out = mot.USMC_Close()
        return out

    def getState(self):
        self.state = USMC_State()
        out = mot.USMC_GetState(self.device,byref(self.state))
        return out

    def getCurPos(self):
        self.state = USMC_State()
        out = mot.USMC_GetState(self.device,byref(self.state))
        pos = self.state.CurPos*(360./4800.)
        return pos

    def getMode(self):
        out = mot.USMC_GetMode(self.device,byref(self.mode))
        return out

    def setMode(self):
        out = mot.USMC_SetMode(self.device,byref(self.mode))
        return out

    def getParams(self):
        out = mot.USMC_GetParameters(self.device,byref(self.params))
        return out

    def setParams(self):
        out = mot.USMC_SetParameters(self.device,byref(self.params))
        return out

    def getStartParams(self):
        out = mot.USMC_GetStartParameters(self.device,byref(self.stparams))
        return out

    def Start(pos,stparams):
        device = 0
        speed = ctypes.c_float(100)
        out = mot.USMC_Start(device,pos,byref(speed),byref(stparams))
        return (speed,stparams,out)

    def Stop(self):
        out = mot.USMC_Stop(self.device)
        return out

    def initializeSettings(self):
        ''' mode '''
        mode = USMC_Mode()
        mode.PMode = 1 # Turn off buttons (1 = buttons disabled)
        mode.PReg = 1 # current reduction regime
        mode.ResetD = 1 # turn power off and make a whole step (True = apply)
        mode.EMReset = 0 # Quick power off
        mode.Tr1T = 0 # limit switch 1 True state
        mode.Tr2T = 0 # limit switch 2 True state
        mode.RotTrT = 0 # Rotary Transducer True state
        mode.TrSwap = 0 # if True, limit switches are to be swapped
        mode.Tr1En = 0 # if True, limit switch 1 enabled
        mode.Tr2En = 0 # if True, limit switch 2 enabled
        mode.RotTeEn = 0 # if True, rotary Transducer Operation enabled
        mode.RotTrOp = 0 # Rotary Transducer Operation Select (stop on error if True)
        mode.Butt1T = 0 # Button 1 True state
        mode.Butt2T = 0 # Button 2 True state
        mode.ResetRT = 0 # Reset Rotary Transducer
        mode.SyncOutEn = 0 # if True output synchornization enabled
        mode.SyncOUTR = 0 # if True output synchronization counter will be reset
        mode.SyncINOp = 0 # Synchronization input mode
        mode.SyncCount = 0 # number of steps after which synchronization output signal occurs
        mode.SyncInvert = 0 # Set this bit to True to invert output synchornization polarity
        mode.EncoderEn = 0 # Enable Encoder on pins (syncin,rottr)
        mode.EncoderInv = 0 # Invert Encoder Counter Direction
        mode.ResBEnc = 0 # Reset Encoder
        mode.ResEnc = 0 # Reset Encoder
        mode.Reserved # not used
        ''' params '''
        params = USMC_Parameters()
        params.AccelT = 980 # acceleration time in ms
        params.DecelT = 980 # deceleration time in ms
        params.PTimeout = 100 # Time (in ms) after which current will be reduced to 60% of normal
        params.BTimeout1 = 1000 # Time (in ms) after which speed of motor rotation will be equal to the one specified in BT01P
        params.BTimeout2 = 1000
        params.BTimeout3 = 1000
        params.BTimeout4 = 1000
        params.BTimeoutR # Time (in ms) after which reset command will be performed
        params.BTimeoutD # This field reserved
        params.MinP # Speed (steps/sec) while performing reset
        params.BTO1P = 2 # Speed after btimeout1
        params.BTO2P = 8
        params.BTO3P = 32
        params.BTO4P = 128
        params.MaxLoft = 1024 # Value in full steps that will be used in backlash operation
        params.StartPos = 1 # Current position saved to FLASH
        params.RTDelta = 600 # Revolution distance -- number of full steps per one full revolution
        params.RTMinError = 2 # number of full steps missed to raise error flag
        params.MaxTemp = 50 # maximum allowed temp
        params.SynOutP = 5 # duration of output synchronization pulse
        params.LoftPeriod = 3 # Speed (steps/sec) of the lst phase of the backlash operation
        params.EncMult = 1.0 # Encoder step multiplier
        params.Reserved = 0 # NA
        ''' start parameters '''
        stparams = USMC_StartParameters()
        stparams.SDivisor = 4 # Step is divided by this factor (1,2,4,8)
        stparams.DefDir = 0 # Direction for backlash operation (relative)
        stparams.LoftEn = 0 # Enable automatic backlash operation (works if slow start/stop mode is off)
        stparams.SlStart = 1 # if True slow start/stop enabled
        stparams.WSyncIN = 0 # if True, controller will wait for input synchronization to start
        stparams.SyncOutR = 0 # If True, output synchronization counter will be reset
        stparams.ForceLoft = 0 # if True and destination position is equal to the current position, backlash will be performed
        return (mode,params,stparams)
