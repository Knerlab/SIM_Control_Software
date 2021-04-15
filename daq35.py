# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 13:47:11 2013
Code for controlling the NI DAQ board
@copyright, Peter Kner, University of Georgia, 2019

"""
import PyDAQmx as daq
#from PyDAQmx import Task
import PyDAQmx.DAQmxConstants as cnst
import PyDAQmx.DAQmxTypes as tp
import numpy as N
import ctypes as ct

def CCDTrig_open(freq, data):
    taskHandle = tp.TaskHandle()
    taskCounter = tp.TaskHandle()
    nsamples = len(data)
    idata = data.astype(tp.uInt32)
    startdata = N.array((5,5),dtype=tp.uInt32)
    try:
        daq.DAQmxCreateTask("",taskCounter)
        daq.DAQmxCreateCOPulseChanFreq(taskCounter,"Dev1/ctr0","",daq.DAQmx_Val_Hz,daq.DAQmx_Val_Low,0.0,freq,0.50)
        daq.DAQmxCfgImplicitTiming(taskCounter,daq.DAQmx_Val_ContSamps,1000)
#        print "taskCounter Value", taskCounter.value
    ##  Create Digital Channel
        daq.DAQmxCreateTask("",taskHandle)
#        print "taskHandle Value", taskHandle.value
        daq.DAQmxCreateDOChan(taskHandle,"Dev1/port0/line0:7","",daq.DAQmx_Val_ChanForAllLines)
        daq.DAQmxCfgSampClkTiming(taskHandle,"Ctr0InternalOutput",freq,daq.DAQmx_Val_Rising,daq.DAQmx_Val_FiniteSamps,nsamples)
    ##  Register done events (error checking)
#        q = daq.DAQmxRegisterDoneEvent(taskHandle,0,None,None)
#        print q,'c'
#        q = daq.DAQmxRegisterDoneEvent(taskCounter,0,None,None)
#        print q,'d'
    ##  write startup values
        daq.DAQmxWriteDigitalU32(taskHandle,2,0,10.0,cnst.DAQmx_Val_GroupByChannel,startdata,None,None)
        daq.DAQmxStartTask(taskHandle)
        daq.DAQmxStopTask(taskHandle)
    ##  load trigger sequence
        daq.DAQmxWriteDigitalU32(taskHandle,nsamples,0,10.0,cnst.DAQmx_Val_GroupByChannel,idata,None,None)
        print("Data Written\n")
    ##  start counter
        daq.DAQmxStartTask(taskCounter)
        return taskHandle,taskCounter
    except:
        errBuff=tp.create_string_buffer(b"",2048)
        daq.DAQmxGetExtendedErrorInfo(errBuff,2048)
        print(errBuff.value)

def CCDTrig_run(taskHandle,taskCounter):
    try:
        daq.DAQmxStartTask(taskHandle)
        #print("Tasks Started\n")
        daq.DAQmxWaitUntilTaskDone(taskHandle,10)
        daq.DAQmxStopTask(taskHandle)
        return 0
    except:
        errBuff=tp.create_string_buffer(b"",2048)
        daq.DAQmxGetExtendedErrorInfo(errBuff,2048)
        print(errBuff.value)
        
def CCDTrig_close(taskHandle,taskCounter):
    try:
        if (taskHandle != 0 ):
            #print("Stopping Tasks\n")
            daq.DAQmxStopTask(taskHandle)
            daq.DAQmxClearTask(taskHandle)
            daq.DAQmxStopTask(taskCounter)
            daq.DAQmxClearTask(taskCounter)
            return 0
    except:
        errBuff=tp.create_string_buffer(b"",2048)
        daq.DAQmxGetExtendedErrorInfo(errBuff,2048)
        print(errBuff.value)

def DAQreset():
    try:
        q = daq.DAQmxResetDevice("Dev1")
        return q
    except:
        errBuff=tp.create_string_buffer(b"",2048)
        daq.DAQmxGetExtendedErrorInfo(errBuff,2048)
        print(errBuff.value)
        
def getVoltage():
    try:
        taskHandle = tp.TaskHandle()
        daq.DAQmxCreateTask("",taskHandle)
#        print "taskHandle Value", taskHandle.value
        val = tp.float64()
        daq.DAQmxCreateAIVoltageChan(taskHandle,"Dev1/ai0","",cnst.DAQmx_Val_RSE,0.0,10.0,cnst.DAQmx_Val_Volts,"")
        daq.DAQmxStartTask(taskHandle)
        daq.DAQmxReadAnalogScalarF64(taskHandle,5.0,daq.byref(val),None)
        if not taskHandle == 0 :
#            print "Stopping Tasks\n"
            daq.DAQmxStopTask(taskHandle)
            daq.DAQmxClearTask(taskHandle)
        return val.value
    except:
        errBuff=tp.create_string_buffer(b"",2048)
        daq.DAQmxGetExtendedErrorInfo(errBuff,2048)
        print(errBuff.value)


def setVoltage(val):
    try:
        taskHandle = tp.TaskHandle()
        daq.DAQmxCreateTask("",taskHandle)
#        print "taskHandle Value", taskHandle.value
        daq.DAQmxCreateAOVoltageChan(taskHandle,"Dev1/ao0","",0.0,10.0,cnst.DAQmx_Val_Volts,"")
        daq.DAQmxStartTask(taskHandle)
        daq.DAQmxWriteAnalogScalarF64(taskHandle,1,5.0,tp.float64(val),None)
        if not taskHandle == 0 :
#            print "Stopping Tasks\n"
            daq.DAQmxStopTask(taskHandle)
            daq.DAQmxClearTask(taskHandle)
        return 0
    except:
        errBuff=tp.create_string_buffer(b"",2048)
        daq.DAQmxGetExtendedErrorInfo(errBuff,2048)
        print(errBuff.value)

def setVoltage_2(val):
    try:
        taskHandle = tp.TaskHandle()
        daq.DAQmxCreateTask("",taskHandle)
#        print "taskHandle Value", taskHandle.value
        daq.DAQmxCreateAOVoltageChan(taskHandle,"Dev1/ao1","",0.0,10.0,cnst.DAQmx_Val_Volts,"")
        daq.DAQmxStartTask(taskHandle)
        daq.DAQmxWriteAnalogScalarF64(taskHandle,1,5.0,tp.float64(val),None)
        if not taskHandle == 0 :
#            print "Stopping Tasks\n"
            daq.DAQmxStopTask(taskHandle)
            daq.DAQmxClearTask(taskHandle)
        return 0
    except:
        errBuff=tp.create_string_buffer(b"",2048)
        daq.DAQmxGetExtendedErrorInfo(errBuff,2048)
        print(errBuff.value)
