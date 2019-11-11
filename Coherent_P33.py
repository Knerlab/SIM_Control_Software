# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 15:45:04 2016
Class to control Coherent obis lasers

@author: kner
"""

import serial
import time
from serial import SerialException

class obis(object):
    
    def __init__(self,portS='COM7'):
        try:
            self.ser = serial.Serial(
                port=portS,
                baudrate=9600,
                parity=serial.PARITY_NONE,
                stopbits=1,
                bytesize=8
#            stopbits=serial.STOPBITS_ONE,
#            bytesize=serial.EIGHTBITS
            )
            self.res=0
            res=self.ser.portstr       # check which port was really used
            if res==portS:
                print('Port Connected')
            else:
                raise Exception('Port Not Connected')
        except SerialException:
            self.ser.close()
            print('port already open')
        print('Model: ')
        print(self.RW('SYST:INF:MOD?'))
        print(self.RW('SYST:INF:TYP?'))
        print(self.RW('SYST:INF:PNUM?'))
        print('Wavelength: ')
        print(self.RW('SYST:INF:WAV?'))
        print('Operating Hours: ')
        print(self.RW('SYST:HOUR?'))
        print('Maximum Power: ')
        print(self.RW('SOUR:POW:LIM:HIGH?'))
        print(self.IsLaserOn())
        self.maxpower = 0.130
        time.sleep(0.1)
        
    def __del__(self):
        print('closing laser')
        self.SetLaserOff()
        self.ser.close()
        
    def GetPower(self):
        output = self.RW('SOUR:POW:LEV?')
        self.power = float(output)
        return output        
        
    def SetPowerLevel(self,power):
        if power>self.maxpower:
            raise Exception('Max power is %f!' % self.maxpower)
        else:
            self.RW('SOUR:POW:LEV:IMM:AMPL %3.3f' % power)
            return self.err
    
    def IsLaserOn(self):
        state = self.RW('SOUR:AM:STAT?')
        if state == 'OFF':
            self.laseron = False
        else:
            self.laseron = True
        return state
        
    def SetLaserOn(self):
        if not self.laseron:
            self.RW('SOUR:AM:STAT ON')
            self.laseron = True
        return True
    
    def SetLaserOff(self):
        if self.laseron:
            self.RW('SOUR:AM:STAT OFF')
            self.laseron = False
        return True            
        
    def SetCWPowerMode(self):
        self.RW('SOUR:AM:INT CWP')
        return self.err
    
    def SetDigitalModMode(self):
        self.RW('SOUR:AM:EXT DIG')
        return self.err
        
    def QueryLaserMode(self):
        out = self.RW('SOUR:AM:SOUR?')
        return out
        
    def GetLaserStatus(self):
        ''' 'C8001088'
            '80000000' '''
        stat = self.RW('SYST:STAT?')
        fault = self.RW('SYST:FAUL?')
        return (stat, fault)
        
    def RW(self,command):
        self.ser.flushInput()
        ser_command = '%s\r' %command
        self.ser.write(ser_command.encode())
        time.sleep(0.1)
        res=self.ser.read(self.ser.inWaiting()).decode()
        n1 = res.find('<')
        n2 = res.find('\r')
        self.err = res[n2+1:].rstrip('\r\n').lstrip('\r\n')
        return res[n1+1:n2]
    
    