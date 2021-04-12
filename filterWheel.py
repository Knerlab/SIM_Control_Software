# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 21:01:25 2014

@author: Kayvan Forouhesh

@copywrite, Peter Kner, University of Georgia, 2019
"""

import serial
import time

error_code = {'R':'Recieved','D':'Done','O':'Value Out of range','E':'Wheel Timed out did not reach the location specified'}

class filter_wheel():
    def __init__(self,portS='COM5'):
        self.timeout=500;
        self.ser = serial.Serial()
        self.ser.port=portS
        self.ser.close()
        self.ser.baudrate=9600
        self.ser.parity=serial.PARITY_NONE
        self.ser.stopbits=1
        self.ser.bytesize=8
        if(not self.ser.isOpen()):     
            self.ser.open()

        self.res=0
        res=self.ser.portstr       # check which port was really used
        if res==portS:
            try:
#                self.ser.open()
                while not self.ser.isOpen():
                    pass
                print('Filter wheel initialized')
                time.sleep(2)
                print('Moving to position 2')
                err = self.RW('Z',4)
                print(error_code[err])
            except:
                print('error')
        else:
            print('Port Not Connected')
            
    def RW(self,command,length=2):
        '''Enter command without '>'
        List of commands:
            >I  --> Hand shaking, returns  <C (long return mode)
            >A  --> Returns Postion
            >Z  --> Goes to base point (filter 2)
            >Gdn--> Move relative. d is direction (N or P), n is the number of positions.
            >Mn --> Move absolute to position n
            >R  --> Run, position is not updated
            >S  --> Hard Stop
            
        Returns 'R' when a command is recieved,
                'D' when the process is done
                'E' when there is an error
                'O' Out of range
        length is the length of responce requested,
            2 for short - e.g. only 'R'
            4 for long - e.g. 'R' plus 'D' or 'E'
        '''
        if length==2:
            timeout = 20
        else:
            timeout = self.timeout
        self.ser.flushInput()
        command_string = '>%s\n' % command
        self.ser.write(command_string.encode())
        i=0
        while ((self.ser.inWaiting()<length)and(i<timeout)):
            i+=1
            time.sleep(0.01)
        if not i==self.timeout:
            time.sleep(0.1)
            self.res=self.ser.read(self.ser.inWaiting()).decode()
            if length == 2:
                n = self.res.find('<')
                return self.res[n+1]
            elif length == 4:
                n = self.res.find('*')
                return self.res[n+1]
        else:
            print('Timed out')
            return 'E'

    def __del__(self):
        self.ser.close()
        if(self.ser.isOpen()):
            self.ser.close()
        print('Port Closed')

    def getPosition(self):
        val = (self.RW('A',4))
        return int(val)

    def setPositionq(self,pos):
        err = self.RW("M%s" %pos,2)
        return error_code[err]

    def setPositionf(self,pos):
        err = self.RW("M%s" %pos,4)
        return error_code[err]

    def zero(self):
        err = self.RW('Z',4)
        return error_code[err]


