# -*- coding: utf-8 -*-

import time
from daq_p35 import daq35 as daq

class zstage(object):
    def __init__(self):
        try:
            daq.setVoltage(0.0)
        except:
            raise Exception("Data Acquisition Card Not available!")
        self.cal = (10.0/200.0) # units microns
        self.tol = 0.05
        self.voltage = 0.0
        self.pos = self.getPosition()

    def __del__(self):
        daq.setVoltage(0.0)

    def getPosition(self):
        v = daq.getVoltage()
        self.pos = v/self.cal
        return (self.pos)

    def setPositionq(self,pos):
        ''' doesn't wait for piezo to move '''
        self.voltage = self.cal*pos
        daq.setVoltage(self.voltage)
        return self.getPosition()

    def setPositionf(self,pos):
        ''' waits for piezo to reach position '''
        self.voltage = self.cal*pos
        daq.setVoltage(self.voltage)
        diff = abs(pos-self.getPosition())
        for m in range(250):
            time.sleep(0.001)
            diff = abs(pos-self.getPosition())
            if (diff<self.tol):
                break
        return self.getPosition()
