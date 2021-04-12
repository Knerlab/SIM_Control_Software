# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 17:03:40 2020

@author: rl74173
"""

from daq_p35 import daq35 as daq

class attenuate(object):
    def __init__(self):
        self.max = 5.0
        self.voltage = 0.0
        try:
            daq.setVoltage_2(0.0)
        except:
            raise Exception("Data Acquisition Card Not available!")

    def __del__(self):
        daq.setVoltage_2(0.0)

    def setattenuation(self, val):
        if (val<self.max):
            self.voltage = val
            daq.setVoltage_2(self.voltage)
            print('Attenuator set')
        else:
            raise Exception("Voltage too large")