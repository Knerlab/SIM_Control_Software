"""
Python code to control the QXGA SLM
@copywrite, Ruizhe Lin and Peter Kner, University of Georgia, 2019
"""

import ctypes as ct
r11 = ct.windll.LoadLibrary('C:/Users/Public/Documents/python_code/QXGA/R11CommLib-1.7-x64.dll')

import subprocess
filepath = r'C:/Program Files/MetroCon-3.3/RepTools/RepSender.exe'

class Dev(ct.Structure):
    pass

Dev._fields_ = [("id", ct.c_char_p), ("next", ct.POINTER(Dev))]

NULL = ct.POINTER(ct.c_int)()
RS485_DEV_TIMEOUT = ct.c_uint16(1000)
RS485_BAUDRATE = ct.c_uint32(256000)
RS232_BAUDRATE = ct.c_uint32(115200)

def repsend(fns):
    output = subprocess.Popen([filepath, '-z', fns, '-d', '0175000547'],stdout=subprocess.PIPE).communicate()[0]
    print(output)
    
def initiate():
    ver = ct.create_string_buffer(8)
    maxlen = ct.c_uint8(10)
    res = r11.R11_LibGetVersion(ver, maxlen)
    if (res==0):
        print('Software version: %s' % ver.value)
    else:
        raise Exception('Libarary not loaded')
    guid = ct.c_char_p(b"54ED7AC9-CC23-4165-BE32-79016BAFB950")
    devcount = ct.c_uint16(0)
    devlist = ct.POINTER(Dev)()
    res = r11.FDD_DevEnumerateWinUSB(guid, ct.pointer(devlist), ct.byref(devcount))
    if (res==0):
        port = devlist.contents.id.decode()
        print('Dev port: %s' % port)
    else:
        raise Exception('Cannot find the port')        
        
def open_usb_port():
    port = ct.c_char_p(b'\\\\?\\usb#vid_19ec&pid_0503#0175000547#{54ed7ac9-cc23-4165-be32-79016bafb950}')
    re = r11.FDD_DevOpenWinUSB(port, RS485_DEV_TIMEOUT)
    if (re==0):
        print('Open Dev port successfully')
        dispTemp = ct.c_uint16(0)
        r11.R11_RpcSysGetDisplayTemp(ct.byref(dispTemp))
        print('Display temperature: %s' %dispTemp.value)
    else:
        raise Exception(' Fail to open the port ')

#boardtype = ct.c_uint8(0)
#res = r11.R11_RpcSysGetBoardType(ct.byref(boardtype))

def getordernum():
    rocount = ct.c_uint16(0)
    res = r11.R11_RpcRoGetCount(ct.byref(rocount))
    if (res==0):
        num = rocount.value
    else:
        raise Exception('Fail to get the order number')
    return(num)

def selecteorder(n):
    roindex = ct.c_uint16(n)
    res = r11.R11_RpcRoSetSelected(roindex)
    if (res==0):
        index = ct.c_uint16(0)
        res = r11.R11_RpcRoGetSelected(ct.byref(index))
        if (res==0)&(index.value==n):
            print('Order is set to #%s' % n)
        else:
            raise Exception('Order set is wrong')
    else:
        raise Exception('Fail to set the order')

def activate():
    res = r11.R11_RpcRoActivate(ct.c_void_p())
    if (res==0):
        print('Activate QXGA successfully')
    else:
        raise Exception('Fail to activate QXGA')

def deactivate():
    res = r11.R11_RpcRoDeactivate(ct.c_void_p())
    if (res==0):
        print('Deactivate QXGA successfully')
    else:
        raise Exception('Fail to deactivate')

def close():
    re = r11.FDD_DevClose()
    if (re==0):
        print('Port closed successfully')
    else:
        raise Exception('Fail to close QXGA')