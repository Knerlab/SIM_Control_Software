basepath = 'C:\\Documents and Settings\\kner\\My Documents\\python\\'
import os,sys,time
try:
    sys.path.index(basepath)
except:
    sys.path += [basepath]

#from Priithon.all import F, Y, N, U
import serial

def prior_test():
    ser = serial.Serial(
        port='COM1',
        baudrate=9600,
        bytesize=8,
        parity=serial.PARITY_NONE,
        stopbits=1,
        timeout=2,
        xonxoff=0,
        rtscts=0)
    try:
        ser.open()
        ser.isOpen()
    finally:
        print("closing port")
        ser.close()
    return True

class prior(object):

    priorcodes = {
    1 : 'NO STAGE',
    2 : 'NOT IDLE',
    3 : 'NO DRIVE',
    4 : 'STRING PARSE',
    5 : 'COMMAND NOT FOUND',
    6 : 'INVALID SHUTTER',
    7 : 'NO FOCUS',
    8 : 'VALUE OUT OF RANGE',
    9 : 'INVALID WHEEL',
    10 : 'ARG1 OUT OF RANGE',
    11 : 'ARG2 OUT OF RANGE',
    12 : 'ARG3 OUT OF RANGE',
    13 : 'ARG4 OUT OF RANGE',
    14 : 'ARG5 OUT OF RANGE',
    15 : 'ARG6 OUT OF RANGE',
    16 : 'INCORRECT STATE',
    17 : 'WHEEL NOT FITTED',
    18 : 'QUEUE FULL',
    19 : 'COMPATIBILITY MODE SET',
    20 : 'SHUTTER NOT FITTED',
    21 : 'INVALID CHECKSUM',
    60 : 'ENCODER ERROR',
    61 : 'ENCODER RUN OFF'
    }

    def __init__(self):
        self.ser = serial.Serial(port='COM1',baudrate=9600,bytesize=8,parity=serial.PARITY_NONE,stopbits=1,timeout=2,xonxoff=0,rtscts=0)
        try:
#            self.ser.open()
            print(self.ser.isOpen())
        except:
            print("Prior XY Stage not available!")
            raise 
        # setup some stuff
        self.limits = [-80000,80000,-50000,50000] # soft limits to prevent bad things [xmin,xmax,ymin,ymax]
        #self.limits=[-180000, -125000, -65000, -10000]
        # get current position
        a=self.getPosition()
        print('Current Position is (%d,%d)' % (self.loc[0],self.loc[1]))

    def __del__(self):
        self.ser.close()

    def getPosition(self):
        self.ser.flushInput()
        self.ser.write(b'PS\r')
        time.sleep(0.1)
        t = self.ser.read(self.ser.inWaiting()).decode()
        if self.error(t):
            raise Exception('Command returned an error!')
        self.loc = list(map(float,(t.rstrip('\r')).split(',')))
        return self.loc

    def MoveRelX(self,d=5):
        ''' move by 5 microns '''
        if self.checksoftlimit(d,0):
            print('Move will put you outside soft limits!')
            d = 0
        comstr = 'GR,%d,0,\r' % d
        self.ser.flushInput()
        self.ser.write(comstr.encode())
        time.sleep(0.1)
        t = self.ser.read(self.ser.inWaiting()).decode()
        if self.error(t):
            raise Exception('Command returned an error!')
        return t.strip('\r')

    def MoveRelY(self,d=5):
        ''' move by 5 microns '''
        if self.checksoftlimit(0,d):
            print('Move will put you outside soft limits!')
            d = 0
        comstr = 'GR,0,%d,\r' % d
        self.ser.flushInput()
        self.ser.write(comstr.encode())
        time.sleep(0.1)
        t = self.ser.read(self.ser.inWaiting()).decode()
        if self.error(t):
            raise Exception('Command returned an error!')
        return t.strip('\r')

    def MoveTo(self,x,y):
        ''' Go to point A(x,y) '''
        if self.checksoftlimitAbs(x,y):
            raise Exception('Out of bounds!')
        else:
            cur_pos=self.getPosition()
            dx=x-cur_pos[0]
            dy=y-cur_pos[1]
            comstr = 'GR,%d,%d,\r' % (dx,dy)
            self.ser.flushInput()
            self.ser.write(comstr.encode())
            time.sleep(0.1)
            t = self.ser.read(self.ser.inWaiting()).decode()
            if self.error(t):
                raise Exception('Command returned an error!')
            return t.strip('\r')

    def MoveAbs(self,dx,dy):
        if self.checksoftlimit(dx,dy):
            raise Exception('Out of bounds!')
        else:
            comstr = 'GR,%d,%d,\r' % (dx,dy)
            self.ser.flushInput()
            self.ser.write(comstr.encode())
            time.sleep(0.1)
            t = self.ser.read(self.ser.inWaiting()).decode()
            if self.error(t):
                raise Exception('Command returned an error!')
            return t.strip('\r')

    def checkhardlimits(self):
        # Check hard limits
        self.ser.flushInput()
        self.ser.write(b'LMT\r')
        time.sleep(0.1)
        t = self.ser.read(self.ser.inWaiting()).decode()
        if self.error(t.decode()):
            raise Exception('Command returned an error!')
        x = float('0x'+t.strip('\r'))
        if x==0:
            hardlimit = False
        else:
            hardlimit = True
            hardlimitcode = x
        return hardlimit

    def checksoftlimit(self,dx,dy):
        out = False
        x = self.loc[0] + dx
        y = self.loc[1] + dy
        if (x<self.limits[0]):
            out = True
        if (x>self.limits[1]):
            out = True
        if (y<self.limits[2]):
            out = True
        if (y>self.limits[3]):
            out = True
        return out

    def checksoftlimitAbs(self,x,y):
        out = False
        if (x<self.limits[0]):
            out = True
        if (x>self.limits[1]):
            out = True
        if (y<self.limits[2]):
            out = True
        if (y>self.limits[3]):
            out = True
        return out

    def error(self,op):
        op = (op.rstrip('\r'),)
        if op[0]=='E':
            err = True
            code = int((t.rstrip('\r'))[t.find(',')+1:])
            print(priorcodes[code])
        else:
            err = False
        return err



