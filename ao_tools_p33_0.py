# -*- coding: utf-8 -*-
"""
Code of Adaptive Optics Implementation
@copywrite, Ruizhe Lin and Peter Kner, University of Georgia, 2019
"""
import time, os
import numpy as N
import Utility as U
import pylab
import zernike
from DAQ27 import daq
import pr_batch_P33_0 as prb
import Phase_Retrieval_b as phaseret
import tifffile as tf
import qxga_exec_p27 as slm
import scipy.ndimage.filters as filters
import metric_tools_p27 as mt
import csv

FLAT_FILE = r'C:/Users/Public/Documents/python_code/SIM_P27/20190408140547_flatfile.txt'

fns = r'C:/Users/Public/Documents/python_code/QXGA/repzfiles/32pix_ao.repz11'

infl_fcn_file = r'C:/Users/Public/Documents/python_code/SIM_P27/influ_func_04012019.tif'
infl_fcn = tf.imread(infl_fcn_file)
Sm = prb.Smat(infl_fcn)

with open(FLAT_FILE, "r") as file:
    cmd_best = eval(file.readline())
#cmd_best = [0.] * 69
mod = N.arange(13)
zmv = N.zeros((13))

results = []

def scan_dm(om,dm,push=0.1,actlist=None,path=None,center=None):
    if actlist==None:
        actlist = range(69)
    if path==None:
        path = om.path
    if center==None:
        center = om.zst.getPosition()
#    cmd_flat = dm.cmd #dm.read_command_file(dm.flat_file)
    cmd_arr = [0.] * 69
#    # flat
    dm.Send(cmd_arr)
    time.sleep(0.1)
    om.StackExt(center-3.0,center+3.0,0.2,verbose=False)
    filename = path + '\\stack_flat_start.tif'
    om.saveTiff(fn=filename)
    for m in actlist:
        print("actuator " + str(m))
        cmd_arr[m] = push
        dm.Send( cmd_arr )
        time.sleep(0.1)
        om.StackExt(center-3.0,center+3.0,0.2,verbose=False)
        filename = path + '\\stack_act%02d.tif' % m
        om.saveTiff(fn=filename)
        cmd_arr[m] = 0.0
    # flat
    dm.Send(cmd_arr)
    time.sleep(0.1)
    om.StackExt(center-3.0,center+3.0,0.2,verbose=False)
    filename = path + '\\stack_flat_end.tif'
    om.saveTiff(fn=filename)
    return True
    
def scan_zermode(om,dm,mode=10,amp=0.1,S=Sm.S,path=None,center=None):
    if path==None:
        path = om.path
    if center==None:
        center = om.zst.getPosition()
#    cmd_flat = dm.cmd #dm.read_command_file(dm.flat_file)
    cmd = [0.]*69
    # flat
    dm.Send(cmd)
    time.sleep(0.1)
    om.StackExt(center-3.0,center+3.0,0.2,verbose=False)
    filename = path + '\\stack_flat_start.tif'
    om.saveTiff(fn=filename)
    
    #zernike mode
    phiin = amp*zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
    dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
    for i in range(69):
        cmd[i] = dmarr[i]
    if (all(i <= 1.0 for i in dmarr)):
        dm.Send(cmd)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    time.sleep(0.1)
    om.StackExt(center-3.0,center+3.0,0.2,verbose=False)
    filename = path + '\\stack_act_zernikemode_%02d.tif' % mode
    om.saveTiff(fn=filename)
        
    # flat
    cmd = [0.]*69
    dm.Send(cmd)
    time.sleep(0.1)
    om.StackExt(center-3.0,center+3.0,0.2,verbose=False)
    filename = path + '\\stack_flat_end.tif'
    om.saveTiff(fn=filename)
    return True

def testSm(mode):
    S = Sm.S
    phiin = zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
    # set mirror with new dm shape
    dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
    phiout = N.dot(Sm.R,dmarr).reshape(33*33)
    return phiout

def setZernike(om, dm, mode, amp, S):
    phiin = N.zeros((33,33))
    cmd = [0.] * 69
    for j,m in enumerate(mode):
            phiin += amp[j]*zernike.Zm(m,rad=17,orig=None,Nx=33)
#    phiin = amp*zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
    # set mirror with new dm shape
    dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
    for i in range(69):
        cmd[i] = dmarr[i] + cmd_best[i]
    if (all(i <= 1.0 for i in dmarr)):
        dm.Send(cmd)
        print('PushDM_DONE')
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    return cmd_best 
    
def ao_zernike_scan(om, dm, mode, amp, S):
    phiin = amp*zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
    # set mirror with new dm shape
    dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
    cmd = [0.] * 69
    for i in range(69):
        cmd[i] = dmarr[i] + cmd_best[i]
    if (all(i <= 1.0 for i in dmarr)):
        dm.Send(cmd)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    return True
    
def setZernikeMode(om, dm, mode, amp, S):
    phiin = amp*zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
    # set mirror with new dm shape
    dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
    for i in range(69):
        cmd_best[i] = dmarr[i] + cmd_best[i]
    if (all(i <= 1.0 for i in dmarr)):
        dm.Send(cmd_best)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    return cmd_best
    
def ao_optimize(om,dm,mf,modes=N.arange(4,11,1),amprange=N.arange(-0.12,0.12,0.02),S=None,gain=50,blaser=False,rlaser=True,LED=1):
    center = om.zst.getPosition()
    results = []
    results.append(('Mode','Amp','Metric'))
    dm.Send(cmd_best)
    dt = N.zeros(amprange.shape)
    om.open_Acq(exposure=0.1,emgain=gain,laser=blaser,Rlaser=rlaser,
                         cntLaser=False, LED12=LED,FTM=False,conv=0,ccd=True,
                         trig=1,regen=False)
    om.get_img()
    for mode in modes:
        for k, amp in enumerate(amprange):
            print(k, amp)
            phiin = amp*zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
            # set mirror with new dm shape
            dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
            cmd = [0.] * 69
            for i in range(69):
                cmd[i] = dmarr[i] + cmd_best[i]
            if (all(i <= 1.0 for i in cmd)):
                dm.Send(cmd)
            else:
                raise Exception(' Error: push value greater than 1.0 ')
            time.sleep(0.02)
            om.zst.setPositionf(center)
            om.get_img()
            dt[k]= mf(om.data) # metric is peak intensity
            results.append((mode,amp,dt[k]))
        pmax = peak(amprange, dt)
        zmv[mode] += pmax
        print('setting mode %d at value of %f' % (mode, pmax))
        phiin = pmax*zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
        # set mirror with new dm shape
        dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
        for i in range(69):
            cmd_best[i] = cmd_best[i] + dmarr[i]
    if (all(i <= 1.0 for i in cmd_best)):
        dm.Send(cmd_best)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    om.ccd.AbortAcquisition()
    om.ccd.SetShutterMode(2)
    return results, cmd_best

def ao_optimize_si(om,dm,mf,modes=[4],amprange=N.arange(-1.0,1.25,0.25),S=None,gain=50,blaser=False,rlaser=True):
    center = om.zst.getPosition()
    dm.Send(cmd_best)
    dt = N.zeros(amprange.shape)
    results = []
    om.open_Acq(exposure=0.1,emgain=gain,laser=blaser,Rlaser=rlaser,
                         cntLaser=False, LED12=1,FTM=False,conv=0,ccd=True,
                         trig=1,regen=False)
    om.get_img()
    slm.selecteorder(16)
    slm.activate()
    for mode in modes:
        for k, amp in enumerate(amprange):
            print(k, amp)
            phiin = amp*zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
            # set mirror with new dm shape
            dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
            cmd = [0.] * 69
            for i in range(69):
                cmd[i] = dmarr[i] + cmd_best[i]
            if (all(i <= 1.0 for i in cmd)):
                dm.Send(cmd)
            else:
                raise Exception(' Error: push value greater than 1.0 ')
            time.sleep(0.02)
            om.zst.setPositionf(center)
            om.get_img()
            dt[k]= mf(om.data) # metric is peak intensity
            results.append((mode,amp,dt[k]))
        pmax = peak(amprange, dt)
        zmv[mode] += pmax
        print('setting mode %d at value of %f' % (mode, pmax))
        phiin = pmax*zernike.Zm(mode,rad=17, orig=None,Nx=33) # was 28
        # set mirror with new dm shape
        dmarr = 0.5*N.dot(S,phiin.reshape(33*33))
        for i in range(69):
            cmd_best[i] = cmd_best[i] + dmarr[i]
    if (all(i <= 1.0 for i in cmd_best)):
        dm.Send(cmd_best)
    else:
        raise Exception(' Error: push value greater than 1.0 ')

    slm.deactivate()
    om.ccd.AbortAcquisition()
    om.ccd.SetShutterMode(2)
    return results, cmd_best

def peak(x,y):
    a,b,c = N.polyfit(x, y, 2)
    zmax = -1*b/a/2.0
    if (a>0):
        print('no maximum')
        return 0.
    elif (zmax>=x.max()) or (zmax<=x.min()):
        print('maximum exceeding range')
        return 0.
    else:
        return zmax
    
def writedmfile(cmd,r):
    t = time.strftime("%Y%m%d%H%M%S")
    fns = t+'_flatfile.txt'
    with open(fns, 'w') as file:
        file.write(str(cmd))
    fn = t+'_modesvalue.txt'
    N.savetxt(fn, (mod,zmv))
    fns1 = t+'_metric.csv'
    with open(fns1, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(r)
    return True
    
def ao_optimize_pr(om,dm,S,center=None,path=None):
    ''' use phase retrieval to correct wavefront '''
    rad = 32 # radius for S matrix
    offset = 110#125
    if path == None:
        path = om.path
    if center == None:
        center = om.zst.getPosition()
#    if S == None:
#        raise Exception('No S matrix')
    # start
#    cmd_flat = dm.cmd
    # take stack
    om.setCCD_EM_ext(exposure=0.1, emgain=100, laser=False, Rlaser=False, LED12=2)    
    om.StackExt(center-1.6, center+1.6, 0.2, verbose=False)
    # do phase retrieval
    p = phaseret.mre(stack_cut(om.data))
    p.zcenter = 2.0
    p.dz = 0.2
    p.offset = offset
    p.run(64)
    if p.geterr()>300:
        print('convergence of phase retrieval is not great')
    # determine new dm settings
    amp,phase = p.get_amp_phase()
    pylab.figure()
    pylab.imshow(phase, interpolation='nearest')
    ph0 = phase[32:96,32:96]
    if False:
        amp0 = amp[32:96,32:96]
        msk = (amp0>0)
        target = 1.0*zernike.Zm(4,rad=30, orig=None,Nx=64)
        ph0 = (ph0 - target)*msk
        pylab.imshow(ph0,interpolation='nearest')
    dmarr = 0.5*N.dot(S,ph0.reshape(4*rad**2))
#    cmd = cmd_flat - dmarr
    cmd = [0.] * 69
    for i in range(69):
        cmd[i] = dmarr[i] + cmd_best[i]
    if (all(i <= 1.0 for i in cmd)):
        dm.Send(cmd)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
#    dm.set_mirror(cmd)
    # take stack
    #om.setCCD_Conv_ext(exposure=0.100,laser=False,Rlaser=False,LED12=1)
    om.StackExt(center-2.0, center+2.0, 0.2, verbose=False)
    # do phase retrieval
    p = phaseret.mre(stack_cut(om.data))
    p.zcenter = 2.0
    p.dz = 0.2
    p.offset = offset
    p.run(32)
    if p.geterr()>300:
        print('convergence of phase retrieval is not great')
    amp,phase = p.get_amp_phase()
    pylab.figure()
    pylab.imshow(phase, interpolation='nearest')
    return (cmd, cmd_best, p)
    
def stack_cut(stack,hsz=64,bgr = 900.):
    ''' cut out stack around maximum '''
    cz,cx,cy = N.unravel_index(stack.argmax(), stack.shape)
    stackout = N.fft.fftshift(stack[:,(cx-hsz):(cx+hsz),(cy-hsz):(cy+hsz)]-bgr, (1,2))
    return stackout

class ga_optimize(object):
    
    def __init__(self,om,dm,S):
        self.popN = 10
        Nz = 15
        self.Nz = Nz
        self.range = 2.0
        #self.aberr = 10.0*np.random.randn(N)
        self.population = N.zeros((self.popN,Nz))
        self.score = N.zeros(self.popN)
        self.om = om
        self.dm = dm
        self.dmflat = dm.cmd
        self.S = S
    
    def __del__(self):
        pass
    
    def create_pop(self):
        for m in range(self.popN):
            self.population[m,:] = self.range*N.random.randn(self.Nz)
    
    def score_gen(self):
        for m in range(self.popN):
            phiin = self.get_phi(self.population[m])
            # set mirror with new dm shape
            dmarr = 0.5*N.dot(self.S,phiin.reshape(64*64))
            self.dm.set_mirror(self.dmflat+dmarr)
            self.om.get_img()
            self.score[m] = self.om.data.max()
            
    def next_gen(self):
        ''' higher propability of mating with better score
            replace bottom with new members '''
        next_gen = N.zeros((self.popN, self.Nz))
        i = N.argsort(self.score)[::-1] # largest first
        top = int(self.popN/2)
        select = int(self.popN/4)
        for m in range(top):
            p1 = N.random.randint(select)
            p2 = N.random.randint(select)
            next_gen[m,:] = self.spawn(self.population[i[p1]], self.population[i[p2]])
        for m in range(top,self.popN):
            next_gen[m,:] = self.range*N.random.randn(self.Nz)
        self.population[:,:] = next_gen
        self.best = self.score[i[0]]
    
    def spawn(self,mom,dad):
        mask = N.random.randint(2, size=self.Nz)
        kid = mask*mom + (1-mask)*dad
        kid = kid + 0.1*N.random.randn(self.Nz)
        return kid
        
    def run(self,gens=100, exposure=0.05):
        self.create_pop()
        self.om.open_Acq(exposure,emgain=0,laser=False,Rlaser=False,
                         cntLaser=False, LED12=1,FTM=False,conv=True,ccd=True,
                         trig=1,regen=False)
        self.score_gen()
        for m in range(gens):            
            self.next_gen()
            self.score_gen()
            #if (m%100==0):
            print (m, self.best)
        self.om.ccd.AbortAcquisition()
        self.om.ccd.SetShutterMode(2)
            
    def set_best(self):
        i = N.argsort(self.score)[::-1]
        phi = self.get_phi(self.population[i[0]])
        dmarr = 0.5*N.dot(self.S,phi.reshape(64*64))
        self.dm.set_mirror(self.dmflat+dmarr)
                    
    def get_phi(self,gene):
        phi = N.zeros((64,64))
        for m in range(1,self.Nz):
            phi = phi + gene[m]*zernike.Zm(m, rad=32, orig=None, Nx=64)
        return phi