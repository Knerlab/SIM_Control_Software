# -*- coding: utf-8 -*-
"""
routines for ao optimization, etc.
stuff that was formerly in dm_setup

updated functions on 03/04/2020 by ruizhe

Created on Thu Feb  9 16:46:24 2017
@author: kner
"""
import time, os
import tifffile as tf
import numpy as N
# import Zernike36 as zernike
import pr_batch_P35 as prb
import circle as c
import csv

FLAT_FILE = r'C:/Users/Public/Documents/python_code/SIM_AO_p36/20210203150151_flatfile.txt'

r = 8
d = 14
a = 0.1

infl_fcn_file = r'C:/Users/Public/Documents/python_code/SIM_AO_p36/influ_func_02032021.tif'
infl_fcn = tf.imread(infl_fcn_file)
Sm = prb.Smat(infl_fcn, radius=r)

nz = 22
zern = N.zeros((nz,d,d))
zarr = N.zeros((nz,d,d))
for ii in range(nz):
    zarr[ii] = c.zernike_noll(ii + 1, d)
# orthogonalize
for ii in range(nz):
    zern[ii] = zarr[ii]
    for jj in range(ii):
        h = (zarr[ii]*zern[jj]).sum()
        zern[ii] = zern[ii] - h*zern[jj]
    zern[ii] = zern[ii]/N.sqrt((zern[ii]**2).sum())

with open(FLAT_FILE, "r") as file:
    cmd_best = eval(file.readline())
#cmd_best = [0.] * 69
cmd_zero = [0.] * 69
mod = N.arange(40)
zmv = N.zeros((40))

results = []

def dmnull(dm):
    cmd_best = cmd_zero
    if (all(i <= 1.0 for i in cmd_best)):
        dm.Send(cmd_best)
        print(' All actuators to 0 ')
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    return True

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
    # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d) # was 28
    # phiin = amp*c.zernike_noll(mode,d)
    phiin = amp*zern[mode]
    dmarr = a*N.dot(S,phiin.reshape((d*d)))
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

def testSm(mode, amp):
    S = Sm.S
    # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d)
    # phiin = amp*c.zernike_noll(mode,d)
    phiin = amp*zern[mode]
    # set mirror with new dm shape
    dmarr = a*N.dot(S,phiin.reshape(d*d))
    phiout = N.dot(Sm.R,dmarr).reshape(d,d)
    return phiout
    
def ao_zernike_set(dm, mode, amp, S):
    # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d)
    # phiin = amp*c.zernike_noll(mode,d)
    phiin = amp*zern[mode]
    # set mirror with new dm shape
    dmarr = a*N.dot(S,phiin.reshape((d*d)))
#    cmd = [0.] * 69
    for i in range(69):
        cmd_best[i] = dmarr[i] + cmd_best[i]
    if (all(i <= 1.0 for i in dmarr)):
        dm.Send(cmd_best)
        print('Zernike mode set on DM successfully')
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    return True

def set_corr_wf(om, dm, f, wf, zn, S, exclude=True):
    amp = N.zeros(zn)
    amp = f.decomwf(wf, zn)
    if (exclude):
        amp[:4] = 0
        # amp[10] = 0.
    phiin = N.zeros((d,d))
    phiin = f.recomwf(amp,d)
    results = []
    results.append(('Mode','Amp'))
    for i in range(zn):
        # phiin += amp[i] * zernike.Zm(z,rad=r,orig=None,Nx=d)
        results.append((i,amp[i]))
    phiin = -phiin
    dmarr = a*N.dot(S,phiin.reshape((d*d)))
    for i in range(69):
        cmd_best[i] = cmd_best[i] + dmarr[i]
    if (all(i <= 1.0 for i in dmarr)):
        dm.Send(cmd_best)
        print('Set new DM successfully')
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    return results, cmd_best

def set_new_wf(om, dm, f, fns, zn, S, exclude=True):
    amp = N.zeros(zn)
    wf = tf.imread(fns)
    amp = f.decomwf(wf, zn)
    if (exclude):
        amp[:4] = 0
    phiin = N.zeros((d,d))
    phiin = f.recomwf(amp,d)
    results = []
    results.append(('Mode','Amp'))
    for i in range(zn):
        # phiin += amp[i] * zernike.Zm(z,rad=r,orig=None,Nx=d)
        results.append((i,amp[i]))
    phiin = -phiin
    dmarr = a*N.dot(S,phiin.reshape((d*d)))
    for i in range(69):
        cmd_best[i] = cmd_best[i] + dmarr[i]
    if (all(i <= 1.0 for i in dmarr)):
        dm.Send(cmd_best)
        print('Set new DM successfully')
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    return results, cmd_best
   
def ao_optimize_snr(om,dm,path,center,mf,lpr,hpr,modes=N.arange(4,11,1),amprange=N.arange(-0.12,0.12,0.02),S=None,gain=50,Blaser=False,Rlaser=False,Ylaser=False,LED=1):
    t=time.localtime()    
    x = N.array([1e4,1e2,1])    
    t1 = int((t[0:3]*x).sum())
    t2 = int((t[3:6]*x).sum())
    fnt = "%s%s" %(t1,t2)
    newfold = path + '/' + fnt + '_ao_iteration' + '_snr' + '/'
    try:
        os.mkdir(newfold)
    except:
        print('Directory already exists')
    # center = om.zst.getPosition()
    results = []
    results.append(('Mode','Amp','Metric'))
    t=time.localtime()
    
    dm.Send(cmd_best)
    dt = N.zeros(amprange.shape)
    om.open_Acq(exposure=0.1,emgain=gain,Blaser=Blaser,Rlaser=Rlaser,Ylaser=Ylaser,
                         LED12=LED,FTM=False,conv=0,ccd=True,trig=1)
    om.get_img()
    for mode in modes:
        for k, amp in enumerate(amprange):
            # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = amp*c.zernike_noll(mode,d)
            phiin = amp*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
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
            fn = "zm%0.2d_amp%.4f" %(mode,amp)
            fn1 = os.path.join(newfold,fn+'.tif')
            tf.imsave(fn1,om.data)
            if (Rlaser==True):
                wl = 0.670
            if (Blaser==True):
                wl = 0.515
            if (Ylaser==True):
                wl = 0.570
            dt[k]= mf(wl,om.data,lpr,hpr) # metric is snr
            results.append((mode,amp,dt[k]))
            print(k, amp, dt[k])
        pmax = peak(amprange, dt)
        if (pmax!=0.0):
            zmv[mode] += pmax
            print('---------------------------------------------setting mode %d at value of %f' % (mode, pmax))
            # phiin = pmax*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = pmax*c.zernike_noll(mode,d)
            phiin = pmax*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
            for i in range(69):
                cmd_best[i] = cmd_best[i] + dmarr[i]
        else:
            print('---------------------------------------------mode %d value equals %f' % (mode, pmax))
    if (all(i <= 1.0 for i in cmd_best)):
        dm.Send(cmd_best)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    om.ccd.AbortAcquisition()
    om.ccd.SetShutterMode(2)
    return results, cmd_best
    
def ao_optimize_hf(om,dm,path,center,mf,modes=N.arange(4,11,1),amprange=N.arange(-0.12,0.12,0.02),S=None,gain=50,Blaser=False,Rlaser=False,Ylaser=False,LED=1):
    t=time.localtime()    
    x = N.array([1e4,1e2,1])    
    t1 = int((t[0:3]*x).sum())
    t2 = int((t[3:6]*x).sum())
    fnt = "%s%s" %(t1,t2)
    newfold = path + '/' + fnt + '_ao_iteration' + '_hpf' + '/'
    try:
        os.mkdir(newfold)
    except:
        print('Directory already exists')
    # center = om.zst.getPosition()
    results = []
    results.append(('Mode','Amp','Metric'))
    t=time.localtime()
    
    dm.Send(cmd_best)
    dt = N.zeros(amprange.shape)
    om.open_Acq(exposure=0.1,emgain=gain,Blaser=Blaser,Rlaser=Rlaser,Ylaser=Ylaser,
                         LED12=LED,FTM=False,conv=0,ccd=True,trig=1)
    om.get_img()
    for mode in modes:
        for k, amp in enumerate(amprange):
            # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = amp*c.zernike_noll(mode,d)
            phiin = amp*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
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
            fn = "zm%0.2d_amp%.4f" %(mode,amp)
            fn1 = os.path.join(newfold,fn+'.tif')
            tf.imsave(fn1,om.data)
            if (Rlaser==True):
                wl = 0.647
            if (Blaser==True):
                wl = 0.515
            if (Ylaser==True):
                wl = 0.570
            dt[k]= mf(wl,om.data) # metric is hpf
            results.append((mode,amp,dt[k]))
            print(k, amp, dt[k])
        pmax = peak(amprange, dt)
        if (pmax!=0.0):
            zmv[mode] += pmax
            print('---------------------------------------------setting mode %d at value of %f' % (mode, pmax))
            # phiin = pmax*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = pmax*c.zernike_noll(mode,d)
            phiin = pmax*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
            for i in range(69):
                cmd_best[i] = cmd_best[i] + dmarr[i]
        else:
            print('---------------------------------------------mode %d value equals %f' % (mode, pmax))
    if (all(i <= 1.0 for i in cmd_best)):
        dm.Send(cmd_best)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    om.ccd.AbortAcquisition()
    om.ccd.SetShutterMode(2)
    return results, cmd_best

def ao_optimize_sharp(om,dm,path,center,mf,modes=N.arange(4,11,1),amprange=N.arange(-0.12,0.12,0.02),S=None,gain=50,Blaser=False,Rlaser=False,Ylaser=False,LED=1):
    t=time.localtime()    
    x = N.array([1e4,1e2,1])    
    t1 = int((t[0:3]*x).sum())
    t2 = int((t[3:6]*x).sum())
    fnt = "%s%s" %(t1,t2)
    newfold = path + '/' + fnt + '_ao_iteration' + '_sharp' + '/'
    try:
        os.mkdir(newfold)
    except:
        print('Directory already exists')
    # center = om.zst.getPosition()
    results = []
    results.append(('Mode','Amp','Metric'))
    t=time.localtime()
    
    dm.Send(cmd_best)
    dt = N.zeros(amprange.shape)
    om.open_Acq(exposure=0.1,emgain=gain,Blaser=Blaser,Rlaser=Rlaser,Ylaser=Ylaser,
                         LED12=LED,FTM=False,conv=0,ccd=True,trig=1)
    om.get_img()
    for mode in modes:
        for k, amp in enumerate(amprange):
            # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = amp*c.zernike_noll(mode,d)
            phiin = amp*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
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
            fn = "zm%0.2d_amp%.4f" %(mode,amp)
            fn1 = os.path.join(newfold,fn+'.tif')
            tf.imsave(fn1,om.data)
            dt[k]= mf(om.data) # metric is sharpness
            results.append((mode,amp,dt[k]))
            print(k, amp, dt[k])
        pmax = peak(amprange, dt)
        if (pmax!=0.0):
            zmv[mode] += pmax
            print('---------------------------------------------setting mode %d at value of %f' % (mode, pmax))
            # phiin = pmax*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = pmax*c.zernike_noll(mode,d)
            phiin = pmax*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
            for i in range(69):
                cmd_best[i] = cmd_best[i] + dmarr[i]        
        else:
            print('---------------------------------------------mode %d value equals %f' % (mode, pmax))
    if (all(i <= 1.0 for i in cmd_best)):
        dm.Send(cmd_best)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    om.ccd.AbortAcquisition()
    om.ccd.SetShutterMode(2)
    return results, cmd_best
    
def ao_optimize_max(om,dm,path,center,mf,modes=N.arange(4,11,1),amprange=N.arange(-0.12,0.12,0.02),S=None,gain=50,Blaser=False,Rlaser=False,Ylaser=False,LED=1):
    t=time.localtime()    
    x = N.array([1e4,1e2,1])    
    t1 = int((t[0:3]*x).sum())
    t2 = int((t[3:6]*x).sum())
    fnt = "%s%s" %(t1,t2)
    newfold = path + '/' + fnt + '_ao_iteration' + '_max' + '/'
    try:
        os.mkdir(newfold)
    except:
        print('Directory already exists')
    # center = om.zst.getPosition()
    results = []
    results.append(('Mode','Amp','Metric'))
    t=time.localtime()
    
    dm.Send(cmd_best)
    dt = N.zeros(amprange.shape)
    om.open_Acq(exposure=0.1,emgain=gain,Blaser=Blaser,Rlaser=Rlaser,Ylaser=Ylaser,
                         LED12=LED,FTM=False,conv=0,ccd=True,trig=1)
    om.get_img()
    for mode in modes:
        for k, amp in enumerate(amprange):
            # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = amp*c.zernike_noll(mode,d)
            phiin = amp*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
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
            fn = "zm%0.2d_amp%.4f" %(mode,amp)
            fn1 = os.path.join(newfold,fn+'.tif')
            tf.imsave(fn1,om.data)
            dt[k]= mf(om.data) # metric is peak intensity
            results.append((mode,amp,dt[k]))
            print(k, amp, dt[k])
        pmax = peak(amprange, dt)
        if (pmax!=0.0):
            zmv[mode] += pmax
            print('---------------------------------------------setting mode %d at value of %f' % (mode, pmax))
            # phiin = pmax*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = pmax*c.zernike_noll(mode,d)
            phiin = pmax*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
            for i in range(69):
                cmd_best[i] = cmd_best[i] + dmarr[i]
        else:
            print('---------------------------------------------mode %d value equals %f' % (mode, pmax))
    if (all(i <= 1.0 for i in cmd_best)):
        dm.Send(cmd_best)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    om.ccd.AbortAcquisition()
    om.ccd.SetShutterMode(2)
    return results, cmd_best
    
def ao_optimize_sc(om,dm,path,center,modes=N.arange(4,11,1),amprange=N.arange(-0.12,0.12,0.02),S=None,gain=50,Blaser=False,Rlaser=False,Ylaser=False,LED=1):
    t=time.localtime()    
    x = N.array([1e4,1e2,1])    
    t1 = int((t[0:3]*x).sum())
    t2 = int((t[3:6]*x).sum())
    fnt = "%s%s" %(t1,t2)
    newfold = path + '/' + fnt + '_ao_iteration' + '/'
    try:
        os.mkdir(newfold)
    except:
        print('Directory already exists')    
    # center = om.zst.getPosition()
    dm.Send(cmd_best)
    om.open_Acq(exposure=0.1,emgain=gain,Blaser=Blaser,Rlaser=Rlaser,Ylaser=Ylaser,
                         LED12=LED,FTM=False,conv=0,ccd=True,trig=1)
    om.get_img()
    for mode in modes:
        print('mode', mode)
        for k, amp in enumerate(amprange):
            print(k, amp)
            # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = amp*c.zernike_noll(mode,d)
            phiin = amp*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape((d*d)))
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
            fn = "zm%0.2d_amp%.4f" %(mode,amp)
            fn1 = os.path.join(newfold,fn+'.tif')
            tf.imsave(fn1,om.data)
            dm.Send(cmd_best)
            time.sleep(0.02)
    om.ccd.AbortAcquisition()
    om.ccd.SetShutterMode(2)
    return True

def ao_optimize_si(om,dm,path,center,mf,slm,modes=N.arange(4,17,1),amprange=N.arange(-0.12,0.12,0.02),S=None,gain=50,Blaser=False,Rlaser=False,Ylaser=False,LED=1):
    t=time.localtime()    
    x = N.array([1e4,1e2,1])    
    t1 = int((t[0:3]*x).sum())
    t2 = int((t[3:6]*x).sum())
    fnt = "%s%s" %(t1,t2)
    newfold = path + '/' + fnt + '_ao_iteration' + '/'
    try:
        os.mkdir(newfold)
    except:
        print('Directory already exists')   
    # center = om.zst.getPosition()
    results = []
    results.append(('Mode','Amp','Metric'))
    dm.Send(cmd_best)
    dt = N.zeros(amprange.shape)
    om.open_Acq(exposure=0.1,emgain=gain,Blaser=Blaser,Rlaser=Rlaser,Ylaser=Ylaser,
                         LED12=LED,FTM=False,conv=0,ccd=True,trig=1)
    om.get_img()
    for m, mode in enumerate (modes):
        slm.selecteorder(30+m)
        slm.activate()   
        for k, amp in enumerate(amprange):
            print(k, amp)
            # phiin = amp*zernike.Zm(mode,rad=r,orig=None,Nx=d)
            # phiin = amp*c.zernike_noll(mode,d)
            phiin = amp*zern[mode]
            # set mirror with new dm shape
            dmarr = a*N.dot(S,phiin.reshape(d*d))
            cmd = [0.] * 69
            for i in range(69):
                cmd[i] = dmarr[i] + cmd_best[i]
            if (all(i <= 1.0 for i in cmd)):
                dm.Send(cmd)
            else:
                raise Exception(' Error: push value greater than 1.0 ')         
            om.zst.setPositionf(center)
            time.sleep(0.02)
            om.get_img()
            fn = "zm%0.2d_amp%.4f" %(mode,amp)
            fn1 = os.path.join(newfold,fn+'.tif')
            tf.imsave(fn1,om.data)
            if (Rlaser==True):
                wl = 0.647
            if (Blaser==True):
                wl = 0.515
            if (Ylaser==True):
                wl = 0.570
            dt[k]= mf(wl,om.data) # metric is peak intensity
            results.append((mode,amp,dt[k]))
        pmax = peak(amprange, dt)
        if (pmax!=0.0):
            zmv[mode] += pmax
            print('setting mode %d at value of %f' % (mode, pmax))
        else:
            print('mode %d value equals %f' % (mode, pmax))
        # phiin = pmax*zernike.Zm(mode,rad=r,orig=None,Nx=d)
        # phiin = pmax*c.zernike_noll(mode,d)
        phiin = pmax*zern[mode]
        # set mirror with new dm shape
        dmarr = a*N.dot(S,phiin.reshape(d*d))
        for i in range(69):
            cmd_best[i] = cmd_best[i] + dmarr[i]
        slm.deactivate() 
    if (all(i <= 1.0 for i in cmd_best)):
        dm.Send(cmd_best)
    else:
        raise Exception(' Error: push value greater than 1.0 ')
    om.ccd.AbortAcquisition()
    om.ccd.SetShutterMode(2)
    return results, cmd_best

def peak(x,y):
    a,b,c = N.polyfit(x, y, 2)
#    xp = N.linspace(x.min(), x.max(), 2000)
#    p = N.poly1d(z[0])
#    d = p(xp)
#    ind = N.unravel_index(N.argmax(d), d.shape)
#    zmax = xp[ind]
#    a = z[0]
    zmax = -1*b/a/2.0
    if (a>0):
        print('no maximum')
        return 0.
    elif (zmax>=x.max()) or (zmax<=x.min()):
        print('maximum exceeding range')
        return 0.
    else:
        return zmax
    
def writedmfile(pth,cmd,re):
    t = time.strftime("%Y%m%d%H%M%S")
    fns = t+'_flatfile.txt'
    fns = os.path.join(pth,fns)
    with open(fns, 'w') as file:
        file.write(str(cmd))
    fn = t+'_modesvalue.txt'
    fn = os.path.join(pth,fn)
    N.savetxt(fn, (mod,zmv))
    fns1 = t+'_metric.csv'
    fns1 = os.path.join(pth,fns1)
    with open(fns1, "w", encoding='utf-8',newline='') as f:
        writer = csv.writer(f)
        writer.writerows(re)
    return True
