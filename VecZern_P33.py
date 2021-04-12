''' orthonormal vector polynomials on unit circle
    Zhao and Burge, Optics Express 15(26),p. 18014 (2007) '''
import os, sys
#path,fnp = os.path.split(__file__)
#Drive = path[0]
#sys.path.append(Drive+':\\python')
##from Priithon.all import N,F,Y
import numpy as N
#import fftfuncs as F
import tifffile as T
import Zernike36 as zernike
from scipy.special import factorial as fac
from operator import itemgetter

Nzern = 37 

def mult(x,y):
    if x == 0:
        return 0
    else:
        return x*y

#def mkl(mx):
#    a = []
#    for n in range(mx):
#        for m in range(-n,n+1,2):
#            a.append((n+N.abs(m),n,m))
#    a.sort(lambda a,b: int(a[0]-b[0] or a[1]-b[1] or b[2]-a[2]))
#    c = []
#    for j,t in enumerate(a):
#        c.append(t[1:])
#    b = N.array(c)
#    #return dict(zip(c,N.arange(mx)))
#    return (b,dict(zip(c,N.arange(len(a)))))

def jnm(mx):
    ''' Noll Zernike ordering !!! '''
    a = []
    for n in range(mx):
        for m in range(-n,n+1,2):
            a.append((n,m,abs(m)))
    #a.sort(lambda a,b: int(a[0]-b[0] or a[2]-b[2] or a[1]-b[1]))
    a.sort(key=itemgetter(0,2,1))
    c = []
    for j,t in enumerate(a):
        c.append(t[:2])
    b = N.array(c)
    #return dict(zip(c,N.arange(mx)))
    return (b,dict(zip(c,N.arange(len(a)))))
    
global nj, nb
nb,nj = jnm(15)

def getS(j, verbose=False):
    n,m = nb[j]
    if verbose: print(n,m)
    sq05 = N.sqrt(1./2.)
    s1 = [2*s-1 for s in range(1,21,2)]
    s2 = [2*s-1 for s in range(2,21,2)]
    am = abs(m)
    if (n==1) and (am==1):
        if m>0:
            if verbose: print("x (%d,%d)" % (n-1,-am+1))
            if verbose: print("y NA")
            ca1,ca2,cb1,cb2 = 1,0,0,0
        else:
            if verbose: print("x NA")
            if verbose: print("y (%d,%d)" % (n-1,am-1))
            ca1,ca2,cb1,cb2 = 0,0,1,0
    elif (m==0):
        if verbose: print("x 0.7(%d,%d)" % (n-1,am+1))
        if verbose: print("y 0.7(%d,%d)" % (n-1,-am-1))
        ca1,ca2,cb1,cb2 = 0,sq05,0,sq05
    elif (m==1):
        if verbose: print("x 0.7(%d,%d) 0.5(%d,%d)" % (n-1,am-1,n-1,am+1))
        if verbose: print("y 0.5(%d,%d)" % (n-1,-am-1))
        ca1,ca2,cb1,cb2 = sq05,0.5,0,0.5
    elif (m==-1):
        if verbose: print("x 0.5(%d,%d)" % (n-1,-am-1))
        if verbose: print("y 0.7(%d,%d) -0.5(%d,%d)" % (n-1,am-1,n-1,am+1))
        ca1,ca2,cb1,cb2 = 0,0.5,sq05,-0.5
    elif (m==n):
        if verbose: print("x 0.5(%d,%d)" % (n-1,am-1))
        if verbose: print("y -0.5(%d,%d)" % (n-1,-am+1))
        ca1,ca2,cb1,cb2 = sq05,0,-sq05,0
    elif (m==-n):
        if verbose: print("x 0.5(%d,%d)" % (n-1,-am+1))
        if verbose: print("y 0.5(%d,%d)" % (n-1,am-1))
        ca1,ca2,cb1,cb2 = sq05,0,sq05,0
    elif (m>0):
        if verbose: print("x 0.5(%d,%d) 0.5(%d,%d)" % (n-1,am-1,n-1,am+1))
        if verbose: print("y -0.5(%d,%d) 0.5(%d,%d)" % (n-1,-am+1,n-1,-am-1))
        ca1,ca2,cb1,cb2 = 0.5,0.5,-0.5,0.5
    elif (m<0):
        if verbose: print("x 0.5(%d,%d) 0.5(%d,%d)" % (n-1,-am+1,n-1,-am-1))
        if verbose: print("y 0.5(%d,%d) -0.5(%d,%d)" % (n-1,am-1,n-1,am+1))
        ca1,ca2,cb1,cb2 = 0.5,0.5,0.5,-0.5
    ######################
    if (m<0):
        sx = lambda nx,ny,radius: (mult(ca1,NormZern(nx,ny,radius,n-1,-am+1))
                + mult(ca2,NormZern(nx,ny,radius,n-1,-am-1)))
        sy = lambda nx,ny,radius: (mult(cb1,NormZern(nx,ny,radius,n-1,am-1))
                + mult(cb2,NormZern(nx,ny,radius,n-1,am+1)))
    else:
        sx = lambda nx,ny,radius: (mult(ca1,NormZern(nx,ny,radius,n-1,am-1))
                + mult(ca2,NormZern(nx,ny,radius,n-1,am+1)))
        sy = lambda nx,ny,radius: (mult(cb1,NormZern(nx,ny,radius,n-1,-am+1))
                + mult(cb2,NormZern(nx,ny,radius,n-1,-am-1)))
    return (sx,sy)

def getphi(j,verbose=False):
    n,m = nb[j]
    if verbose: print(n,m)
    jp = nj.get((n-2,m),0)
    if jp==0:
        phi = lambda nx,ny,radius: (1./N.sqrt(4*n*(n+1)))*NormZern(nx,ny,radius,n,m)
    else:
        phi = lambda nx,ny,radius: (1./N.sqrt(4*n*(n+1)))*(NormZern(nx,ny,radius,n,m)
                                - N.sqrt((n+1)/(n-1))*NormZern(nx,ny,radius,n-2,m))
    return phi

def testvec(j):
    nx = 256
    ny = 256
    rad = 100
    sx,sy = getS(j)
    T.imshow(F.zzernikeArr(shape=(nx,ny),no=j,radius=rad))
    T.imshow(sx(nx,ny,rad))
    T.imshow(sy(nx,ny,rad))

def testortho(j1,j2,verbose=False):
    nx = 256
    ny = 256
    rad = 100
    sx1,sy1 = getS(j1)
    sx2,sy2 = getS(j2)
    dp = (sx1(nx,ny,rad)*sx2(nx,ny,rad) + sy1(nx,ny,rad)*sy2(nx,ny,rad))
    if verbose: T.imshow(dp)
    return (dp.sum()/(N.pi*rad**2)) #(nx*ny))

def orthomat(jm):
    q = F.zeroArrF((jm,jm))
    for m in range(1,jm):
        for n in range(1,jm):
            q[m,n] = testortho(m,n)
    T.imshow(q)
    return q

def testZernOrtho(zn):
    nx = 256
    rad = 100
    qn = F.zeroArrF(zn,zn)
    msk = NormZern(nx,nx,rad,0,0)
    norm = msk.sum()
    for j1 in range(zn):
        for j2 in range(zn):
            n1,m1 = nb[j1]
            t1 = NormZern(nx,nx,rad,n1,m1)
            n2,m2 = nb[j2]
            t2 = NormZern(nx,nx,rad,n2,m2)
            qn[j1,j2] = (t1*t2).sum()/norm
    return qn

def NormZern(nx,ny,radius,n,m):
    global nj
    if abs(m)>n:
        return N.nan
    else:
        if m==0:
            fac = N.sqrt(n+1)
        else:
            fac = N.sqrt(2*(n+1))
        j = nj.get((n,m),0)
#    out = fac*F.zzernikeArr(shape=(nx,ny),no=j,radius=radius)
    #out = fac*F.ringArr((nx,ny),radius1=0,radius2=radius)
#        print n,m
        out = zernike.Z(m,n,radius,None,nx)
        return out

def convj(normzernarr):
    nl = len(normzernarr)
    zernarr = F.zeroArrF(nl)
    for j in range(nl):
        n,m = nb[j]
        if m==0:
            fac = N.sqrt(n+1)
        else:
            fac = N.sqrt(2*(n+1))
        zernarr[j] = (fac)*normzernarr[j]
    return zernarr

def diff(phi,rad):
    #dx = N.diff(phi,1,0)
    #dy = N.diff(phi,1,1)
    tpi = 2*N.pi
    nx,ny = phi.shape
    cntr=True
#    msk = F.zzernikeArr(shape=phi.shape,no=0,crop=1,radius=rad-1) #rad-1
    msk = discArray(phi.shape,radius=rad-1)#F.ringArr(phi.shape,radius1=0,radius2=rad-1)
    if cntr:
        ind1 = N.arange(-1,nx-1)%nx
        ind2 = N.arange(1,nx+1)%nx
        dx = (phi[ind2,:] - phi[ind1,:])
        dy = (phi[:,ind2] - phi[:,ind1])
        dy = msk*(N.mod(dy+N.pi,tpi)-N.pi)*rad/2.
        dx = msk*(N.mod(dx+N.pi,tpi)-N.pi)*rad/2.
    else:
        ind1 = N.arange(nx)
        ind2 = N.arange(1,nx+1)%nx
        dy = -1*(phi[ind1,:] - phi[ind2,:])
        dx = phi[:,ind2] - phi[:,ind1]
        dy = msk*(N.mod(dy+N.pi,tpi)-N.pi)*rad
        dx = msk*(N.mod(dx+N.pi,tpi)-N.pi)*rad
    return (dx,dy) # was (dx,dy) !

def getZc(t):
    nt = len(t)#/2 #! was 32
    g = N.zeros(nt)
    for j in range(1,nt):
        n,m = nb[j]
        jp = nj.get((n+2,m),0)
        if jp>(len(t)-1):
            jp = 0
        g[j] = t[j]*bb(j,j) + t[jp]*bb(jp,j)
    return g

def bb(j1,j2):
    n1,m1 = nb[j1]
    n2,m2 = nb[j2]
    if (n1==n2) and (m1==m2):
        if abs(m1)==n1:
            mel = 1./N.sqrt(2*n1*(n1+1))
        else:
            mel = 1./N.sqrt(4*n1*(n1+1))
    elif (n2==n1-2) and (m1==m2):
        mel = -1./N.sqrt(4*n1*(n1-1))
    else:
        mel = 0.0
    return mel
    
def discArray(shape=(128,128),radius=64,origin=None,dtype=N.float64):
    nx = shape[0]
    ny = shape[1]
    ox = nx/2
    oy = ny/2
    x = N.linspace(-ox,nx-ox,nx)
    y = N.linspace(-oy,ny-oy,ny)
    X,Y = N.meshgrid(x,y)
    rho = N.sqrt(X**2 + Y**2)
    disc = (rho<radius).astype(dtype)
    if not origin==None:
        s0 = origin[0]-int(nx/2)
        s1 = origin[1]-int(ny/2)
        disc = N.roll(N.roll(disc,s0,0),s1,1)
    return disc

def buildphiZ(carr,shape=(256,256),rad=128):
    ''' build phi from Zernike coeffs '''
    nx,ny = shape
    phi = N.zeros(shape)
    nz = len(carr)
    for j in range(nz):
        n,m = nb[j]
        #phi += carr[j]*NormZern(nx,ny,rad,n,m)
        phi += carr[j]*zernike.Z(m,n,rad,None,nx)
    return phi
    
def VecZernDecomp(bpp,nx,radius,verbose=False,phase=True):
    ''' the output is the amplitude of the different Zernike components in radians
        where the Zernikes are normalized to an RMS amplitude of 1 '''
#    if p==None:
#        nx = bpp.shape[0]
#        radius = nx/2
#    else:
#        nx = p.Nx#params['Nx']
#        dx = p.dx#params['dx']
#        wl = p.wl#params['wl']
#        nap = p.na #params['na']
#        n2 = p.n2# params['n2']
#        dp = 1/(nx*dx)
#        radius = (2*nap/wl)/2/dp
#        factor = nx/2./radius
    #########################
    if phase:
        phi = bpp
    else:
        phi = N.angle(bpp)
    dx,dy = diff(phi,radius)
    if verbose:
        T.imshow(dx, vmax = dx.max(), vmin=dx.min())
        T.imshow(dy, vmax = dx.max(), vmin=dx.min())
    scoeff = [0]
    for j in range(1,Nzern):
        sx,sy = getS(j)
        t = (sx(nx,nx,radius)*dx - sy(nx,nx,radius)*dy).sum()/(N.pi*(radius)**2) #changed to minus
        if verbose: print(j, t)
        scoeff.append(t)
    return N.array(getZc(scoeff))