#!/usr/bin/python
import numpy as np
from numpy import linalg as LA

import sys
def printf(format, *args):
      sys.stdout.write(format % args)

def genRandomMatrix(size):
  return 2*np.random.random((size,size))-1
def genRandomEV(size,eps):
  return 1+ (2*np.random.random([size])-1)*eps

def genRandomPDMatrix(size,eps):
  P = genRandomMatrix(size)
  P = np.dot(P,P.transpose())
  ev = 1 + (2*np.random.random([size])-1)*eps
  return np.dot(P,np.dot(np.diag(ev),LA.inv(P)));

def isPositive(ev):
  return all(np.isreal(ev)) and all(ev > 0)

def transposeVec(A):
 return np.einsum('...ji',A)
#  return np.transpose(A,(0,2,1))

def matmulVec(A,B):
 return np.einsum('...ij,...jk',A,B)
#  msize = A[0][0].size
#  niter = A.size/(msize*msize)
#  return np.sum(transposeVec(A).reshape(niter,msize,msize,1)*B.reshape(niter,msize,1,msize),-3)

def find_preconditioner(mm, niter, eps):
  # numpy acceleration ideas taken from
  # https://jameshensman.wordpress.com/2010/06/14/multiple-matrix-multiplication-in-numpy/
  size = mm[0].size;
  A = 2*np.random.random((niter,size,size))-1
  B = transposeVec(A)
  Q = matmulVec(A,B)
  Qinv = LA.inv(Q)
  l = np.diag(1 + (2*np.random.random([size])-1)*eps)
  P = matmulVec(np.dot(Q,l),Qinv); 
  Pm = np.dot(P,mm)
  #Pm = np.dot(P,mm+np.identity(size))
  EV = LA.eigvals(Pm)
  REV = np.logical_and(np.all(np.isreal(EV),1), np.all(EV > 0,1))

  pmat = None
  evR  = 1e10;
  for i in range(REV.size):
    if REV[i] == True:
      ev = np.absolute(EV[i])
      evT  = np.max(ev)/np.min(ev)
      if (evT < evR):
        pmat = P[i]
        evR  = evT

  return pmat;


mm = np.array([
  [1.,  0.36602540378443865],
  [-1.3660254037844386,  1.]
  ]);
if True:
  mm = np.array([
    [ 1.1111111111111111111,  0.4485512379056265999, -0.0808317913155015634],
    [-1.5596623490167377110,  0.4444444444444444444,  0.4485512379056265999],
    [ 0.6363873468710571190, -1.5596623490167377110,  1.1111111111111111111]
    ])
if True:
  mm = np.array([
    [1.172029553823032,  0.4483838992489365,  -0.1216722269123576, 0.02804689929765592],
    [ -1.736281890383897,  0.4113037795103013,  0.6330178882990265,  -0.1216722269123576],
    [ 0.8262368847139853,  -1.285163043161573,  0.4113037795103013,  0.4483838992489365],
    [ -0.3759017444351098, 0.8262368847139853,  -1.736281890383897,  1.172029553823032]
    ])
if False:
  mm = np.array([
    [1.206348790012935 , 0.4391463979268352 , -0.1306066001318439, 0.0487844200342935 , -0.01226495874790658],
    [-1.845259182395017 , 0.4347623210981763 , 0.6736742077868462,  -0.2051201585243705, 0.0487844200342935],
    [0.9987488459399110,  -1.292927564706024 , 0.2844444444444444 , 0.6736742077868462,  -0.1306066001318439],
    [-0.5326716355661115, 0.6837488290237370,  -1.292927564706024 , 0.4347623210981763 , 0.4391463979268352],
    [0.2491918438040957,  -0.5326716355661115, 0.9987488459399110 , -1.845259182395017,  1.206348790012935]
    ])

if True:  #  Chebyshev
  mm = np.array([
  [3.129901731909032, 0.6574252126145215, -0.1978177134420483,  0.14352519875080175,  -0.2903266038658256],
  [   -4.5314121905019755,  0.2925171707277256, 0.844425286270601,  -0.4944796688377792,  0.9402786498365101],
  [   2.1442539109393475, -1.341687376966689, 0.04000000000000042,  1.2428019900214273, -1.8853686041100468],
  [   -1.0813931430457429,  0.5744796630598727, -0.9433106754704905,  -0.12805900410862886, 4.032526994191361],
  [   0.37032657876534397,  -0.2846397412794731,  0.4567031233749075, -1.1563106002056767,  -1.5343601365022363]
  ]
  )


niter = 1000000
niter = 1000
eps = 0.01

P = None
evR = 1e10;
kk = 0;
while True:
#  if (kk == 1):
#    break
  kk += 1
  update = False;
  p = find_preconditioner(mm,niter,eps)
  if p != None:
    mmp = np.dot(p,mm)
    ev = LA.eigvals(mmp)
    evT = np.max(ev)/np.min(ev)
    if evT < evR:
      evR = evT
      P   = p
      update = True;

  if update == True:
    print " ########################### "
    print " ########################### "
    p = P
    size = p[0].size
    mmp = np.dot(p,mm)
  #  mmp = np.dot(p,mm+np.identity(size))
    print "kk =" , kk
    print "p=",   np.sort(LA.eigvals(p))
    ev = LA.eigvals(mmp)
    print "mmp=", np.sort(ev)
    print "range= ", np.max(ev)/np.min(ev)
    print " -------------- "
    printf("{\n");
    for i in range(0,size):
      printf("{");
      for j in range (0,size):
        printf("%.15f ", p[i][j])
        printf(",") if (j < size-1) else printf(" ");
      printf("},\n") if (i!=j) else printf("}\n");
    printf("};\n");
