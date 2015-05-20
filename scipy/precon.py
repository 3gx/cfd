#!/usr/bin/python
import numpy as np
from numpy import linalg as LA

import sys
def printf(format, *args):
      sys.stdout.write(format % args)

def genRandomMatrix(size):
  return 2*np.random.random((size,size))-1
def genRandomEV(size,eps=0.1):
  return 1+ (2*np.random.random([size])-1)*eps

def genRandomPDMatrix(size,eps=0.1):
  P = genRandomMatrix(size)
  P = np.dot(P,P.transpose())
  ev = 1 + (2*np.random.random([size])-1)*eps
  return np.dot(P,np.dot(np.diag(ev),LA.inv(P)));

def isPositive(ev):
  return all(np.isreal(ev)) and all(ev > 0)

def transposeVec(A):
  return np.transpose(A,(0,2,1))

def matmulVec(A,B):
  msize = A[0].size
  niter = A.size/msize
  return np.sum(transposeVec(A).reshape(niter,size,size,1)*B.reshape(niter,size,1,size),-3)

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
if False:
  mm = np.array([
    [1.172029553823032,  0.4483838992489365,  -0.1216722269123576, 0.02804689929765592],
    [ -1.736281890383897,  0.4113037795103013,  0.6330178882990265,  -0.1216722269123576],
    [ 0.8262368847139853,  -1.285163043161573,  0.4113037795103013,  0.4483838992489365],
    [ -0.3759017444351098, 0.8262368847139853,  -1.736281890383897,  1.172029553823032]
    ])
if True:
  mm = np.array([
    [1.206348790012935 , 0.4391463979268352 , -0.1306066001318439, 0.0487844200342935 , -0.01226495874790658],
    [-1.845259182395017 , 0.4347623210981763 , 0.6736742077868462,  -0.2051201585243705, 0.0487844200342935],
    [0.9987488459399110,  -1.292927564706024 , 0.2844444444444444 , 0.6736742077868462,  -0.1306066001318439],
    [-0.5326716355661115, 0.6837488290237370,  -1.292927564706024 , 0.4347623210981763 , 0.4391463979268352],
    [0.2491918438040957,  -0.5326716355661115, 0.9987488459399110 , -1.845259182395017,  1.206348790012935]
    ])

size = mm[0].size;
niter = 1000000
pmat = []
evR  = []
eps = 0.01

#Q = np.empty([niter,size,size])
#E = np.empty([niter,size,size])
#for i in range(0,niter):
#  Q[i] = genRandomMatrix(size)
#  E[i] = np.diag(genRandomEV(size,eps))

#Q = [np.mat(2*np.random.random((size,size))-1) for i in range(niter)]
#niter=100
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

pmat = []
evR  = []
for i in range(REV.size):
  if REV[i] == True:
    ev = np.absolute(EV[i])
    evR  += [[np.max(ev)/np.min(ev),len(pmat)]]
    pmat += [P[i]]

evR.sort()
print " -------------- "
print evR[0:3]
print " -------------- "
print "pmat.size= ", len(pmat)
if len(pmat) > 0:
  k = evR[0][1]
  print "k=", k
  p = pmat[k]
  mmp = np.dot(p,mm)
#  mmp = np.dot(p,mm+np.identity(size))
  print "p=",   np.sort(LA.eigvals(p))
  print "mmp=", np.sort(LA.eigvals(mmp))
  print " -------------- "
  printf("{\n");
  for i in range(0,size):
    printf("{");
    for j in range (0,size):
      printf("%.15f ", p[i][j])
      printf(",") if (j < size-1) else printf(" ");
    printf("},\n") if (i!=j) else printf("}\n");
  printf("};\n");
