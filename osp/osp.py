import numpy as np
from numpy import linalg as LA
from  scipy import optimize 
import scipy as sp
import sys



def scaled_chebyshev_basis(s,p,zmin,zmax,z):
  b = np.zeros((s+1,p+1))
  m1 = 2/(zmax-zmin)
  m0 = -(1+zmin*m1);

  # T_0 = 1
  b[0][0] = 1;

  # T_1 = m1*x + m0
  b[1][0] = m0;
  b[1][1] = m1;

  # T_{n+1} = 2*(m1*x + m0)*T_{n} - T_{n-1}
  # L_{n} = (2n-1)/n*(m1*x+m0)*L_{n} - (n-1)/n*L_{n-1}
  zero_arr = np.zeros(1)
  for k in range(0,s-1):
#    b[k+2][:] = (2.0*(k+2)-1)/(k+2)*(m1*np.concatenate((zero_arr,b[k+1][0:-1])) + m0*b[k+1][:]) - (k+2-1.0)/(k+2)*b[k][:]
    b[k+2][:] = 2*(m1*np.concatenate((zero_arr,b[k+1][0:-1])) + m0*b[k+1][:]) - b[k][:]

  c= np.zeros((s+1,len(z)),dtype=z.dtype)
  c[0][:] = 1;
  c[1][:] = m1*z+m0;

  for k in range(0,s-1):
#    c[k+2][:] = (2.0*(k+2)-1)/(k+2)*(m1*z + m0)*c[k+1][:] - (k+2.0-1)/(k+2)*c[k][:]
    c[k+2][:] = 2*(m1*z + m0)*c[k+1][:] - c[k][:]

  return [b,c]

def minimizePoly(s,p,h,ev_space,eta,tol,maxiter=128,verbose=False,poly_guess=None):
  if verbose:
    print "============================================="
  hval  = h*np.real(ev_space) + 1j*np.imag(ev_space)
  hval = h*ev_space;
  [b,c] = scaled_chebyshev_basis(s,p,min(np.real(hval)),0,hval)


  fixed_coefficients = np.ones(p+1)/sp.misc.factorial(np.linspace(0,p,p+1))

  def func(x):
    return x[-1];
  def func_jac(x):
    return np.concatenate((np.zeros(len(x)-1),[1]))

  def cfunc1(x,b,coeff):
    return np.dot(x[:-1],b) - coeff
  def cfunc1_jac(x,J):
    return J;

  if True:
    def cfunc2(x,c):
      g = np.dot(x[:-1],c);
      f = np.real(g*np.conj(g)) - 1
      imax = np.argmax(f)
      return np.real(x[-1] - f[imax])
    def cfunc2_jac(x,c):
      g  = np.dot(x[:-1],c[:-1]);
      f = np.real(g*np.conj(g)) - 1
      imax = np.argmax(f)
      ct = c.T;
      df = -ct[imax]*np.conj(g[imax]);
      df[-1] = 0.5;
      return np.real(df + np.conj(df)).T
  else:
    def cfunc2(x,c):
      g = np.dot(x[:-1],c);
      return np.real(x[-1]*np.ones(len(g)) - (g*np.conj(g) - 1))
    def cfunc2_jac(x,c):
      g  = np.dot(x[:-1],c[:-1]);
      df = -c*np.conj(g);
      df[-1] = 0.5;
      return np.real(df + np.conj(df)).T


  m,n = b.shape
  m,k = c.shape
  J1 = np.zeros((s+2,n))
  J2 = np.zeros((s+2,k),dtype=c.dtype)

  for i in range(s+1):
    J1[i] = b[i]
    J2[i] = c[i]
  J1 = J1.T

  x0 = np.zeros(s+2);
  cons = ({'type': 'eq',
    'fun' : lambda x: cfunc1(x,b,fixed_coefficients),
    'jac' : lambda x: cfunc1_jac(x,J1)},
    {'type': 'ineq',
    'fun' : lambda x: cfunc2(x,c)})
  #  'jac' : lambda x: cfunc2_jac(x,J2)})
    
  res=optimize.minimize(func, x0, constraints=cons, jac=func_jac,
      method='SLSQP', options={'disp': verbose, 'maxiter': maxiter}, tol=1e-13)

  if verbose:
    print "------------------------------------"
    print res
    print "------------------------------------"
    print 'Value= ', res.fun-eta
    print "coeff= "
    for x in np.dot(res.x[:-1],b):
      sys.stdout.write("%.16g," % x)
    print ""
    for x in fixed_coefficients:
      sys.stdout.write("%.16g," % x)
    print "\n============================================="


  if res.success:
    return [True, res.x[:-1], res.fun, res.nit]
  elif abs(res.fun-eta) < tol:
    return [True, res.x[:-1], res.fun, res.nit]
  else:
    return [False, res.x, res.fun, res.nit]

def maximizeH(s,p,ev_space):
  h_min = 0; #60 #0.00*max(abs(ev_space))
  h_max = 2.01*s*s*max(abs(ev_space))

  max_iter = 1280;
  max_steps = 1000
  tol_bisect = 0.001
  tol_feasible = 1.0e-12
  eta = 0.0

  print "max_iter= %d " % max_iter
  print "max_steps= %d " % max_steps
  print "tol_bisect= %g " % tol_bisect
  print "tol_feasible= %g " % tol_feasible
  print "eta= %g" % eta

  poly = None
  v    = None
  converged = False;
  for step in range(max_steps):
    if ((h_max-h_min < tol_bisect*h_min) or (h_max < tol_bisect)):
      if converged:
        break;
      else:
        h = h_min
    else:
      h = 0.25*h_max + 0.75*h_min
      h = 0.5*h_max + 0.5*h_min


#    h = 0.5*h_max + 0.5*h_min

    [conv, poly, v, nit] = minimizePoly(s,p,h,ev_space,eta,tol_feasible,max_iter,verbose=False)
    print "%5d  h_min= %g   h_max= %g  -- h= %g nit= %d  v= %g " % (step, h_min, h_max, h, nit, v-eta)

    if False:
      if not conv:
        converged = False
        print " >>>> Failed to converge "
        h_max = h;
      else:
        converged = True
        if abs(v-eta) <= tol_feasible:
          h_min = h;
        else:
          h_max = h;
    else:
      if not conv or abs(v-eta) > tol_feasible:
        converged = False
        print " >>>> Failed to converge "
        h_max = h;
      else:
        converged = True
        h_min = h;


  if True or converged:
    print " Converged with h= %g  h/s^2= %g" % (h, h/s**2)
    [conv, poly, v, nit] = minimizePoly(s,p,h,ev_space,eta,tol_feasible,max_iter,verbose=True)
    return [poly, h]
  else:
    return [None, None]


if False:
  npts = 1000;

  ev_space = -np.linspace(0,1,npts);

#  ev_space = np.random.random(npts);
#  ev_space = np.sort(np.concatenate(([0],ev_space,[1])))
#  ev_space = -0.5*(1 + np.cos(ev_space*np.pi))

#  ev_space = np.linspace(0,np.pi,npts)
#  ev_space = -0.5*(1 + np.cos(ev_space))
  ev_space = -np.linspace(0,1,npts);

  if True:
    kappa=1;
#    beta =5.0;
    beta = 0.2;

    imag_lim = beta;
    l1 = 1j*np.linspace(0,imag_lim,npts/2);
    l2 = 1j*imag_lim + np.linspace(-kappa,0,npts);
    l3 = -kappa + 1j*np.linspace(0,imag_lim,npts/2);
    ev_space = np.concatenate((l1,l2,l3));


  s = 30;
  p = 8;

  print "npts= %d  s= %d  p= %d " % (npts, s, p)

  [poly, h] = maximizeH(s,p,ev_space);
  if h != None:
    print "------ Polynomial coefficients -------- "
    print "coeff= {"
    for x in poly[:-1]:
      sys.stdout.write("%.16g," % x)
    sys.stdout.write("%.16g};\n" % poly[-1])
    print "h= %.16g  h/s^2= %g " % (h, h/(s*s))
  else:
    print " ------- Not converged ------ "

  sys.exit(-1)




if True:
  npts = 1000;
  kappa=1;
  beta = 0.5;
  beta= 0.5;

  imag_lim = beta;
  l1 = 1j*np.linspace(0,imag_lim,50);
  l2 = 1j*imag_lim + np.linspace(-kappa,0,npts);
  l3 = -kappa + l1 #1j*np.linspace(0,imag_lim,50);
  ev_space = np.concatenate((l1,l2,l3));
#  ev_space =  np.linspace(-kappa,0,npts);
  eta=0
  tol=1.0e-12
  h=35.58

  p=4
  s=30

  maxiter=1024
  [conv, poly, v, nit] = minimizePoly(s,p,h,ev_space,eta,tol,maxiter,verbose=True)

  if not conv:
    print " ------- Not converged ------ "

  print "------ Polynomial coefficients -------- "
  print "coeff= {"
  for x in poly[:-1]:
    sys.stdout.write("%.16g," % x)
  sys.stdout.write("%.16g};\n" % poly[-1])
  print "h= %.16g  h/s^2= %g " % (h, h/(s*s))
