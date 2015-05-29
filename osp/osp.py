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
  zero_arr = np.zeros(1)
  for k in range(0,s-1):
    b[k+2][:] = 2*(m1*np.concatenate((zero_arr,b[k+1][0:-1])) + m0*b[k+1][:]) - b[k][:]

  c= np.zeros((s+1,len(z)))
  c[0][:] = 1;
  c[1][:] = m1*z+m0;

  for k in range(0,s-1):
    c[k+2][:] = 2*(m1*z + m0)*c[k+1][:] - c[k][:]

  return [b,c]

def minimizePoly(s,p,h,ev_space,maxiter=128,verbose=False):
  if verbose:
    print "============================================="
  hval = h*ev_space;
  [b,c] = scaled_chebyshev_basis(s,p,min(np.real(hval)),0,hval)

  fixed_coefficients = np.ones(p+1)/sp.misc.factorial(np.linspace(0,p,p+1))

  def func(x,c):
    return LA.norm(np.dot(x,c),ord=np.inf)-1;

  cons = ({'type': 'eq',
          'fun' : lambda x: np.dot(x,b) - fixed_coefficients})

  x0 = np.zeros(s+1)
  x0 = np.ones(s+1)
  res=optimize.minimize(func, x0, args=(c,),constraints=cons,
#      method='SLSQP', options={'disp': True, 'maxiter': 1024, 'ftol': 1e-13, 'eps':1e-13},tol=1e-15)
      method='SLSQP', options={'disp': verbose, 'maxiter': maxiter}, tol=1e-13)

  if verbose:
    print "------------------------------------"
    print res
    print "------------------------------------"
    print 'Value= ', func(res.x,c)
    print "coeff= "
    for x in np.dot(res.x,b):
      sys.stdout.write("%.16g," % x)
    print ""
    for x in fixed_coefficients:
      sys.stdout.write("%.16g," % x)
    print "\n============================================="

  if res.success:
    return [res.x,func(res.x,c)]
  else:
    return [None, None]

def maximizeH(s,p,ev_space):
  h_min = 0;
  h_max = 2.01*s*s*max(abs(ev_space))

  max_iter = 128;
  max_steps = 1000
  tol_bisect = 1e-3
  tol_feasible = 1e-12

  poly = None
  v    = None
  converged = False;
  for step in range(max_steps):
    if ((h_max-h_min < tol_bisect*h_min) or (h_max < tol_bisect)) and converged:
      break;

    h = 0.5*(h_max + h_min)

    [poly, v] = minimizePoly(s,p,h,ev_space,max_iter,verbose=False)

    if v == None:
      print "%5d  h_min= %g   h_max= %g  -- h= %g   " % (step, h_min, h_max, h)
    else:
      print "%5d  h_min= %g   h_max= %g  -- h= %g   v= %g " % (step, h_min, h_max, h, v)

    if v == None:
      converged = False
      print " >>>> Failed to converge "
      h_max = h;
    else:
      converged = True
      if v <= tol_feasible:
        h_min = h;
      else:
        h_max = h;


  if converged:
    print " Converged with h= %g" % h
    [poly, v] = minimizePoly(s,p,h,ev_space,max_iter,verbose=True)
    return [poly, h]
  else:
    return [None, None]


if True:
  npts = 10000;
  ev_space = -np.linspace(0,1,npts);
  s = 20;
  p = 8;

  [poly, h] = maximizeH(s,p,ev_space);
  if h != None:
    print "------ Polynomial coefficients -------- "
    for x in poly:
      sys.stdout.write("%.16g," % x)
    print ""
  else:
    print " ------- Not converged ------ "




if False:
  npts = 1000;
  ev_space = -np.linspace(0,1,npts);
  h = 3.165161132812499e+01;
  s = 10;
  h= 1.613250732421875e+03*0.99;
  s = 100;
  p = 8;
  [poly, v] = minimizePoly(s,p,h,ev_space,maxiter=128,verbose=True)

  if v != None:
    print "------ Polynomial coefficients -------- "
    for x in poly:
      sys.stdout.write("%.16g," % x)
    print ""
  else:
    print " ------- Not converged ------ "
