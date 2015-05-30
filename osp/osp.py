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

  c= np.zeros((s+1,len(z)),dtype=z.dtype)
  c[0][:] = 1;
  c[1][:] = m1*z+m0;

  for k in range(0,s-1):
    c[k+2][:] = 2*(m1*z + m0)*c[k+1][:] - c[k][:]

  return [b,c]

def minimizePoly(s,p,h,ev_space,eta,tol,maxiter=128,verbose=False,poly_guess=None):
  if verbose:
    print "============================================="
  hval  = h*np.real(ev_space) + 1j*np.imag(ev_space)
#  hval = h*ev_space;
  [b,c] = scaled_chebyshev_basis(s,p,min(np.real(hval)),0,hval)



  fixed_coefficients = np.ones(p+1)/sp.misc.factorial(np.linspace(0,p,p+1))

  def func(x,c,eta):
    g = np.dot(x,c);
    f = g*np.conj(g)
    imax = np.argmax(f-(1-eta))
    return np.real(f[imax])-(1-eta);

  def func_jac(x,c,eta):
    g = np.dot(x,c)
    f = g*np.conj(g)
    imax = np.argmax(f-(1-eta))
    ct = c.T
    df = ct[imax]*np.conj(g[imax]) + np.conj(ct[imax])*g[imax];
    return np.real(df)

  def cfunc(x,b,coeff):
    return np.dot(x,b) - coeff

  def cfunc_jac(x,b):
    return b.T

  cons = ({'type': 'eq',
    'fun' : lambda x: cfunc(x,b,fixed_coefficients),
    'jac' : lambda x: cfunc_jac(x,b)})
  
  x0 = np.zeros(s+1)
  x0 = np.ones(s+1)
  res=optimize.minimize(func, x0, args=(c,eta),constraints=cons,jac=func_jac,
      method='SLSQP', options={'disp': verbose, 'maxiter': maxiter}, tol=1e-13)

  if verbose:
    print "------------------------------------"
    print res
    print "------------------------------------"
    print 'Value= ', res.fun-eta
    print "coeff= "
    for x in np.dot(res.x,b):
      sys.stdout.write("%.16g," % x)
    print ""
    for x in fixed_coefficients:
      sys.stdout.write("%.16g," % x)
    print "\n============================================="


  if res.success:
    return [True, res.x, res.fun, res.nit]
  elif abs(res.fun-eta) < tol:
    return [True, res.x, res.fun, res.nit]
  else:
    return [False, res.x, res.fun, res.nit]

def maximizeH(s,p,ev_space):
  h_min = 0; #60 #0.00*max(abs(ev_space))
  h_max = 2.01*s*s*max(abs(ev_space))

  max_iter = 128;
  max_steps = 1000
  tol_bisect = 1e-3
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


#    h = 0.5*h_max + 0.5*h_min

    [conv, poly, v, nit] = minimizePoly(s,p,h,ev_space,eta,tol_feasible,max_iter,verbose=False)
    print "%5d  h_min= %g   h_max= %g  -- h= %g nit= %d  v= %g " % (step, h_min, h_max, h, nit, v-eta)

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


  if True or converged:
    print " Converged with h= %g  h/s^2= %g" % (h, h/s**2)
    [conv, poly, v, nit] = minimizePoly(s,p,h,ev_space,eta,tol_feasible,max_iter,verbose=True)
    return [poly, h]
  else:
    return [None, None]


if True:
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
    beta =0.1;
#    beta = 10;

    imag_lim = beta;
    l1 = 1j*np.linspace(0,imag_lim,50);
    l2 = 1j*imag_lim + np.linspace(-kappa,0,npts);
    l3 = -kappa + 1j*np.linspace(0,imag_lim,50);
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




if False:
  npts = 1000;
  ev_space = -np.linspace(0,1,npts);
  h = 3.165161132812499e+01;
  s = 10;
  h= 1.613250732421875e+03*0.99;
  s = 100;
  p = 8;
  [conv, poly, v, nit] = minimizePoly(s,p,h,ev_space,maxiter=128,verbose=True)

  if conv:
    print "------ Polynomial coefficients -------- "
    for x in poly:
      sys.stdout.write("%.16g," % x)
    print ""
  else:
    print " ------- Not converged ------ "
