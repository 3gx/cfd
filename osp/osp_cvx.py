from cvxpy import *
import numpy as np
from scipy import misc

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

  return [b.T,c.T]


npts = 100
ev_space = -np.linspace(0,1,npts);
h = 200; #3.165161132812499e+01;
s = 20;
p=1;

hval = h*ev_space;
[B,C] = scaled_chebyshev_basis(s,p,min(np.real(hval)),0,hval)

fixed_coefficients = np.ones(p+1)/misc.factorial(np.linspace(0,p,p+1))

x = Variable(s+1)
objective  = Minimize(max_entries(abs(C*x)-1))
constraints = [B*x-fixed_coefficients == 0]
prob = Problem(objective, constraints)
print "Optimal value", prob.solve(solver=CVXOPT,verbose=True, reltol=1e-15,abstol=1e-15, refinement=4, kktsolver)
#print "Optimal var"
#print x.value

poly = np.array(x.value).T
print poly[0]
print np.dot(B,poly[0])
