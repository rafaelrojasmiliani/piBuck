
import sympy as sp
import scipy as scp
import numpy as np

def powInMonomial(s,d):
  if(isinstance(d,sp.Symbol)):
#    print 'is symbol',d
    if(d==s):
      return 1
    else:
      return 0
  elif(isinstance(d,sp.Pow)):
#    print 'is pow',d.args
    if(s==d.args[0]):
      return d.args[1]
    else:
      return 0
  elif(isinstance(d,sp.Mul)):
    #print d.args
    for m in d.args:
     # print m
      res=powInMonomial(s,m);
      if(res!=0):
        return res
  return 0 

def null(A, eps=1e-15):
  u, s, vh = np.linalg.svd(A)
  null_mask = (s <= eps)
  null_space = np.compress(null_mask, vh, axis=0)
  return scp.transpose(null_space)

def rank(A, atol=1e-13, rtol=0):
  A = np.atleast_2d(A)
  s = np.linalg.svd(A, compute_uv=False)
  tol = max(atol, rtol * s[0])
  result = int((s >= tol).sum())
  return result
