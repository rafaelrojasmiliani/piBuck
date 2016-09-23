
import sympy as sp

def powInMonomial(s,d):
  if(isinstance(d,sp.Symbol)):
    print 'is symbol',d
    if(d==s):
      return 1
    else:
      return 0
  elif(isinstance(d,sp.Pow)):
    print 'is pow',d.args
    if(s==d.args[0]):
      return d.args[1]
    else:
      return 0
  elif(isinstance(d,sp.Mul)):
    print d.args
    for m in d.args:
      print m
      res=powInMonomial(s,m);
      if(res!=0):
        return res
  return 0 
