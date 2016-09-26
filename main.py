from pibuck import *
from sympy import *
#define equations of motion

f=[0,0,0];

v,z,t=symbols(r'v z t',real=True);
F=symbols(r'F',real=True)
m,c,k=symbols(r'm c k',positive=True)

f[0]=v
f[1]=1/m*(-k*z-c*v+F)

pb=cPiBuck();

fd=pb.getFunDimDic();

pb.addVar(v,fd['length']/fd['time'])
pb.addVar(z,fd['length'])
pb.addVar(F,fd['force'])
pb.addVar(t,fd['time'])
pb.addConst(m,fd['mass'])
pb.addConst(c,fd['mass']/fd['time'])
pb.addConst(k,fd['mass']/fd['time']**2)

pb.info()

print pb.getDimMat()
