from pibuck import *
from sympy import *
#define equations of motion

f=[0,0,0];

t=symbols(r't',real=True);
z=symbols(r'z',cls=Function)(t);
F=symbols(r'F',cls=Function)(t)
m,c,k=symbols(r'm c k',positive=True)

eq=m*diff(z,t,2)+c*diff(z,t)+k*z==F
pb=cPiBuck();

fd=pb.getFunDimDic();

pb.addVar(z,fd['length'])
pb.addVar(t,fd['time'])
pb.addConst(m,fd['mass'])
pb.addConst(c,fd['mass']/fd['time'])
pb.addConst(k,fd['mass']/fd['time']**2)
pb.addConst(F,fd['force'])

pb.info()

print(pb.getDimMat())

res=pb.getPiPars([m,k,F]);
print(res)
