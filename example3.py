

from pibuck import *
from sympy import *
#define equations of motion
# for the electromagnetic energy harvester

f=[0,0,0];

t=symbols(r't',real=True);
z=symbols(r'z',cls=Function)(t);
y=symbols(r'y',cls=Function)(t)
Q=symbols(r'Q',cls=Function)(t)
R1=symbols(r'R_1',cls=Function)(t)
R2=symbols(r'R_2',cls=Function)(t)
p=symbols(r'p',cls=Function)(t)
m,c,k,yinf,epsilon,A,d,V0=symbols(r'm c k y_inf epsilon A d V_0',positive=True)

eq1=Eq(m*diff(z,t,2)+c*diff(z,t)+k*z,-m*diff(y,t,2)+Rational(1,2)*Q**2/(epsilon*A))
eq2=Eq(diff(Q,t),-1/R1*(Q*(d-z)/(epsilon*A)-V0)-1/R2*(Q*(d-z)/(epsilon*A)))

pb=cPiBuck();

fd=pb.getFunDimDic();

pb.addVar(t,fd['time'])
pb.addVar(z,fd['length'])
pb.addVar(y,fd['length'])
pb.addVar(Q,fd['electric charge'])
pb.addVar(R1,fd['electric resistance'])
pb.addVar(R2,fd['electric resistance'])
pb.addConst(m,fd['mass'])
pb.addConst(c,fd['mass']/fd['time'])
pb.addConst(k,fd['mass']/fd['time']**2)
pb.addConst(yinf,fd['length'])
pb.addConst(A,fd['length']**2)
pb.addConst(d,fd['length'])
pb.addConst(V0,fd['electric potential'])
pb.addConst(epsilon,fd['dielectric coefficient'])

pb.addVar(p,fd['power'])

pb.info()

#print pb.getDimMat()
#
tau=symbols(r'tau',positive=True)
xi =symbols(r'xi',cls=Function)(tau)
eta=symbols(r'eta',cls=Function)(tau)
q=symbols(r'q',cls=Function)(tau)
r1=symbols(r'r_1',cls=Function)(tau)
r2=symbols(r'r_2',cls=Function)(tau)

zeta=symbols(r'zeta',positive=True)
v0=symbols(r'v_0',positive=True)
d0=symbols(r'd_0',positive=True)
A0=symbols(r'A_0',positive=True)
epsilon0=symbols(r'epsilon_0',positive=True)


pi=symbols(r'pi',cls=Function)(tau)

subsDict=   {c:zeta,t:tau,z:xi,y:eta,Q:q,\
             R1:r1,R2:r2,epsilon:epsilon0,d:d0,A:A0,p:pi}

adimMapPars=[m,k,yinf,V0]
res=pb.getPiPars(adimMapPars,subsDict)
eq1bar=pb.normalizeExp(eq1,t,tau,adimMapPars,subsDict);
eq2bar=pb.normalizeExp(eq2,t,tau,adimMapPars,subsDict);

u1,u2=symbols(r'u_1 u_2',positive=True)
eq2bar=eq2bar.subs({r1:1/u1,r2:1/u2})

eq1dl=tools.simplifyEq(eq1bar)
eq2dl=tools.simplifyEq(eq2bar)


