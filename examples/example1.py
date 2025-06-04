# Example 1: piezoelectric energy harvester
from pibuck import *
from sympy import *

# Define the variables of the problem
# piezoelectric energy harvester

t = symbols(r't', real=True)      # time
z = symbols(r'z', cls=Function)(t)  # sprung mass displacement
y = symbols(r'y', cls=Function)(t)  # base mass displacement
V = symbols(r'V', cls=Function)(t)  # circuit voltage
R = symbols(r'R', cls=Function)(t)  # circuit resistance
p = symbols(r'p', cls=Function)(t)  # harvested power

m = symbols(r'm', positive=True)
c = symbols(r'c', positive=True)
k = symbols(r'k', positive=True)
yinf = symbols(r'y_inf', positive=True)
kp = symbols(r'k_p', positive=True)
Cp = symbols(r'Cp', positive=True)

# System equations

eq1 = Eq(m*diff(z, t, 2) + c*diff(z, t) + k*z, -m*diff(y, t, 2) - kp*V)
eq2 = Eq(Cp*diff(V, t) + V/R, kp*diff(z, t))
eq3 = Eq(p, -V**2 / R)

pb = cPiBuck()
fd = pb.getFunDimDic()

pb.addVar(t, fd['time'])
pb.addVar(z, fd['length'])
pb.addVar(y, fd['length'])
pb.addVar(V, fd['electric potential'])
pb.addVar(R, fd['electric resistance'])
pb.addVar(p, fd['power'])

pb.addConst(m, fd['mass'])
pb.addConst(c, fd['mass']/fd['time'])
pb.addConst(k, fd['mass']/fd['time']**2)
pb.addConst(yinf, fd['length'])
pb.addConst(kp, fd['force']*(fd['electric potential']**-1))
pb.addConst(Cp, fd['electric capacitance'])

pb.info()

# Dimensionless variables
zeta = symbols(r'zeta', positive=True)
cp = symbols(r'cp', positive=True)
tau = symbols(r'tau', positive=True)
xi = symbols(r'xi', cls=Function)(tau)
eta = symbols(r'eta', cls=Function)(tau)
v = symbols(r'v', cls=Function)(tau)
r = symbols(r'r', cls=Function)(tau)
pi = symbols(r'pi', cls=Function)(tau)

subsDict = {
    t: tau,
    z: xi,
    y: eta,
    V: v,
    R: r,
    p: pi,
    c: zeta,
    Cp: cp,
}

adimMapPars = [m, k, yinf, kp]
res = pb.getPiPars(adimMapPars, subsDict)
print("\n---the relation between the dimensional variables and the pi pars is---\n")
for Vv, j in res.items():
    print(Vv, " = ", j)

print("\n---The original expressions are---\n")
eq1bar = pb.normalizeExp(eq1, t, tau, adimMapPars, subsDict)
eq2bar = pb.normalizeExp(eq2, t, tau, adimMapPars, subsDict)

print(simplify(eq1bar))
print(simplify(eq2bar))
