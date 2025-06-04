
#   Example of using PiBuck symbolic library
#   for transforming a set of equations 
#   into a dimensionless one.
#
#
from pibuck import *
from sympy import *

# FIRST we define the variables of our problem
# an piezoelectric energy harvester
t=symbols(r't',real=True)      # time
z=symbols(r'z',cls=Function)(t)# sprung mass displacement
y=symbols(r'y',cls=Function)(t)# basement mass displacement
V=symbols(r'V',cls=Function)(t)# Circuit voltage
R=symbols(r'R',cls=Function)(t)# Circuit resistance
p=symbols(r'p',cls=Function)(t)# Harveted Power

m=symbols(r'm',positive=True)# Sprung mass
c=symbols(r'c',positive=True)# Damping coeficient
k=symbols(r'k',positive=True)# Stiffness
yinf=symbols(r'y_inf',positive=True)# Upper bounf of the ground displacement
kp=symbols(r'k_p',positive=True)# Piezzoelectric coefficient
Cp=symbols(r'Cp',positive=True)# Inductance of the circuit

#SECOND we define out system equations
eq1=Eq(m*diff(z,t,2)+c*diff(z,t)+k*z,-m*diff(y,t,2)-kp*V)
eq2=Eq(Cp*diff(V,t)+V/R,kp*diff(z,t))
eq3=Eq(p,-V**2/R)

#THIRTH call the construct of cPiBuck class and get a new instance
pb=cPiBuck();


# 4 Get the dictionary of Fundamental Dimension
fd=pb.getFunDimDic();

# 5 Add the variable to our system
pb.addVar(t,fd['time'])
pb.addVar(z,fd['length'])
pb.addVar(y,fd['length'])
pb.addVar(V,fd['electric potential'])
pb.addVar(R,fd['electric resistance'])
pb.addVar(p,fd['power'])
# 5 Add the constant to out system. We differentiate
#   between variable and constant, because we will 
#   construct the dimensionless monomials with 
#   constants
pb.addConst(m,fd['mass'])
pb.addConst(c,fd['mass']/fd['time'])
pb.addConst(k,fd['mass']/fd['time']**2)
pb.addConst(yinf,fd['length'])
pb.addConst(kp,fd['force']*(fd['electric potential']**-1))
pb.addConst(Cp,fd['electric capacitance'])

pb.info()


# After getting pb.info() we known that
# we need only 8 variables for our
# dimensionless system. We know define
# our desired set of dimensionless variables:

zeta= symbols(r'zeta',positive=True)    # dimensionless damping
cp   = symbols(r'cp',positive=True)       # dimensionless inductance
tau = symbols(r'tau',positive=True)     # dimensionless Time
xi  = symbols(r'xi',cls=Function)(tau)  # dimensionless sprung mass displacement
eta = symbols(r'eta',cls=Function)(tau) # dimensionless basement displacement.
v   = symbols(r'v',cls=Function)(tau)   # dimensionless voltage
r   = symbols(r'r',cls=Function)(tau)   # dimensionless resistance
pi  = symbols(r'pi',cls=Function)(tau)  # dimensionless power

# Now we define the set of variable to
# substitute. We make a dictionary that
# associate to each dimensional variable
# its dimensionless substitute.
subsDict=   {t:tau,z:xi,y:eta,V:v,\
             R:r,p:pi,\
             c:zeta,Cp:cp}
# And to construct that set of dimensionless 
# variables we choose a set of constant.
adimMapPars=[m,k,yinf,kp]
#CAUTION: this set of constant must be sufficiently rich.
#         V.e., its dimensional matrix must span all the
#         space.
res=pb.getPiPars(adimMapPars,subsDict)
print("\n---the relation between the dimensional variables and the pi pars is---\n")
for V,j in res.items():
  print(V, " = ", j)

print("\n---The original expressions are---\n")
eq1bar=pb.normalizeExp(eq1,t,tau,adimMapPars,subsDict);
eq2bar=pb.normalizeExp(eq2,t,tau,adimMapPars,subsDict);

print(simplify(eq1bar))
print(simplify(eq2bar))
