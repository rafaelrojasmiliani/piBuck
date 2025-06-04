#M=Dimension(mass=1,symbol='M')
#T=Dimension(time=1,symbol='T')
#L=Dimension(length=1,symbol='L')
#I=Dimension(current=1,symbol='I')
#Th=Dimension(current=1,symbol='\Theta')

from sympy import *
import numpy as np
import scipy

M,L,T,I,Th=symbols('M L T I \Theta')
dimensions=[M,L,T,I,Th]
dimNum=4

m=M
c=M/T
k=M/T**2
y0=L/T**2
Cp=T**4*I**2/L**2/M
kp=T*I/L

V=M*L**2/T**3/I
x=L
t=T
u=I**2*T**3/M/L**2

constants=[m,c,k,y0,Cp,kp];
adc0=[0,2,3,4]
adc1=[1,5];
variables=[x,t,V,u]

#construc dimensional matrix for the constant:
def getDimMat(vec):
    mc=np.zeros((dimNum,len(vec)))
    for j,const in enumerate(vec):
        if(isinstance(const,Symbol)):
            mc[dimensions.index(const),j]=1
        elif(isinstance(const,Pow)):
            mc[dimensions.index(const.args[0]),j]=const.args[1]
        elif(isinstance(const,Mul)):
            print(const, ' is a mul')
            for d in const.args:
                if(isinstance(d,Symbol)):
                    mc[dimensions.index(d),j]=1
                elif(isinstance(d,Pow)):
                    mc[dimensions.index(d.args[0]),j]=d.args[1]
                else:
                    print('error in', j, const, d, type(d))
        else:
            print('error in ', j, const)
    return mc
                
def binomialCoeff(n, k):
    result = 1
    for i in range(1, k+1):
        result = result * (n-i+1) / i
    return int(result)
            
mc=getDimMat(constants);
mv=getDimMat(variables);
md=mc[:,adc0];
mc=mc[:,adc1];


if(rank(mv)!=dimNum):
    print("the problem can't be solved")



