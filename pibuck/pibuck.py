


import sympy as sp
import numpy as np
from . import tools


class cPiBuck:
  def __init__(self):
    self.cnst=[]
    self.cnstDim=[]  # fundim of the variables
    self.var=[]  # variables
    self.varDim=[]
    self.nci=[]
    self.nnci=[]
    self.dims=[]  # fundamental fundim in the system
    M,L,T,I,Th=sp.symbols('M L T I \Theta')
    self.fundim=[M,L,T,I,Th]
    self.fundimdic={ 'mass':M,'length':L,'time':T,\
                'electric current':I,'temperature':Th,\
                'force':M*L/T**2,\
                'power':(M)*(L**2)*(T**-3),\
                'energy':(M)*(L**2)*(T**-2),\
                'density':M*(L**-3),\
                'electric resistance':M*(L**2)*(T**-3)*(I**-2),\
                'electric conductance':(M**-1)*(L**-2)*(T**3)*(I**2),\
                'electric capacitance':(M**-1)*(L**-2)*(T**4)*(I**2),\
                'electric potential':(M)*(L**2)*(T**-3)*(I**-1),\
                'electric charge':(T)*(I),\
                'electric inductance':(M)*(L**2)*(T**-2)*(I**-2),\
                'dielectric coefficient':(I**2)*(T**4)*(M**-1)*(L**-3)}
  def getFunDimArray(self):
    return self.fundim;
  def getFunDimDic(self):
    return self.fundimdic;
  # add variable
  def addVar(self,v,d):
    self.var.append(v);
    self.varDim.append(d);
    for dim in self.fundim:
      if(not(dim in self.dims) and tools.powInMonomial(dim,d)!=0):
        self.dims.append(dim);
  #add constant
  def addConst(self,c,d):
    self.cnst.append(c);
    self.cnstDim.append(d);
    for dim in self.fundim:
      if(not(dim in self.dims) and tools.powInMonomial(dim,d)!=0):
        self.dims.append(dim);
  #initialize the problem, tells the user what he can 
  #do
  def info(self):
    print('problem constants: %i' % len(self.cnst))
    for i, j in zip(self.cnst, self.cnstDim):
      print('[' + str(i) + '] = ' + str(j) + ' ', end='')

    print('\nproblem variables: %i' % len(self.var))
    for i, j in zip(self.var, self.varDim):
      print('[' + str(i) + '] = ' + str(j) + ' ', end='')

    s=''
    for i in self.dims:
      s=s+' '+str(i)
    print('\nproblem physical magnitudes '+ s)
    print('by the pi Buckingham theorem we have ' +
    str(len(self.cnst)+len(self.var))+' - '+
    str(len(self.dims))+' = '+str(len(self.cnst)+len(self.var)-len(self.dims))+' variables')
    dm=self.getDimMat();
    print('you need '+str(tools.rank(dm))+' parameters to construct a dimensionless sytem')
    
  #tell to the class which variables you whish to
  #use to make the system dimensionless
  def setNormalizingConst(self,a):
    for i,c in enumerate(self.cnst):
        if(c in a):
          self.nci.append(i);
        else:
          self.nnci.append(i);
  def getDimMat(self):
      mc=np.zeros((len(self.dims),len(self.cnst)))
# for on the dimension of cosntatns
      for j,const in enumerate(self.cnstDim):
        for i,dim in enumerate(self.dims):
          mc[i,j]=tools.powInMonomial(dim,const);
      mv=np.zeros((len(self.dims),len(self.var)))
# for on the dimension of variables
      for j,const in enumerate(self.varDim):
        for i,dim in enumerate(self.dims):
          mv[i,j]=tools.powInMonomial(dim,const);
      return np.concatenate((mc,mv),axis=1);

  def getPiPars(self,v,subsDict):
    result={}
    dm=self.getDimMat();
    minpars=int(dm.shape[1]-tools.rank(dm));
    A=np.zeros((len(self.dims),len(v)));
    j1=0
    mask=[]
# a for in the const that access to dm, this is
# possible because constants are the first elements of dm
    allPars=self.cnst[:];
    allPars=allPars+self.var
    #print 'aux = '+str(aux)
#dimPar= dimentional Paramenter
    for j,dimPar in enumerate(allPars):
      if(dimPar in v):
        A[:,j1]=dm[:,j];
        j1=j1+1;
      else:
        mask.append(j)
# masl contain all the indexes of the
# dimentional parameters that we can
# substitute which are present in allPars array
    auxVec=0;
    j1=0;
    #print A
    for j in mask:
      p=subsDict[allPars[j]]
      auxVec=np.linalg.solve(A,dm[:,j])
      for k,s in enumerate(auxVec):
        p=p*(v[k]**sp.Rational(s))
      result[allPars[j]]=p
    return result;

    
  def normalizeExp(self,exp,t,t0,v,subsDict):
    result=exp;
    pipars=self.getPiPars(v,subsDict)
    dt0dt=sp.simplify(t0/pipars[t]);
    for par,pip in pipars.items():
      if result.has(par) and (par in self.var) and par!=t:
        result=tools.subsVars(\
          result,{par:pip},[t,t0],dt0dt)

    for par,pip in pipars.items():
      if result.has(par) and (par in self.cnst):
        result=result.subs(par,pip)
    return result;


    
