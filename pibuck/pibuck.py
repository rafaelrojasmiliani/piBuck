


import sympy as sp
import numpy as np
from . import tools


class cPiBuck:
  def __init__(self):
    self.cnst=[];
    self.cnstDim=[];#fundim of the variables
    self.var=[];#variables
    self.varDim=[];
    self.nci=[];
    self.dims=[];#fundamental fundim in the system
    M,L,T,I,Th=sp.symbols('M L T I \Theta')
    self.fundim=[M,L,T,I,Th]
    self.fundimdic={ 'mass':M,'length':L,'time':T,\
                'current':I,'temperature':Th,\
                'force':M*L/T**2,\
                'power':(M)*(L**2)*(T**-3),\
                'energy':(M)*(L**2)*(T**-2),\
                'density':(M)*(L**3),\
                'electrical resistance':M*(L**2)*(T**-3)*(I**-2),\
                'electrical conductance':(M**-2)*(L**-2)*(T**3)*(I**2),\
                'electrical capacitance':(M**-1)*(L**-2)*(T**4)*(I**2),\
                'electric potential':(M)*(L**2)*(T**-3)*(I**-1)}
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
    print 'problem constants: %i'%len(self.cnst);
    for i,j in zip(self.cnst,self.cnstDim):
      print '['+str(i)+'] = '+str(j)+' ',
    
    print '\nproblem variables: %i'%len(self.var);
    for i,j in zip(self.var,self.varDim):
      print '['+str(i)+'] = '+str(j)+' ',
    
    s=''
    for i in self.dims:
      s=s+' '+str(i)
    print '\nproblem physical magnitudes '+ s
    print 'by the pi Buckingham theorem we have '+\
    str(len(self.cnst)+len(self.var))+' - '+\
    str(len(self.dims))+' = '+str(len(self.cnst)+len(self.var)-len(self.dims))+' variables'
    dm=self.getDimMat();
    minpars=int(dm.shape[1]-tools.rank(dm));
    print 'you need '+str(minpars)+' parameters to construct a dimensionless sytem'
    
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
      for j,const in enumerate(self.cnstDim):
        for i,dim in enumerate(self.dims):
          mc[i,j]=tools.powInMonomial(dim,const);
      mv=np.zeros((len(self.dims),len(self.var)))
      for j,const in enumerate(self.cnstDim):
        for i,dim in enumerate(self.dims):
          mv[i,j]=tools.powInMonomial(dim,const);
      return np.concatenate((mc,mv),axis=1);
  def getPiPars(self,v):
    dm=self.getDimMat();
    minpars=int(dm.shape[1]-tools.rank(dm));



    
