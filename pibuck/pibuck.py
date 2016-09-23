


import sympy as sp
import numpy as np
from . import tools
import scipy


class cPiBuck:
  def __init__(self):
    self.cnst=[];
    self.cnstDim=[];#fundim of the variables
    self.var=[];#variables
    self.varDim=[];
    self.nci=[];
    self.dims=[];#fundamental fundim in the system
    M,L,T,I,Th=sp.symbols('M L T I \Theta')
    fundim=[M,L,T,I,Th]
  def getFunDimArray(self):
    return fundim;
  # add variable
  def addVar(self,v,d):
    self.var.append(v);
    self.varDim.append(d);
    for dim in fundim:
      if(not(dim in self.dims) and tools.powInMonomial(dim,d)!=0):
        self.dims.append(dim);
  #add constant
  def addConst(self,c,d):
    self.cnst.append(v);
    self.cnstDim.append(d);
    for dim in fundim:
      if(not(dim in self.dims) and powInMonomial(dim,d)!=0):
        self.dims.append(dim);
  #tell to the class which variables you whish to
  #use to make the system dimensionless
  def setNormalizingConst(self,a):
    for i,c in enumerate(self.cnst):
        if(c in a):
          self.nci.append(i);
        else:
          self.nnci.append(i);
#  def getDimMat(self):
#      mc=np.zeros((len(self.dims),len(self.cnst)))
#      for j,const in enumerate(self.cnstDim):
#        for i,dim in enumerate(self.dims):
#          mc(i,j)=powInMonomial(const);
#      mv=np.zeros((len(self.dims),len(self.var)))
#      for j,const in enumerate(self.cnstDim):
#        for i,dim in enumerate(self.dims):
#          mv(i,j)=powInMonomial(const);
#      md=mc[:,self.nci];
#      cd=self.cnst[self.nci]
#      mc=mc[:,self.nnci];
#      mdv

    
