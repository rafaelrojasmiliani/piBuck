


"""Core functionality for computing dimensionless groups."""

from __future__ import annotations

import numpy as np
import sympy as sp

from . import tools


class cPiBuck:
    """Utility to obtain dimensionless parameters using Buckingham's Pi."""

    def __init__(self):
        self.cnst = []
        self.cnstDim = []  # fundamental dimensions of the constants
        self.var = []
        self.varDim = []
        self.nci = []
        self.dims = []  # fundamental dimensions in the system

        M, L, T, I, Th = sp.symbols("M L T I \\Theta")
        self.fundim = [M, L, T, I, Th]
        self.fundimdic = {
            "mass": M,
            "length": L,
            "time": T,
            "electric current": I,
            "temperature": Th,
            "force": M * L / T ** 2,
            "power": (M) * (L ** 2) * (T ** -3),
            "energy": (M) * (L ** 2) * (T ** -2),
            "density": (M) * (L ** 3),
            "electric resistance": M * (L ** 2) * (T ** -3) * (I ** -2),
            "electric conductance": (M ** -1) * (L ** -2) * (T ** 3) * (I ** 2),
            "electric capacitance": (M ** -1) * (L ** -2) * (T ** 4) * (I ** 2),
            "electric potential": (M) * (L ** 2) * (T ** -3) * (I ** -1),
            "electric charge": (T) * (I),
            "electric inductance": (M) * (L ** 2) * (T ** -2) * (I ** -2),
            "dielectric coefficient": (I ** 2) * (T ** 4) * (M ** -1) * (L ** -3),
        }
    def getFunDimArray(self):
        return self.fundim

    def getFunDimDic(self):
        return self.fundimdic

    # add variable
    def addVar(self, v, d):
        self.var.append(v)
        self.varDim.append(d)
        for dim in self.fundim:
            if dim not in self.dims and tools.powInMonomial(dim, d) != 0:
                self.dims.append(dim)

    # add constant
    def addConst(self, c, d):
        self.cnst.append(c)
        self.cnstDim.append(d)
        for dim in self.fundim:
            if dim not in self.dims and tools.powInMonomial(dim, d) != 0:
                self.dims.append(dim)

    def info(self):
        """Print a short summary about the problem."""

        print(f"problem constants: {len(self.cnst)}")
        for name, dim in zip(self.cnst, self.cnstDim):
            print(f"[{name}] = {dim} ", end="")

        print(f"\nproblem variables: {len(self.var)}")
        for name, dim in zip(self.var, self.varDim):
            print(f"[{name}] = {dim} ", end="")

        s = "".join(f" {i}" for i in self.dims)
        print("\nproblem physical magnitudes" + s)
        print(
            "by the pi Buckingham theorem we have "
            + f"{len(self.cnst)+len(self.var)} - {len(self.dims)} = "
            + f"{len(self.cnst)+len(self.var)-len(self.dims)} variables"
        )
        dm = self.getDimMat()
        print(
            f"you need {tools.rank(dm)} parameters to construct a dimensionless sytem"
        )

    def setNormalizingConst(self, arr):
        for i, c in enumerate(self.cnst):
            if c in arr:
                self.nci.append(i)
            else:
                self.nnci.append(i)

    def getDimMat(self):
        mc = np.zeros((len(self.dims), len(self.cnst)))
        for j, const in enumerate(self.cnstDim):
            for i, dim in enumerate(self.dims):
                mc[i, j] = tools.powInMonomial(dim, const)
        mv = np.zeros((len(self.dims), len(self.var)))
        for j, const in enumerate(self.varDim):
            for i, dim in enumerate(self.dims):
                mv[i, j] = tools.powInMonomial(dim, const)
        return np.concatenate((mc, mv), axis=1)

    def getPiPars(self, v, subsDict):
        result = {}
        dm = self.getDimMat()
        A = np.zeros((len(self.dims), len(v)))
        j1 = 0
        mask = []

        allPars = self.cnst[:] + self.var
        for j, dimPar in enumerate(allPars):
            if dimPar in v:
                A[:, j1] = dm[:, j]
                j1 += 1
            else:
                mask.append(j)

        for j in mask:
            p = subsDict[allPars[j]]
            auxVec = np.linalg.solve(A, dm[:, j])
            for k, s in enumerate(auxVec):
                p = p * (v[k] ** sp.Rational(s))
            result[allPars[j]] = p
        return result

    
    def normalizeExp(self, exp, t, t0, v, subsDict):
        result = exp
        pipars = self.getPiPars(v, subsDict)
        dt0dt = sp.simplify(t0 / pipars[t])
        for par, pip in pipars.items():
            if result.has(par) and (par in self.var) and par != t:
                result = tools.subsVars(result, {par: pip}, [t, t0], dt0dt)

        for par, pip in pipars.items():
            if result.has(par) and (par in self.cnst):
                result = result.subs(par, pip)
        return result


    
