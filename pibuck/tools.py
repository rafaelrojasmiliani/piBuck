
"""Utility helpers used by :mod:`pibuck`."""

from __future__ import annotations

import numpy as np
import scipy as scp
import sympy as sp


def powInMonomial(symbol, monomial):
    """Return the exponent of ``symbol`` within ``monomial``."""

    if isinstance(monomial, sp.Symbol):
        return 1 if monomial == symbol else 0
    if isinstance(monomial, sp.Pow):
        return monomial.args[1] if monomial.args[0] == symbol else 0
    if isinstance(monomial, sp.Mul):
        for arg in monomial.args:
            res = powInMonomial(symbol, arg)
            if res != 0:
                return res
    return 0

def null(A, eps=1e-15):
    """Return the null space of a matrix."""

    u, s, vh = np.linalg.svd(A)
    null_mask = s <= eps
    null_space = np.compress(null_mask, vh, axis=0)
    return scp.transpose(null_space)

def rank(A, atol=1e-13, rtol=0):
    """Compute the rank of a matrix with a tolerance."""

    A = np.atleast_2d(A)
    s = np.linalg.svd(A, compute_uv=False)
    tol = max(atol, rtol * s[0])
    result = int((s >= tol).sum())
    return result


def subsVars(exp, mapping, timesubs, dnewdold):
    """Substitute variables and their derivatives."""

    t, tau = timesubs
    res = exp
    for var, var0 in mapping.items():
        for order in reversed(range(0, 4)):
            if exp.has(sp.diff(var, t, order)):
                res = res.subs(
                    sp.diff(var, t, order),
                    sp.diff(var0, tau, order) * dnewdold ** order,
                )
    return res
    
def simplifyEq(expr):
    """Simplify both sides of an equality by cancelling common factors."""

    lhs = expr.lhs.factor()
    rhs = expr.rhs.factor()
    res_lhs = lhs
    res_rhs = rhs
    if type(lhs) == type(rhs):
        for el in lhs.args:
            if el in rhs.args:
                res_lhs = res_lhs / el
                res_rhs = res_rhs / el
    return sp.Eq(res_lhs, res_rhs)

