"""Detector scores and finite optimization."""
from __future__ import annotations
from dataclasses import dataclass
from math import log
from functools import lru_cache
from .groups import *

@dataclass(frozen=True)
class DetectorResult:
    kernel: frozenset[Element]
    rho: float
    word_length: int
    coefficient: float
    quotient_cyclic: bool


@lru_cache(maxsize=None)
def _ell_cyclic(group, active_colors, psi, N):
    return signed_word_length(group,active_colors,psi,N), quotient_is_cyclic(group,N)

def detector_result(group, active_colors, effective_activities, psi, N):
    psi=group.norm(psi)
    if psi in N: return None
    visible=[i for i,c in enumerate(active_colors) if group.norm(c) not in N]
    rho=sum(float(effective_activities[i]) for i in visible)
    if not visible or not (rho < 1.0): return None
    ell,cyc=_ell_cyclic(group,tuple(active_colors),psi,N)
    return DetectorResult(N,rho,ell,ell*log(1.0/rho),cyc)


def optimize_detectors(group, active_colors, effective_activities, psi):
    vals=[]
    for N in enumerate_subgroups(group):
        r=detector_result(group,active_colors,effective_activities,psi,N)
        if r is not None: vals.append(r)
    vals.sort(key=lambda r:(r.coefficient,-r.rho,r.word_length),reverse=True)
    return tuple(vals)


def best_detector(group, active_colors, effective_activities, psi, cyclic_only=False):
    vals=optimize_detectors(group,active_colors,effective_activities,psi)
    vals=tuple(r for r in vals if (r.quotient_cyclic or not cyclic_only))
    return vals[0] if vals else None
