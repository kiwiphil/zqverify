"""Signed multicolour currents and the v6 deterministic quotient extraction."""
from __future__ import annotations
from dataclasses import dataclass
from math import factorial
from typing import Sequence
from .groups import FiniteAbelianGroup, Element, quotient_equal
from .lattice import OrientedBox

@dataclass(frozen=True, order=True)
class CurrentType:
    plaquette: int
    color: int
    eps: int

@dataclass(frozen=True)
class MulticolorProblem:
    box: OrientedBox
    group: FiniteAbelianGroup
    active_colors: tuple[Element,...]
    def __post_init__(self):
        object.__setattr__(self,'active_colors',tuple(self.group.norm(c) for c in self.active_colors))
    @property
    def types(self):
        # +1 before -1, deterministic within plaquette/color.
        return tuple(CurrentType(p,c,e) for p in range(len(self.box.plaquettes)) for c in range(len(self.active_colors)) for e in (1,-1))

@dataclass(frozen=True)
class Current:
    counts: tuple[int,...]
    def __post_init__(self):
        if any(n<0 for n in self.counts): raise ValueError
    @property
    def size(self): return sum(self.counts)
    def leq(self,other): return all(a<=b for a,b in zip(self.counts,other.counts))
    def sub(self,other):
        if not other.leq(self): raise ValueError
        return Current(tuple(a-b for a,b in zip(self.counts,other.counts)))


def zero_current(problem): return Current((0,)*len(problem.types))

def add_one(current,i):
    x=list(current.counts); x[i]+=1; return Current(tuple(x))

def flux(problem,current):
    g=problem.group; out=[g.zero for _ in problem.box.plaquettes]
    for idx,(tau,n) in enumerate(zip(problem.types,current.counts)):
        if not n: continue
        out[tau.plaquette]=g.add(out[tau.plaquette], g.smul(tau.eps*n,problem.active_colors[tau.color]))
    return tuple(out)

def source(problem,current):
    g=problem.group; fl=flux(problem,current); out=[]
    B=problem.box.boundary2_matrix
    for e,row in enumerate(B):
        s=g.zero
        for p,inc in enumerate(row):
            if inc: s=g.add(s,g.smul(-inc,fl[p]))
        out.append(s)
    return tuple(out)

def source_is_cycle(problem,j):
    g=problem.group; B1=problem.box.boundary1_matrix
    for row in B1:
        s=g.zero
        for e,inc in enumerate(row):
            if inc: s=g.add(s,g.smul(inc,j[e]))
        if s!=g.zero: return False
    return True

def current_weight(problem,current,betas):
    out=1.0
    for tau,n in zip(problem.types,current.counts):
        if n: out *= (float(betas[tau.color])/2.0)**n/factorial(n)
    return out

def incidence_contribution(problem,tau,edge_index):
    inc=problem.box.boundary2_matrix[edge_index][tau.plaquette]
    return problem.group.smul(-tau.eps*inc,problem.active_colors[tau.color])

@dataclass(frozen=True)
class ExtractionStep:
    edge:int; type_index:int; defect_before:Element; defect_after:Element
@dataclass(frozen=True)
class ExtractionResult:
    extracted:Current; remainder:Current; steps:tuple[ExtractionStep,...]


def extract(problem,current,target,N):
    g=problem.group
    if len(current.counts)!=len(problem.types): raise ValueError
    q=zero_current(problem); steps=[]
    while True:
        sq=source(problem,q)
        defective=None
        for e,(a,b) in enumerate(zip(sq,target)):
            if not quotient_equal(g,a,b,N): defective=e; break
        if defective is None: break
        e=defective
        defect_before=g.sub(target[e],sq[e])
        chosen=None
        for ti,tau in enumerate(problem.types):
            if q.counts[ti] >= current.counts[ti]: continue
            if problem.box.boundary2_matrix[e][tau.plaquette]==0: continue
            c=incidence_contribution(problem,tau,e)
            if c in N: continue
            chosen=ti; break
        if chosen is None: raise AssertionError("availability lemma failed")
        q2=add_one(q,chosen)
        defect_after=g.sub(target[e],source(problem,q2)[e])
        steps.append(ExtractionStep(e,chosen,defect_before,defect_after)); q=q2
        if len(steps)>sum(current.counts): raise AssertionError("termination bound failed")
    r=current.sub(q)
    sr=source(problem,r)
    if not all(x in N for x in sr): raise AssertionError("remainder source not N-valued")
    return ExtractionResult(q,r,tuple(steps))
