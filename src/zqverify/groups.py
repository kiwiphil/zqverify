"""Exact finite Abelian group and quotient helpers for the detector paper."""
from __future__ import annotations
from dataclasses import dataclass
from collections import deque
from heapq import heappush, heappop
from itertools import product
from math import inf
from typing import Iterable, Sequence
from functools import lru_cache

Element = tuple[int, ...]

@dataclass(frozen=True)
class FiniteAbelianGroup:
    moduli: tuple[int, ...]
    def __post_init__(self):
        if not self.moduli or any(n < 2 for n in self.moduli):
            raise ValueError("moduli must all be >=2")
    @property
    def zero(self) -> Element:
        return (0,) * len(self.moduli)
    @property
    def elements(self) -> tuple[Element, ...]:
        return tuple(product(*(range(n) for n in self.moduli)))
    @property
    def order(self) -> int:
        x=1
        for n in self.moduli: x*=n
        return x
    def norm(self, a: Sequence[int]) -> Element:
        if len(a)!=len(self.moduli): raise ValueError("wrong element length")
        return tuple(int(x)%n for x,n in zip(a,self.moduli))
    def add(self,a:Sequence[int],b:Sequence[int])->Element:
        return tuple((int(x)+int(y))%n for x,y,n in zip(a,b,self.moduli))
    def neg(self,a:Sequence[int])->Element:
        return tuple((-int(x))%n for x,n in zip(a,self.moduli))
    def sub(self,a,b)->Element:
        return self.add(a,self.neg(b))
    def smul(self,k:int,a:Sequence[int])->Element:
        return tuple((int(k)*int(x))%n for x,n in zip(a,self.moduli))


def subgroup_generated(group:FiniteAbelianGroup, generators:Iterable[Sequence[int]])->frozenset[Element]:
    gens=tuple(group.norm(g) for g in generators)
    H={group.zero}
    changed=True
    while changed:
        changed=False
        for h in tuple(H):
            for g in gens:
                x=group.add(h,g)
                if x not in H:
                    H.add(x); changed=True
                x=group.sub(h,g)
                if x not in H:
                    H.add(x); changed=True
    return frozenset(H)


@lru_cache(maxsize=None)
def enumerate_subgroups(group:FiniteAbelianGroup)->tuple[frozenset[Element], ...]:
    seen={frozenset((group.zero,))}
    queue=deque(seen)
    elems=group.elements
    while queue:
        H=queue.popleft()
        for g in elems:
            if g in H: continue
            K=subgroup_generated(group, tuple(H)+(g,))
            if K not in seen:
                seen.add(K); queue.append(K)
    return tuple(sorted(seen,key=lambda H:(len(H),tuple(sorted(H)))))


def in_subgroup(group:FiniteAbelianGroup,a:Sequence[int],N:frozenset[Element])->bool:
    return group.norm(a) in N


def quotient_key(group:FiniteAbelianGroup,a:Sequence[int],N:frozenset[Element])->Element:
    a=group.norm(a)
    return min(group.add(a,n) for n in N)


def quotient_equal(group,a,b,N)->bool:
    return group.sub(a,b) in N


def quotient_elements(group:FiniteAbelianGroup,N:frozenset[Element])->tuple[Element,...]:
    return tuple(sorted({quotient_key(group,a,N) for a in group.elements}))


def quotient_is_cyclic(group:FiniteAbelianGroup,N:frozenset[Element])->bool:
    qels=quotient_elements(group,N)
    qsize=len(qels)
    if qsize<=1: return True
    z=quotient_key(group,group.zero,N)
    for g in group.elements:
        seen={z}; cur=z
        for _ in range(qsize):
            cur=quotient_key(group,group.add(cur,g),N)
            seen.add(cur)
        if len(seen)==qsize: return True
    return False


def signed_word_length(group, active_colors, target, N=frozenset()):
    if not N: N=frozenset((group.zero,))
    target=quotient_key(group,target,N); z=quotient_key(group,group.zero,N)
    if target==z: return 0
    steps=[]
    for c in active_colors:
        c=group.norm(c)
        if c in N: continue
        steps.extend((c,group.neg(c)))
    dist={z:0}; q=deque((z,))
    while q:
        x=q.popleft(); d=dist[x]
        for step in steps:
            y=quotient_key(group,group.add(x,step),N)
            if y not in dist:
                if y==target: return d+1
                dist[y]=d+1; q.append(y)
    raise ValueError("visible colors do not generate target quotient")


def weighted_word_length(group, active_colors, costs, target, N=frozenset()):
    if not N: N=frozenset((group.zero,))
    target=quotient_key(group,target,N); z=quotient_key(group,group.zero,N)
    if target==z: return 0.0
    steps=[]
    for c,u in zip(active_colors,costs):
        c=group.norm(c)
        if c in N: continue
        steps.extend(((c,float(u)),(group.neg(c),float(u))))
    dist={z:0.0}; heap=[(0.0,z)]
    while heap:
        d,x=heappop(heap)
        if d!=dist[x]: continue
        if x==target: return d
        for step,u in steps:
            y=quotient_key(group,group.add(x,step),N); nd=d+u
            if nd < dist.get(y,inf)-1e-15:
                dist[y]=nd; heappush(heap,(nd,y))
    raise ValueError("visible colors do not generate target quotient")


def canonical_signed_representatives(group:FiniteAbelianGroup)->tuple[Element,...]:
    """One representative from each nonzero inverse pair {x,-x}."""
    out=[]; used={group.zero}
    for x in group.elements:
        if x in used: continue
        nx=group.neg(x)
        r=min(x,nx)
        out.append(r); used.add(x); used.add(nx)
    return tuple(sorted(set(out)))
