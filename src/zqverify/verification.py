"""Finite theorem-facing checks for projected filling and history levels."""
from __future__ import annotations
from collections import deque, defaultdict
from heapq import heappush, heappop
from math import inf
from .groups import quotient_key, quotient_equal, weighted_word_length
from .currents import incidence_contribution

def _projected_state(problem, edge_values, N):
    return tuple(quotient_key(problem.group,x,N) for x in edge_values)

def rectangle_target(problem, origin, axes, sizes, psi):
    gamma=problem.box.rectangular_loop(origin,axes,sizes)
    return tuple(problem.group.smul(v,psi) for v in gamma)

def minimum_projected_filling_size(problem,target,N,max_states=500000):
    g=problem.group; z=_projected_state(problem,(g.zero,)*len(problem.box.edges),N)
    tgt=_projected_state(problem,target,N)
    if z==tgt: return 0
    moves=[]
    for tau in problem.types:
        if problem.active_colors[tau.color] in N: continue
        vals=tuple(incidence_contribution(problem,tau,e) for e in range(len(problem.box.edges)))
        moves.append(_projected_state(problem,vals,N))
    q=deque([z]); dist={z:0}
    while q:
        x=q.popleft(); d=dist[x]
        for mv in moves:
            y=tuple(quotient_key(g,g.add(a,b),N) for a,b in zip(x,mv))
            if y not in dist:
                if y==tgt: return d+1
                dist[y]=d+1; q.append(y)
                if len(dist)>max_states: raise RuntimeError('BFS state cap exceeded')
    raise ValueError('target not fillable')

def minimum_projected_filling_cost(problem,target,N,costs,max_states=500000):
    g=problem.group; z=_projected_state(problem,(g.zero,)*len(problem.box.edges),N); tgt=_projected_state(problem,target,N)
    moves=[]
    for tau in problem.types:
        if problem.active_colors[tau.color] in N: continue
        vals=tuple(incidence_contribution(problem,tau,e) for e in range(len(problem.box.edges)))
        moves.append((_projected_state(problem,vals,N),float(costs[tau.color])))
    dist={z:0.0}; heap=[(0.0,z)]
    while heap:
        d,x=heappop(heap)
        if d!=dist[x]: continue
        if x==tgt: return d
        for mv,c in moves:
            y=tuple(quotient_key(g,g.add(a,b),N) for a,b in zip(x,mv)); nd=d+c
            if nd < dist.get(y,inf)-1e-14:
                dist[y]=nd; heappush(heap,(nd,y))
                if len(dist)>max_states: raise RuntimeError('Dijkstra state cap exceeded')
    raise ValueError('target not fillable')

def history_level_activity(problem,target,N,betas,depth):
    """Exact multiplicative activity H_depth of the formal history tree."""
    g=problem.group
    start=_projected_state(problem,(g.zero,)*len(problem.box.edges),N)
    tgt=_projected_state(problem,target,N)
    masses={start:1.0}
    for _ in range(depth):
        nxt=defaultdict(float)
        for state,mass in masses.items():
            if state==tgt: continue
            e=next(i for i,(a,b) in enumerate(zip(state,tgt)) if a!=b)
            for ti,tau in enumerate(problem.types):
                if problem.box.boundary2_matrix[e][tau.plaquette]==0: continue
                c=incidence_contribution(problem,tau,e)
                if c in N: continue
                vals=list(state)
                for ee in range(len(vals)):
                    cc=incidence_contribution(problem,tau,ee)
                    vals[ee]=quotient_key(g,g.add(vals[ee],cc),N)
                nxt[tuple(vals)] += mass*(float(betas[tau.color])/2.0)
        masses=dict(nxt)
    return sum(masses.values())
