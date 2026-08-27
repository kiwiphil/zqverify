from zqverify import *
from zqverify.currents import Current, MulticolorProblem, source, extract, source_is_cycle, current_weight

def target_from_gamma(G,gamma,psi):
    return tuple(G.smul(v,psi) for v in gamma)

def make_current(problem, entries):
    counts=[0]*len(problem.types)
    lookup={(t.plaquette,t.color,t.eps):i for i,t in enumerate(problem.types)}
    for key,n in entries:
        counts[lookup[key]] += n
    return Current(tuple(counts))

def test_boundary_squared_zero():
    B=OrientedBox((2,2,1))
    for p in range(len(B.plaquettes)):
        c=[0]*len(B.plaquettes); c[p]=1
        assert all(x==0 for x in B.boundary1(B.boundary2(c)))

def test_extraction_with_invisible_residual():
    B=OrientedBox((1,1)); G=FiniteAbelianGroup((2,2)); e1=(1,0); e2=(0,1); mix=(1,1)
    P=MulticolorProblem(B,G,(e1,e2,mix))
    gamma=B.rectangular_loop((0,0),(0,1),(1,1)); target=target_from_gamma(G,gamma,e1)
    # K=-(e1+e2)S plus K=+e2 S => K=-e1 S, hence source=e1 gamma.
    n=make_current(P,[((0,2,-1),1),((0,1,+1),1)])
    assert source(P,n)==target and source_is_cycle(P,target)
    N=subgroup_generated(G,(e2,))
    res=extract(P,n,target,N)
    assert source(P,res.extracted)!=target
    assert all(x in N for x in source(P,res.remainder))
    assert tuple(a+b for a,b in zip(res.extracted.counts,res.remainder.counts))==n.counts

def test_no_monotonicity_assumed_but_termination():
    B=OrientedBox((1,1)); G=FiniteAbelianGroup((5,)); X=((1,),(2,)); P=MulticolorProblem(B,G,X)
    gamma=B.rectangular_loop((0,0),(0,1),(1,1)); target=target_from_gamma(G,gamma,(2,))
    # Arrange two visible quanta whose net flux is -2 mod 5; deterministic type order may take a nonmatching step first.
    n=make_current(P,[((0,0,-1),2)])
    assert source(P,n)==target
    N=frozenset((G.zero,)); r=extract(P,n,target,N)
    assert len(r.steps)==2 and source(P,r.extracted)==target

def test_factorial_submultiplicativity_numeric():
    B=OrientedBox((1,1)); G=FiniteAbelianGroup((3,)); P=MulticolorProblem(B,G,((1,),))
    a=make_current(P,[((0,0,1),2)]); b=make_current(P,[((0,0,1),3)])
    c=Current(tuple(x+y for x,y in zip(a.counts,b.counts)))
    beta=(0.2,)
    assert current_weight(P,c,beta) <= current_weight(P,a,beta)*current_weight(P,b,beta)+1e-18
