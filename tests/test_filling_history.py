import pytest
from zqverify import *
from zqverify.currents import MulticolorProblem
from zqverify.verification import rectangle_target, minimum_projected_filling_size, minimum_projected_filling_cost, history_level_activity

@pytest.mark.parametrize('q,r,R,T',[(3,1,1,1),(4,2,1,1),(5,2,1,1),(5,2,1,2),(6,3,1,2)])
def test_projected_filling_cyclic(q,r,R,T):
    G=FiniteAbelianGroup((q,)); X=((1,),); B=OrientedBox((R,T)); P=MulticolorProblem(B,G,X); N=frozenset((G.zero,))
    target=rectangle_target(P,(0,0),(0,1),(R,T),(r,))
    got=minimum_projected_filling_size(P,target,N)
    assert got==R*T*min(r,q-r)


def test_projected_filling_product_quotient():
    G=FiniteAbelianGroup((2,2)); X=((1,0),(0,1),(1,1)); B=OrientedBox((1,1)); P=MulticolorProblem(B,G,X)
    N=subgroup_generated(G,((0,1),)); target=rectangle_target(P,(0,0),(0,1),(1,1),(1,0))
    assert minimum_projected_filling_size(P,target,N)==1


def test_weighted_projected_filling():
    G=FiniteAbelianGroup((2,2)); X=((1,0),(0,1)); B=OrientedBox((1,1)); P=MulticolorProblem(B,G,X); N=frozenset((G.zero,))
    target=rectangle_target(P,(0,0),(0,1),(1,1),(1,1)); costs=(0.3,0.7)
    assert abs(minimum_projected_filling_cost(P,target,N,costs)-1.0)<1e-13

@pytest.mark.parametrize('depth',[0,1,2,3,4])
def test_history_level_branch_bound(depth):
    G=FiniteAbelianGroup((3,)); X=((1,),); B=OrientedBox((1,1)); P=MulticolorProblem(B,G,X); N=frozenset((G.zero,))
    target=rectangle_target(P,(0,0),(0,1),(1,1),(1,)); beta=(0.08,); rho=2*0.08
    H=history_level_activity(P,target,N,beta,depth)
    assert H <= rho**depth + 1e-14
