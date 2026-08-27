from zqverify import *
from zqverify.search import generating_active_sets

def test_active_sets_generate():
    for mod in ((2,2),(2,4),(3,3)):
        G=FiniteAbelianGroup(mod)
        for X in generating_active_sets(G,3):
            assert len(subgroup_generated(G,X))==G.order

def test_z6_rescue_example():
    G=FiniteAbelianGroup((6,)); X=((1,),(2,),(4,)); psi=(1,)
    b=best_detector(G,X,(0.03,0.55,0.55),psi)
    assert sum((0.03,0.55,0.55))>=1
    assert b is not None and abs(b.rho-0.03)<1e-15
