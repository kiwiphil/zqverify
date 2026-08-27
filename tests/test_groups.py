import math
from zqverify import *

def test_subgroups_z2x4():
    G=FiniteAbelianGroup((2,4))
    Hs=enumerate_subgroups(G)
    assert frozenset((G.zero,)) in Hs
    assert frozenset(G.elements) in Hs
    assert len(Hs)==8

def test_word_lengths_product():
    G=FiniteAbelianGroup((2,4)); X=((1,0),(0,1)); N=frozenset((G.zero,))
    assert signed_word_length(G,X,(1,1),N)==2
    assert signed_word_length(G,X,(0,3),N)==1

def test_quotient_z6_even_kernel():
    G=FiniteAbelianGroup((6,)); X=((1,),(2,),(4,)); N=subgroup_generated(G,((2,),))
    assert len(N)==3
    assert signed_word_length(G,X,(1,),N)==1
    assert quotient_is_cyclic(G,N)

def test_weighted_zero_cost_pseudometric():
    G=FiniteAbelianGroup((2,2)); X=((1,0),(0,1)); N=frozenset((G.zero,))
    assert weighted_word_length(G,X,(0.0,1.0),(1,0),N)==0.0

def test_noncyclic_detector_example():
    G=FiniteAbelianGroup((2,2,2)); X=((1,0,0),(0,1,0),(0,0,1)); psi=(1,1,1)
    x=0.03
    a=best_detector(G,X,(x,x,x),psi,False)
    c=best_detector(G,X,(x,x,x),psi,True)
    assert a and c
    assert not a.quotient_cyclic
    assert a.word_length==3 and abs(a.rho-0.09)<1e-15
    assert c.word_length==1 and abs(c.rho-0.03)<1e-15
    assert a.coefficient > c.coefficient

def test_mixed_order_noncyclic_example():
    G=FiniteAbelianGroup((2,4)); X=((1,0),(0,1)); psi=(1,1); x=0.03
    a=best_detector(G,X,(x,x),psi,False); c=best_detector(G,X,(x,x),psi,True)
    assert a and c and not a.quotient_cyclic
    assert a.word_length==2 and c.word_length==1
    assert a.coefficient > c.coefficient
