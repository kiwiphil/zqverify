import pytest
from zqverify import *
from zqverify.fourier import normalized_fourier

cyclic_cases=[(q,r) for q in range(2,9) for r in range(q)]
@pytest.mark.parametrize('q,r',cyclic_cases)
def test_cyclic_word_metric_grid(q,r):
    G=FiniteAbelianGroup((q,)); X=((1,),); N=frozenset((G.zero,))
    assert signed_word_length(G,X,(r,),N)==min(r,q-r)

@pytest.mark.parametrize('q',[2,3,4,5,6,7])
def test_cyclic_fourier_nonnegative(q):
    G=FiniteAbelianGroup((q,)); X=((1,),)
    vals=[normalized_fourier(G,X,(0.09,),eta) for eta in G.elements]
    assert min(vals)>=-1e-13
    assert abs(vals[0]-1.0)<1e-13

@pytest.mark.parametrize('moduli,expected', [((2,),2),((3,),2),((4,),3),((2,2),5),((2,4),8),((3,3),6)])
def test_subgroup_counts_known_small_groups(moduli,expected):
    G=FiniteAbelianGroup(moduli)
    assert len(enumerate_subgroups(G))==expected
