from itertools import product
from zqverify import *
from zqverify.fourier import normalized_fourier, exact_2d_rectangle_wilson

def test_fourier_coefficients_positive_small_groups():
    for mod in ((2,),(3,),(4,),(2,2)):
        G=FiniteAbelianGroup(mod)
        X=tuple(tuple(1 if i==j else 0 for i in range(len(mod))) for j in range(len(mod)))
        b=tuple(0.08 for _ in X)
        for eta in G.elements:
            assert normalized_fourier(G,X,b,eta) >= -1e-13

def test_detector_bounds_2d_product_groups():
    levels=(0.03,0.12,0.55); Bd=2.0
    for mod in ((2,2),(2,4),(3,3)):
        G=FiniteAbelianGroup(mod)
        X=tuple(tuple(1 if i==j else 0 for i in range(len(mod))) for j in range(len(mod)))
        for psi in G.elements:
            if psi==G.zero: continue
            for eff in product(levels, repeat=len(X)):
                best=best_detector(G,X,eff,psi)
                if best is None: continue
                betas=tuple(x/Bd for x in eff)
                W=exact_2d_rectangle_wilson(G,X,betas,psi,1)
                bound=best.rho**best.word_length/(1-best.rho)
                assert -1e-13 <= W <= bound+1e-12

def test_factorized_invisible_activity_exact():
    G=FiniteAbelianGroup((2,2)); X=((1,0),(0,1)); psi=(1,0)
    w1=exact_2d_rectangle_wilson(G,X,(0.08,0.1),psi,2)
    w2=exact_2d_rectangle_wilson(G,X,(0.08,2.0),psi,2)
    assert abs(w1-w2)<1e-13
