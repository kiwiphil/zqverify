from zqverify import *
from zqverify.currents import MulticolorProblem
from zqverify.flux import exact_flux_wilson
from zqverify.fourier import exact_2d_rectangle_wilson

def group_target(G,gamma,psi): return tuple(G.smul(v,psi) for v in gamma)

def test_flux_matches_2d_oneplaquette():
    G=FiniteAbelianGroup((3,)); X=((1,),); B=OrientedBox((1,1)); P=MulticolorProblem(B,G,X)
    gamma=B.rectangular_loop((0,0),(0,1),(1,1)); target=group_target(G,gamma,(1,))
    beta=(0.1,)
    assert abs(exact_flux_wilson(P,beta,target)-exact_2d_rectangle_wilson(G,X,beta,(1,),1))<1e-13

def test_one_cube_z2_detector_bound():
    G=FiniteAbelianGroup((2,)); X=((1,),); B=OrientedBox((1,1,1)); P=MulticolorProblem(B,G,X)
    gamma=B.rectangular_loop((0,0,0),(0,1),(1,1)); target=group_target(G,gamma,(1,))
    for beta in (0.02,0.05,0.10,0.20):
        W=exact_flux_wilson(P,(beta,),target)
        rho=4*beta
        assert W <= rho/(1-rho)+1e-12
