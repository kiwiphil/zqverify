"""Finite detector search catalogue used by the release report."""
from __future__ import annotations
from itertools import combinations, product
from .groups import *
from .detector import best_detector

def generating_active_sets(group,max_size=3):
    reps=canonical_signed_representatives(group)
    out=[]
    for k in range(1,min(max_size,len(reps))+1):
        for X in combinations(reps,k):
            if len(subgroup_generated(group,X))==group.order:
                out.append(X)
    return tuple(out)

def run_catalog(activity_levels=(0.03,0.12,0.55)):
    specs=[(2,),(3,),(4,),(5,),(6,),(2,2),(2,2,2),(2,4),(3,3)]
    total=adv=rescued=0
    by_group={}
    for moduli in specs:
        G=FiniteAbelianGroup(moduli); sub=adv_g=res_g=0
        for X in generating_active_sets(G,3):
            for psi in G.elements:
                if psi==G.zero: continue
                for acts in product(activity_levels, repeat=len(X)):
                    total+=1; sub+=1
                    arbitrary=best_detector(G,X,acts,psi,False)
                    cyclic=best_detector(G,X,acts,psi,True)
                    if sum(acts)>=1 and arbitrary is not None:
                        rescued+=1; res_g+=1
                    if arbitrary and cyclic and arbitrary.coefficient>cyclic.coefficient+1e-12 and not arbitrary.quotient_cyclic:
                        adv+=1; adv_g+=1
        by_group[str(moduli)]={'problems':sub,'noncyclic_advantage':adv_g,'rescued':res_g}
    return {'problems':total,'noncyclic_advantage':adv,'crude_rho_ge_1_but_detector':rescued,'by_group':by_group}
