#!/usr/bin/env python3
from __future__ import annotations
import json
from itertools import product
from pathlib import Path
from zqverify import *
from zqverify.currents import MulticolorProblem
from zqverify.fourier import exact_2d_rectangle_wilson
from zqverify.flux import exact_flux_wilson

def basis(mod): return tuple(tuple(1 if i==j else 0 for i in range(len(mod))) for j in range(len(mod)))
def target(G,gamma,psi): return tuple(G.smul(v,psi) for v in gamma)
rows=[]
# 162 genuinely quotient-sensitive one-plaquette d=2 cases.
for mod in ((2,2),(2,4),(3,3)):
    G=FiniteAbelianGroup(mod); X=basis(mod); Bd=2.0
    for psi in G.elements:
        if psi==G.zero: continue
        for eff in product((0.03,0.12,0.55), repeat=len(X)):
            b=best_detector(G,X,eff,psi)
            if b is None: continue
            W=exact_2d_rectangle_wilson(G,X,tuple(x/Bd for x in eff),psi,1)
            bound=b.rho**b.word_length/(1-b.rho)
            rows.append(dict(kind='product-d2',moduli=mod,psi=psi,effective=eff,W=W,bound=bound,ratio=W/bound,kernel_size=len(b.kernel),cyclic=b.quotient_cyclic))
# 48 cyclic d=2 area 1 and 2 cases.
for q in (2,3,4):
    G=FiniteAbelianGroup((q,)); X=((1,),); Bd=2.0
    for r in range(1,q):
        for beta in (0.05,0.10,0.20,0.40):
            eff=(Bd*beta,); b=best_detector(G,X,eff,(r,))
            if b is None: continue
            for area in (1,2):
                W=exact_2d_rectangle_wilson(G,X,(beta,),(r,),area)
                bound=b.rho**(b.word_length*area)/(1-b.rho)
                rows.append(dict(kind='cyclic-d2',moduli=(q,),psi=(r,),beta=beta,area=area,W=W,bound=bound,ratio=W/bound,kernel_size=len(b.kernel),cyclic=b.quotient_cyclic))
# 4 d=3 Z2 one-cube cases.
G=FiniteAbelianGroup((2,)); X=((1,),); B=OrientedBox((1,1,1)); P=MulticolorProblem(B,G,X)
gamma=B.rectangular_loop((0,0,0),(0,1),(1,1)); tgt=target(G,gamma,(1,))
for beta in (0.02,0.05,0.10,0.20):
    eff=(4*beta,); b=best_detector(G,X,eff,(1,)); W=exact_flux_wilson(P,(beta,),tgt)
    bound=b.rho**b.word_length/(1-b.rho)
    rows.append(dict(kind='z2-d3-cube',moduli=(2,),psi=(1,),beta=beta,W=W,bound=bound,ratio=W/bound,kernel_size=len(b.kernel),cyclic=b.quotient_cyclic))
assert len(rows)==214, len(rows)
assert all(r['W'] <= r['bound']+1e-12 for r in rows)
summary={
 'cases':len(rows), 'passing':sum(r['W']<=r['bound']+1e-12 for r in rows),
 'max_ratio':max(r['ratio'] for r in rows),
 'product_cases':sum(r['kind']=='product-d2' for r in rows),
 'nontrivial_kernel_product_cases':sum(r['kind']=='product-d2' and r['kernel_size']>1 for r in rows),
 'noncyclic_quotient_product_cases':sum(r['kind']=='product-d2' and not r['cyclic'] for r in rows),
}
path=Path(__file__).resolve().parents[1]/'reports'/'finite_volume_bound_audit.json'
path.write_text(json.dumps({'summary':summary,'rows':rows},indent=2,default=list))
print(json.dumps(summary,indent=2))
