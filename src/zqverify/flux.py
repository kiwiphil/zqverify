"""Tiny-volume exact flux enumeration."""
from itertools import product
from .fourier import fourier_coefficient

def exact_flux_wilson(problem,betas,target):
    g=problem.group; box=problem.box
    hats={eta:fourier_coefficient(g,problem.active_colors,betas,eta).real for eta in g.elements}
    Z0=0.0; Zj=0.0
    for labels in product(g.elements, repeat=len(box.plaquettes)):
        # source = -partial labels
        src=[]
        for e,row in enumerate(box.boundary2_matrix):
            s=g.zero
            for p,inc in enumerate(row):
                if inc: s=g.add(s,g.smul(-inc,labels[p]))
            src.append(s)
        w=1.0
        for x in labels: w*=hats[x]
        if all(x==g.zero for x in src): Z0+=w
        if tuple(src)==tuple(target): Zj+=w
    return Zj/Z0
