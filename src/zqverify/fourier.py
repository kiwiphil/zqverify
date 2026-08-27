"""Direct finite Fourier calculations for positive-character actions."""
from __future__ import annotations
import cmath, math
from itertools import product

def character(group,h,g):
    phase=sum((int(a)*int(x))/n for a,x,n in zip(h,g,group.moduli))
    return cmath.exp(2j*math.pi*phase)
def plaquette_weight(group,active_colors,betas,g):
    return math.exp(sum(float(b)*character(group,c,g).real for c,b in zip(active_colors,betas)))
def fourier_coefficient(group,active_colors,betas,eta):
    z=0j
    for g in group.elements:
        z += plaquette_weight(group,active_colors,betas,g)*character(group,eta,g).conjugate()
    return z/group.order
def normalized_fourier(group,active_colors,betas,eta):
    return (fourier_coefficient(group,active_colors,betas,eta)/fourier_coefficient(group,active_colors,betas,group.zero)).real
def exact_2d_rectangle_wilson(group,active_colors,betas,psi,area):
    return normalized_fourier(group,active_colors,betas,psi)**int(area)
