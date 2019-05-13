# -*- coding: utf-8 -*-

import os
import sys

sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations")
sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/root")
sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/utils")

sys.path.append(r"C:\Users\daemo\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations")
sys.path.append(r"C:\Users\daemo\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/root")
sys.path.append(r"C:\Users\daemo\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/utils")

from wpg.generators import build_gauss_wavefront 
from utils.gauss_schell_utils import modes2D
from wpg.wavefront import Wavefront
from gvr_utils import *
import numpy as np

def modes2D(m,N):
    """
    return N pairs of Transverse Gauss-Hermite Mode Order pairs with order up to m. 

    """
    A=list()
    for i in range(1,m+1):
       A.append ( list(itertools.product([0,i],repeat=2)))
       
    # combine   
    A=list(itertools.chain.from_iterable(A))
    
    # remove duplicates
    temp = []
    for a,b in A:
        if (a,b) not in temp: #to check for the duplicate tuples
            temp.append((a,b))
    
    return temp[0:N]
    
def Modes2D(m, N):
    A = []
    A.append((0,0))
    while len(A)<N:
        A.append((np.random.randint(0,N), np.random.randint(0,N)))
    return A
        
def eigenvaluePartCoh(p_I, p_mu, n):
    """
    GVR
    
    return eigenvalue normalised to eigenvalue of fundamental
    
    definitions follow Starikov and Wolf, 1982:

    
    p_mu and p_I are the rms widths of the degree of coherence and of the intensity of the source.
    
    beta =  p_mu/p_I is a measure of the "degree of global coherence" of the source.

    When beta >> 1, the source is effectively spatially coherent in the global sense and is then
    found to be well represented by a single mode. 
    When beta << 1, the source is effectively spatially incoherent in the global
    sense, and the number of modes needed to describe its behavior is of the order of beta/3
     
    """
    
    a = 1/(2*p_I**2)
    b = 1/(2*p_mu**2)
    c = (a**2 + 2*a*b)**(1/2)
    l_0=1.0
    
    l_n = l_0*(b/(a+b+c))**n
    
    return l_n

def eigval(b,n):
    bb = b**2/2
    l = (1/(bb+1+b*np.sqrt(bb**2+1)))**n
    return l

### test
def modes(M, N):
    
    mod = []
    mod.append((0,0))
    
    for m in range(M):
        # mod.append((m,m)) this adds even (m,m) modes // ie., (1,1), (2,2) etc.
        for item in itertools.permutations(range(m), 2):
            mod.append(item)
    # modes = list(set(modes)) does not retain order of set
    
    ### long/rough soln'
    modes = []
    

    for item in mod:
            if item in modes:
                pass
            else:
                modes.append(item)
                
    modes = modes[0:N]
        
    return modes    

