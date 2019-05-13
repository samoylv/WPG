import os
import sys


def twg_path():
    sys.path.insert(0,os.path.join('..','..'))
    
    try:
        sys.path.append(r"C:\Users\daemo\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/utils")
        sys.path.append(r"C:\Users\daemo\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/root")    
    except KeyError:
        print("Directory not Found")
    
    try:
        sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/utils")
        sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/root") 
        sys.path.append(r"C:\Users\twguest\OneDrive\Documents\Github\WPG_Wavefront_Simulations")    
    except KeyError:
        print("Directory not Found")
