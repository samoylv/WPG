# -*- coding: utf-8 -*-

import numpy as np

from matplotlib import pyplot as plt
from skimage import draw

class Mask():
    """ A CLASS FOR MAKING 2D DIFFRACTION GRATING IMAGES (SQUARE ONLY)"""
    
    def __init__(self, nx = 256, ny = 256, xdim = 1, ydim = 1, VERBOSE = True):
        
        self.nx = nx
        self.ny = ny
        
        self.xdim = xdim
        self.ydim = ydim
        
        self.res = self.xdim/self.nx
        
        self.array = np.zeros((self.nx, self.ny), dtype=np.int8)
        
        self.center = int(self.nx/2)  # center of array
        
        if VERBOSE == True:
            print("Mask Dimensions: {} x {}".format(self.nx, self.ny))
            print("Mask Size: {} m x {} m".format(self.xdim, self.ydim))
            print("Mask Resolution: {} m".format(self.res))
        
        
    def DoubleSlit(self, width, height, sep):
        
        self.width = width/self.res
        self.height = height/(self.res*2)
        self.sep = sep/(self.res*2)
        
        r1 = (self.center - self.sep, self.center - self.height)
        r2 = (self.center - self.sep + self.width, self.center + self.height)
        
        r,c = draw.rectangle(start = r1, end = r2, shape = self.array.shape)
        
        self.array[r.astype(dtype=np.int16), c.astype(dtype=np.int16)] = 1
        plt.imshow(self.array)
    def DS(self, width, height, period, sep, n):
        
        """WRITE GRATINGS INTO ARRAY
        --- n is the number of gratings either side of the central stop
        --- sep is the width of the central stop"""
    
        self.width = width
        self.height = height
        
        self.period = np.around(period/(self.resolution), decimals = 0)
        self.sep = np.around(sep/(self.resolution*2), decimals = 0)

        self.n = n
        
        rx = np.around(self.width/(self.resolution), decimals = 0)
        ry = np.around(self.height/(self.resolution), decimals = 0)
        
        for i in range(self.n):

            # NEGATIVE GRATINGS
            start = (self.C + self.sep + i*(self.period+rx), self.C-ry)
            end = (self.C + self.sep + i*(self.period+rx)+rx, self.C+ry)
            
            r,c = draw.rectangle(start = start, end = end, shape=self.array.shape)
            
            self.array[r.astype(dtype=np.int16), c.astype(dtype=np.int16)] = 1
        
        self.array = np.rot90(self.array) + np.rot90(np.rot90(np.rot90(self.array))) 
        plt.set_cmap('bone')
        plt.imshow(self.array)
        
        return self.array

    def QS(self, width, height, period, sep, n):
        
        """WRITE GRATINGS INTO ARRAY
        --- n is the number of gratings either side of the central stop
        --- sep is the width of the central stop"""
    
        self.width = width
        self.height = height
        
        self.period = np.around(period/(self.resolution), decimals = 0)
        self.sep = np.around(sep/(self.resolution*2), decimals = 0)

        self.n = n
        
        rx = np.around(self.width/(self.resolution), decimals = 0)
        ry = np.around(self.height/(self.resolution), decimals = 0)
        
        for i in range(self.n):

            # NEGATIVE GRATINGS
            start = (self.C-ry, self.C - self.sep - (i+1)*(self.period+rx)-rx)
            end = (self.C+ry, self.C - self.sep - (i+1)*(self.period+rx)-rx)
            
            r,c = draw.rectangle(start = start, end = end, shape=self.array.shape)
            
            self.array[r.astype(dtype=np.int16), c.astype(dtype=np.int16)] = 1
        
        self.array += np.rot90(np.rot90(self.array))
        self.array += np.rot90(self.array)
        
        plt.set_cmap('bone')
        plt.imshow(self.array)
        
        return self.array
    
if __name__ == "__main__":
    pass
