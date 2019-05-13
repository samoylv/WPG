import numpy as np

from skimage import draw

from matplotlib import pyplot as plt

class Photomask():
    def __init__(self, nx, ny):
        
        self.nx, self.ny = nx, ny
        self.mask = np.zeros((self.nx, self.ny), dtype = np.uint8)    
            
        
    def slit(self, dx, dy, _x = 0):
        """
        slit dimensions dx,dy,_x are in m
        """

        maxi = int(len(self.mask))
        
        dx = int(dx)
        dy = int(dy)
        xshift = int(dx/2)  
        yshift = int(dy/2) 

        _x = int(_x)
        
        xc = int(maxi/2)
        yc = int(len(self.mask)-(maxi/2))

        start = (yc - yshift, _x + xc - xshift)
        end = (yc + yshift, _x + xc + xshift)
        
        r,c = draw.rectangle(start = start, end = end, shape = self.mask.shape)
        self.mask[r,c] = 1

    def nslits(self, dx, dy, _a = 100, n = 2):
    
        for N in range(n):
            
            self.slit(dx,dy,_x = (-30*N*2)-int(_a))
            self.slit(dx,dy,_x = (+30*N*2)+int(_a))
            
    def save_array(self, filename):
        np.save(filename, self.mask)
        
def pixels2rspace(pixelsize, actualsize):
    npixels = actualsize/pixelsize
    return npixels



if __name__ == "__main__":
    
    Mask = Photomask(500,500)
    Mask.nslits(pixels2rspace(1.011e-09, 30e-09), 
                pixels2rspace(1.011e-09, 50e-09), 
                pixels2rspace(1.011e-09, 100e-09),
                3)
    plt.imshow(Mask.mask)
    Mask.save_array("in/hajaa")
    