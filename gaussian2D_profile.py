# -*- coding: utf-8 -*-
"""
Make a simple 2D Gaussian profile.  One example is in the representation of an
ideal laser beam.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import gaussian1D_profile as gp




def gaussian_2D_profile(x_min, x_max, x_step, y_min, y_max, y_step,
                        xc, yc, sigma_x, sigma_y, amplitude):
    """Function to produce a 2D Gaussian profile.
    
    Parameters
    ----------
    
    x_min, x_max, x_step: float, float, float
        Creates a sequence (1D ndarray) of x points over which to compute the Gaussian
    y_min, y_max, y_step: float, float, float
        Creates a sequence (1D ndarray) of y points over which to compute the Gaussian        
    xc,yc: float, float
        The center points in x and y of the gaussian profile
    sigma_x, sigma_y: float, float
        1/e-squared width of beam in the x and y axis respectively
    amplitude: float 
        Amplitude at peak value
        
    Returns
    -------
    
    Z: ndarray
        The 2D gaussian profile amplitude values
        
    """

    idx = np.arange(x_min, x_max, x_step)
    idy = np.arange(y_min, y_max, y_step)

    d0x = 2*sigma_x  #beam diameter (1/e-squared) is 2 * standard deviation
    d0y = 2*sigma_y

    x = amplitude * np.e**(-2*np.power((idx-xc)/d0x,2))
    y = amplitude * np.e**(-2*np.power((idy-yc)/d0y,2))


    plt.plot(idx, x, label='x')
    plt.plot(idy, y, label='y')
    plt.title('Gaussian 1D Profiles in x and y')
    plt.legend()
    plt.show()

    #now 2d...can just multiply the two 1D curves together
    X,Y = np.meshgrid(idx,idy)
    
    #Z = np.e**(-2*(((X-xc)/d0x)**2 * ((Y-yc)/d0y)**2))  #interesting star type pattern

    # The following methods are equivalent e^x * e^y = e^(x+y)
    #Z = np.e**(-2*np.power((X-xc)/d0x,2))*np.e**(-2*np.power((Y-yc)/d0y,2)) #1.92ms per loop
    #Z = np.e**(-2*((X-xc)/d0x)**2) * np.e**(-2*((Y-yc)/d0y)**2) # 1.03ms per loop
    Z = np.e**(-2*(((X-xc)/d0x)**2 + ((Y-yc)/d0y)**2))          # 586us per loop (4x faster)
    
    return Z


def plot_2d_gaussian(Z):
    """Plot the gaussian profile as a 2D image and contour lines.
    
    Parameters
    ----------
    
    Z: ndarray
        2D Gaussian values
        
    """   
    
    plt.figure()
    #plt.axis('off')
    plt.imshow(Z,extent=(-50,50,-50,50))
    plt.contour(Z)
    plt.title('2D Gaussian Profile')
    plt.show()


if __name__ == '__main__':
    Z = gaussian_2D_profile(-50,50,1,-50,50,1, 0, 0, 10, 5, 1)
    plot_2d_gaussian(Z)
    
    
    
    