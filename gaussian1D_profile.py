# -*- coding: utf-8 -*-
"""

Make a simple 1D gaussian profile. 

"""

import numpy as np
import matplotlib.pyplot as plt


def gaussian_1D_profile(x_min, x_max, x_step, center, sigma, amplitude):
    """Function to create a 1D Gaussian distribution. 
    
    Parameters
    ----------
    
    x_min, x_max, x_step: float, float, float
        Creates a sequence (1D ndarray) of points over which to compute the Gaussian
    center: float
        The center point of the gaussian profile
    sigma: float
        1/e-squared width of beam
    amplitude: float 
        Amplitude at peak value
        
    Returns
    -------
    
    x,y: ndarray
        the gaussian profile amplitude values
        
    """
    
    x = np.arange(x_min, x_max,x_step)  #create spatial array
    d = 2*float(sigma)
    y = amplitude*np.e**(-2*np.power((x-center)/d, 2))
    
    return x,y

    # todo: learn how to do proper unit testing...heres some manual checks
    # what if center > max(X)?  still works, just get the tail end
    # what if center, sigma negative?  Since is getting squared, doesn't matter
    # what if amplitude is neg or zero? Straight line at zero
    # what if d = 0? Straight line
    # what if the ndarray goes negative?  Is ok.
    # What if the array is empty or null? should catch an error.

def plot_1d_gaussian(x,y,hold=True):  
    """Plot the gaussian profile.
    
    Parameters
    ----------
    
    x: ndarray
        X axis values
    y: float
        Y axis values 
        
    """   
    plt.hold = hold
    plt.plot(x,y)
    plt.xlabel('X axis')
    plt.ylabel('Amplitude')
    plt.title('Gaussian 1D Profile')
    plt.show() 
    
    # todo: check if the hold true works or not
    

if __name__ == '__main__':
    
    x,y = gaussian_1D_profile(-50,50,.2, 0, 10, 1)
    plot_1d_gaussian(x,y,True)


