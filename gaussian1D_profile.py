# -*- coding: utf-8 -*-
"""

Make a simple 1D gaussian profile. 

"""

import numpy as np
import matplotlib.pyplot as plt


def gaussian_1D_profile(X, center, sigma, amplitude):
    """Function to create a 1D Gaussian distribution. 
    
    Parameters
    ----------
    
    X: ndarray
        The points over which to compute the Gaussian
    center_x: float
        The center point of the gaussian profile
    sigma: float
        1/e-squared width of beam
    amplitude: float 
        Amplitude at peak value
        
    Returns
    -------
    
    y: ndarray
        the gaussian profile amplitude values
        
    """

    d = 2*float(sigma)
    y = amplitude*np.e**(-2*np.power((X-center)/d, 2))
    return y

    # todo: learn how to do proper unit testing...heres some manual checks
    # what if center > max(X)?  still works, just get the tail end
    # what if center, sigma negative?  Since is getting squared, doesn't matter
    # what if amplitude is neg or zero? Straight line at zero
    # what if d = 0? Straight line
    # what if the ndarray goes negative?  Is ok.
    # What if the array is empty or null? should catch an error.

def plot_1d_gaussian(x,y):
    """Plot the gaussian profile.
    
    Parameters
    ----------
    
    x: ndarray
        X axis values
    y: float
        Y axis values 
        
    """   
    
    plt.plot(x,y)
    plt.xlabel('X axis')
    plt.ylabel('Amplitude')
    plt.title('Gaussian 1D Profile')
    plt.show() 
    

if __name__ == '__main__':

    x = np.arange(-50,50,.2)         #create spatial array
    y = gaussian_1D_profile(x, 0, 10, 1)
    plot_1d_gaussian(x,y)


