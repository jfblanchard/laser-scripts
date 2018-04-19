#!/usr/bin/env python3
"""

Functions for computing commonly used laser beam parameters.

"""

import numpy as np
import scipy.constants as const


#define units...maybe use one of those units modules instead
m = 1           
cm = 1e2      
mm = 1e3       
um = 1e6
nm = 1e9
pm = 1e12


def bpp(waist_rad,half_divergence):
    """The Beam Parameter Product (bpp) is defined as the product of the beam 
    waist radius and the half angle divergence.  The BPP for an ideal 
    (diffraction-limited) gaussian beam is the minimum value possible, which is 
    defined as the center wavelength divided by Pi.  This function computes the
    bpp from wavelength and m-squared value. This function uses the most beam 
    radius and half-divergence, which is the most common, but bpp can also be 
    computed from beam diameter and full angle divergence (see below).  The
    units are generally mm-mrad. 
    
    Parameters
    ----------
    waist_rad : float
        The waist radius of the beam in meters.
    half_divergence : float
        The half divergence of the beam in radians.
        
    Returns
    -------
    BPP : float
        The Beam parameter product in mm-mrad
    """
    
    bpp = waist_rad*mm * half_divergence*mm
    return bpp
    


def bpp_wl_m2(wavelength, msquared=1.0):
    """This function computes the beam parameter product of a gaussian beam 
    from its center wavelength and m-squared value. 
    
    
    Parameters
    ----------
    wavelength : float
        The center wavelength of the laser in meters
    msquared : float, optional, default=1.0
        The msquared value of the beam
        
    Returns
    -------
    bpp : float
        The Beam parameter product in mm-mrad
    """
    
    bpp = (wavelength*mm * msquared*mm)/np.pi
    return bpp


def beam_half_divergence(wavelength,waist, msquared=1.0):
    """Compute the half angle divergence from wavelength, waist radius, and 
    m-squared value. 
    
    Parameters
    ----------
    wavelength : float
        The wavelength in meters
    waist : float
        The waist radius in meters     
    msquared : float, optional
        The msquared value of the beam
        
    Returns
    -------
    divergence : float
        The half-angle divergence in rad
    """
    theta = msquared*wavelength/(np.pi * waist)
    return theta


if __name__ == '__main__':
    #do some test cases
    
    

