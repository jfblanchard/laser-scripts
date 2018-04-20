#!/usr/bin/env python3
"""

Functions for computing commonly used laser beam parameters.

Todo: Maybe merge the bpp functions for radius and diameter.

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


def bpp_raduis_halfdiv(waist_rad, half_divergence):
    """The Beam Parameter Product (bpp) is defined as the product of the beam 
    waist size and the far field divergence. This function uses the beam radius
    and half-divergence, but bpp can also be computed from beam diameter and 
    full angle divergence (see below).  The units are generally mm-mrad. 
    
    Parameters
    ----------
    waist_rad : float
        The waist radius of the beam in meters.
    half_divergence : float
        The half divergence of the beam in radians.
        
    Returns
    -------
    bpp : float
        The Beam parameter product in mm-mrad
    """
    
    bpp = waist_rad*mm * half_divergence*mm
    return bpp

    
def bpp_diam_fulldiv(waist_diam, full_divergence):
    """The Beam Parameter Product (bpp) is defined as the product of the beam 
    waist size and the far field divergence.  This function uses the beam 
    diameter and full angle divergence.  The units are generally mm-mrad. 
    
    Parameters
    ----------
    waist_diam : float
        The waist diameter of the beam in meters.
    full_divergence : float
        The half divergence of the beam in radians.
        
    Returns
    -------
    bpp : float
        The Beam parameter product in mm-mrad
    """
    
    bpp = 4* waist_diam*mm * full_divergence*mm
    return bpp


def bpp_wl_m2(wavelength, msquared=1.0):
    """This function computes the beam parameter product of a gaussian beam, 
    which is defined as the center wavelength divided by Pi.  The bpp of a beam 
    with an M-squared > 1 is just the bpp value multiplied by the m-squared
    value.
    
    
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


def beam_waist(wavelength,divergence, msquared=1.0):
    """Compute the beam waist from wavelength, divergence, and 
    m-squared value. 
    
    Parameters
    ----------
    wavelength : float
        The wavelength in meters
    divergence : float
        The half divergenence in rad       
    msquared : float, optional
        The msquared value of the beam
        
    Returns
    -------
    waist : float
        The beam waist radius in meters
    """
    
    waist = msquared*wavelength/(np.pi * divergence)
    return waist

    def rayleigh_range(wavelength,waist):
    """Compute the rayleigh range of a gaussian beam in free space.  This is 
    the distance where the cross section doubles, and radius increases sqrt(2).
    
    Parameters
    ----------
    wavelength : float
        The wavelength in meters
    waist : float
        The waist radius in meters    
        
    Returns
    -------
    zR : float
        The rayleigh range in meters
    """
    
    zR = (np.pi*waist**2)/wavelength
    return zR

if __name__ == '__main__':
    #do some test cases
    
    

