#!/usr/bin/env python3
"""

Functions for calculating basic laser beam parameters.


"""

import numpy as np
import scipy.constants as const


#define units
m = 1           
cm = 1e2      
mm = 1e3       
um = 1e6
nm = 1e9
pm = 1e12



def initialize_beam(waist,wavelength,msquared=1.0):
    """Initialize a beam using waist, wavelength, and (optionally) m-squared.
    
    """  
    
    
def gaussian_profile_1D(P,waist,r):
    """ Return the optical power as measured at a distance r from the axis of
    the beam    
    
    Parameters
    ----------
    P : float
        The total power of the beam in watts
    waist : float
        The beam waist parameter in mmKe
    r: float
        The distance from the axis of the beam       
        
    Returns
    -------
    I_r: the optical power per unit area measured at dist r from the axis
    
    Usage
    -----
    r = np.linspace(-1,1,1000)
    g = np.zeros(len(r))
    for i in range(len(r)):
        g[i] = GaussianBeam(1,.1,r[i])    
    
    """
    I_r = (2*P/np.pi*waist**2)*np.e**(-2*r**2/waist**2)
    
    return I_r
    

def bpp_radius_halfdiv(waist_rad, half_divergence):
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


def beam_size(wavelength,waist,z):
    """Compute the beam size (radius) at location z from wavelength and 
    waist radius.  
    
    Parameters
    ----------
    wavelength : float
        The wavelength in microns
    waist : float
        The waist radius in microns      
    z: float
        The z location in microns at which to compute the beam radius           
        
    Returns
    -------
    w_z : float
        The beam radius at location z
    """  
    zR = get_rayleigh_range(wavelength,waist)
    w_z = waist * np.sqrt(1 + (z/zR)**2)
    return w_z
    
    #make this handle arrays too.
    
    
def q_parameter(wavelength, waist,z):
    """Compute the q-parameter from wavelength, waist radius, and z location.
    
    Parameters
    ----------
    wavelength : float
        The wavelength in microns
    waist : float
        The waist radius in microns      
    z: float
        The z location in microns at which to compute the q-parameter          
        
    Returns
    -------
    q_z : float
        The complex q-parameter at location z. Using q(z) = z + zR*i
    """  
    re = z
    zR = (np.pi*waist**2)/wavelength
    im = zR
    q_z = complex(re,im)
    return q_z
    

def peak_power(avg_power, pulse_width, rep_freq):
    """Compute peak power given average power, pulse width, and rep rate
    
    Parameters
    ----------
    avg_power : float
        The total power of the beam in watts
    pulse_with : float
        The FWHM pulse width of the laser in seconds
    rep_freq: float
        The pulse repetition frequency in Hz  
        
    Returns
    -------
    peak_pwr : the peak power"""
    
    dc = pulse_width * rep_freq    
    peak_pwr = avg_power / dc
    return peak_pwr


def intensity(P,radius):
    """For flat top, this is just the average power divided by the beam area.
    Units are W/cm**2.  Assumes Power is W, and radius is in m
    
    """
    w = radius*100  #convert from m to cm
    I = P/(np.pi*w**2)
    return I


#change to special case of intensity?
def peak_intensity(P,radius):
    """For Gaussian beams, the peak intensity is 2x that of a flat top.  I.e.:
    P/((pi*radius**2)/2).
    Power is in W, radius in m.  Returns W/cm**2
    """
    w = radius*100   #convert from m to cm
    return P/((np.pi*w**2)/2)
    
def fluence(J,radius):
    """The fluence of a laser pulse is the energy (J) divided by the beam area.
    Units are J/cm**2.  Assumes energy is J, and radius is in m.
    
    """
    w = radius*100  #convert from m to cm
    F = J/(np.pi*w**2)
    return F
    

def peak_fluence(J,radius):
    """For Gaussian beams, the peak fluence is J/((pi*radius**2)/2).
    Energy is in J, radius in m.  Returns J/cm**2
    """
    w = radius*100   #convert from m to cm
    return J/((np.pi*w**2)/2)
   
    
def beamsize_power_through_app(P_rat,r_app):
	"""Return the spot size diameter given the power ratio through an aperture 
    and aperture size.  Assumes a gaussian beam with m-squared =1.
	
	Parameters
	----------	
	P_rat : float
		The power ratio through the aperture
	r_app : float
		The aperture size
		
	Returns
	-------
	w0 : float
		Spot size diameter
	"""
	
	w0 = np.sqrt(-2*r_app**2/np.log(1-P_rat))
	return w0

def mode_frequency_separation(L):
    """ Compute the frequency separation of modes due to a laser cavity of length L 
    
	Parameters
	----------	
	L : float
		The optical cavity length in meters

	Returns
	-------
	sep : float
		Frequency separation in Hz
	"""    
 
    return const.c/(2*L)
     
    
def freq_to_wl(freq):
    
    """Calculate the wavelength given the freqency of the EM radiation. Assumes
     vacuum (n=1).  Uses the equation c = wl*freq.
	
	Parameters
	----------	
	freq : float
		The frequency of the electromagnetic wave in Hz

	Returns
	-------
	wl : float
		The wavelength (m) (in vacuum) corresponding the the input frequency
     """
         
    wl = const.c/freq
    return wl
 
 
def wl_to_freq(wl):
    """Calculate the frequency given the wavelength of the EM radiation. Assumes
     vacuum (n=1).  Uses the equation c = wl*freq.
	
	Parameters
	----------	
	freq : float
		The frequency of the electromagnetic wave in Hz

	Returns
	-------
	freq : float
		The wavelength (in vacuum) corresponding the the input frequency
     """
     
    freq = const.c/wl
    return freq
   
    
def wl_to_wavenumber(wl, angular=False):
    """Given wavelength in meters, convert to wavenumber in 1/cm.  The wave 
    number represents the number of wavelengths in one cm.  If angular is true,
    will calculate angular wavenumber.  """
    
    if angular:
        wnum = (2*np.pi)/(wl*100)
    else:  
        wnum = 1/(wl*100)
    
    return wnum


 
def photon_energy(wavelength):
    """Return the photon energy in joules given wavelength in microns"""
    
    wavelength = wavelength/1e6
    energy = (const.h * const.c)/ wavelength
    return energy
 
 
def emitted_frequency(E1,E2):
    """Calculate the frequency of an emitted photon due to a transition between 
    energy levels E1 and E2
    """
    
    freq = (E1 - E2)/const.hbar
    return freq
    
    #related: transition energy ~ 1.24/wl  so, 1um light is about 1.24eV
  

     
def boltzman_statistics(wl,K):
    """Calculate the ratio of ions in two energy states E1 and E2 in thermal 
    equilibrium given the energy difference defined by a photon of wavelength
    wl.  """
    
    energy = (const.h * const.c)/wl   #energy of a photon (joules) of 1.064um light
    ratio = np.exp(-1*energy/(const.k*K))
    print 'The ratio of ions in energy level 2 compared to 1 is: ' + str(ratio)
    
    return ratio


def planks_law(wl,T):
    """The planks law function computes the spectral density of EM radiation
    emitted by a black body in thermal equilibrium at a given temperature T.
    
    Parameters
    ----------	
    wl : float
        The wavlength in meters
    T  : float
        The temperature in Kelvin

    Returns
    -------
    B : float
        The spectral radiance of the body (amount of energy it gives offf as 
        radiation of different frequencies). Units are Ws x sr-1 x m-3 (when
        using wavelength) 
    """
    
    h = const.h
    c = const.c
    k = const.k
    
    p = ((2*h*c**2)/wl**5)*(1/(np.exp((h*c)/(wl*k*T))-1))  #planks law for wl
    return p


#The next two functions deal with the photoelectric effect, maybe move to 
#another module
def KE_escaped_electron(wavelength, p):
    """Compute the escape kinetic energy of an electron.  This is the photon
    energy minus the escape energy of the material (p).  P_gold is 7.68xe-19J
    Wavelength in microns"""
    
    ke = (const.h * const.c)/(wavelength/1e6) - p
    return ke
    
    
def escape_threshold(p):
    """compute the threshold wavelength to release electrons from the material.   
    Returns wavelength in microns.
    p is escape energy of the material (Joules) (p_gold is 7.68e-19J)
    """
    
    wavelength = 1e6*(const.h * const.c)/p
    return wavelength
    

#Time of flight - for laser ranging applications    
def TOF_to_distance(t):
    """Convert time of flight to distance (in meters) for Lidar
    """
    dist = const.c * t/2
    return dist   


if __name__ == '__main__':
    #do some test cases
    
    

