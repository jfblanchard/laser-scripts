# laser-tools

Utilities for performing common laser beam computations.

This module has a several functions for generating and plotting 1D and 2D gaussian profiles, as well as other utilities which can be used for modeling and computing laser beam parameters.  Some of the functions include:

### Beam parameter related:

  - *bpp_radius_halfdiv*   - Compute the beam parameter product using radius and half divergence
  - *bpp_diam_fulldiv*     - Compute the beam parameter product using diameter and full divergence
  - *bpp_wl_m2*            - Compute the beam parameter product from wavelength and m-squared
  - *beam_half_divergence* - Compute the half angle divergence from wavelength, waist radius, and m-squared value.
  - *beam_waist*           - Compute the beam waist from wavelength, divergence, and m-squared value
  - *rayleigh_range*       - Compute the rayleigh range of a gaussian beam in free space
  - *beam_size*            - Compute the beam size (radius) at location z
  - *q_parameter*          - Compute the q-parameter from wavelength, waist radius, and z location
  
### Power related:
- *peak_power*
- *intensity*
- *peak_intensity*
- *fluence*
- *peak_fluence*
- *beamsize_power_through_app* 

### Other
- *mode_frequency_separation*
- *freq_to_wl*
- *wl_to_freq*
- *wl_to_wavenumber*
- *photon_energy*
- *emitted_frequency*
- *boltzman_statistics*
- *planks_law*
- *KE_escaped_electron*
- *escape_threshold*
- *TOF_to_distance*

<br>
<h3>Gaussian Profiles Usage<//h3>

<h4>Gaussian 1D Profile</h4>

Usage (running in jupyter console): 
```python
%run gaussian1D_profile.py
```
will produce the following plot:

![1D Gaussian Image](/images/gaussian1D_image.png)

<h4>Gaussian 2D Profile</h4>

Usage (running in jupyter console): 
```python
%run gaussian2D_profile.py
```
will produce the following plots:

![2D Gaussian lines](/images/gaussian2D_lines.png)
![2D Gaussian profile](/images/gaussian2D_image.png)


