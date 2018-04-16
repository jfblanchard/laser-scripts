# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:14:07 2015

@author: jon.f.blanchard
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import gaussian1D_profile as gp

xc = 50           #peak center location - look into integers vs. floats  is this pixels or what?
yc = 30
d0x = 20.0         #1/e**2 diameter also known as the 2 x sigma
d0y = 10.0         #1/e**2 diameter

idx = np.arange(100)
idy = np.arange(124)

x = np.e**(-2*np.power((idx-xc)/d0x,2))
y = np.e**(-2*np.power((idy-yc)/d0y,2))

x = gp.gaussian_1D_profile(idx,xc,d0x,1)
y = gp.gaussian_1D_profile(idy,yc,d0y,1)

plt.plot(x, label='x')
plt.plot(y, label='y')
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

plt.figure()
plt.imshow(Z)
plt.contour(Z)
plt.title('2D Gaussian Profile')
plt.show()

#works

if __name__ == '__main__':
    print('In Main')
    