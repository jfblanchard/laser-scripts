<h3>Usage of gaussian-beam package <//h3>

<h4>gaussian1D_profile.py</h4>

Usage: 
```python
import gaussian1D_profile as gp
x,y = gp.gaussian_1D_profile(-50,50,.2, 0, 10, 1)
gp.plot_1d_gaussian(x,y,True)
```

will produce the following plot:

![1D Gaussian Image](/images/gaussian1D_image.png)
