#!/usr/bin/env python
# coding: utf-8

# ## **TWIN SLIT EXPERIMENT: DIFFRACTION AND INTERFERENCE**
# 
# Kaylin Shanahan 2023 
# 
# ---
# 
# This is a python notebook which models the two-slit experiment. The following tasks will model 1D intensity patterns for a single slit and a twin slit, then a 2D Intensity pattern, first with just the horizontal, and then as would be viewed on a viewing card of the twin slit experiment.
# 
# The formula used in this notebook is called the Twin Slit formula:
# $$ I = I_{0} \big[\frac {sin(\alpha)} {\alpha}\big]^{2} [cos(\beta)]^{2} $$
# 
# Where, $I_{0}$ = initial intensity, the width of one slit is $\alpha$ = 0.09 mm, the distance from the slits to the screen is $L$ = 480 nm, the distance between the slits is $d$ = 0.4 mm, the wavelength of light is $\lambda$ = 670 nm.
# 
# $$ \alpha = \frac {\pi \alpha} {\lambda} sin \theta $$
# $$ \beta = \frac {\pi d} {\lambda} sin \theta $$
# 
# It is used to give the horizontal variation of the intensity of light at a screen after passing through a double slit.
# 
# ---

# ---
# # Task 1:
# ---
# 
# This program plots the predicted 1-D intensity pattern $ \frac {I}{I_{0}} $ of light observed from two slits. Displayed will be a pattern of peaks and trough as a function of $y$ in milimetres. It should be noted that $ tan \theta = \frac {x}{L} $, and the z-axis is the direction from the slit to the screen.
# 
# ---
# 
# * Input constants and array
# * Calculate theta, alpha, beta and intensity
# * Ouput plot of the 1D intensity pattern of two slits
# 
# ---

# In[89]:


# Numpy is needed for calculations
import numpy as np

# Matplotlib is needed for plot
import matplotlib.pyplot as plt 

# Input constants
a = 0.09   # width of slit in milimetres
L = 480   # distance from slits in milimetres
d = 0.4   # distance between slits in milimetres
wavelength = 0.00067   # Laser wavelength in milimetres

# Create array of x-values
x = np.linspace(-7, 7, 400)

# Calculate theta
theta = np.arctan(x/L)

# Calculate alpha
alpha = ((np.pi * a) / wavelength)*(np.sin(theta))

# Calculate beta
beta = ((np.pi * d) / wavelength)*(np.sin(theta))

# Calculate I/I0
Intensity_2 = ( (np.sin(alpha)/alpha) **2 ) * ( (np.cos(beta)) **2 )

# Plot the 1D intensity pattern of two slits
plt.plot(x, Intensity_2)   # Create plot
plt.title("1D Intensity Pattern of Light From Two Slits")   # Title the plot
plt.xlabel("Detector-slit Position (mm)")   # Label x-axis
plt.ylabel("$I$ / $I_0$")   # Label y-axis
plt.grid()   # Create grid on graph
plt.show()   # Display graph


# ---
# # Task 2:
# ---
# 
# The following program plots the predicted 1-D intensity pattern of light
# observed from just one of the slits, along with the 1-D intensity pattern of light observed from two slits.
# 
# ---
# 
# * Input constants as in task 1
# * Calculate intensity of one slit and two slits
# * Ouput plot of the 1D intensity pattern of two slits and one slit
# 
# ---

# In[90]:


# Calculate I/I0 for two slits
Intensity_2 = ( (np.sin(alpha)/alpha) **2 ) * ( (np.cos(beta)) **2 )

# Calculate I/I0 for one slit
Intensity_1 = ( (np.sin(alpha)/alpha) **2 )

# Plot the 1D intensity pattern of two slits and one slit
plt.plot(x, Intensity_2, linestyle="-", label="Two Slits")   # Create plot
plt.plot(x, Intensity_1, linestyle="--", color="red", label="One Slit")   # Create plot
plt.title("1D Intensity Pattern of Light From One Slit and Two Slits")   # Title the plot
plt.xlabel("Detector-slit Position (mm)")   # Label x-axis
plt.ylabel("$I/I0$")   # Label y-axis
plt.legend()   # Create legend on graph
plt.grid()   # Create grid on graph
plt.show()   # Display graph


# ---
# # Task 3:
# ---
# 
# 
# The following program simulates the appearance of the pattern in 2-D
# on the viewing card as in the experiment. In this case it is assumed that the 1-D pattern is repeated in vertical direction.
# 
# ---
# 
# * Input arrays
# * Calculate theta, alpha, beta and intensity
# * Ouput plot of 2D intensity pattern
# 
# ---

# In[91]:


# Create array of x-values
x = np.linspace(-3, 3, 400)   # Adjust x range, with 3:2 ratio

# Create array of y-values
y = np.linspace(-2, 2, 400)   # Adjust y range, with 3:2 ratio

# Create grid of x and y values
X, Y = np.meshgrid(x, y)

# Create array for the 2D intensity pattern
Intensity_2D_two = np.zeros((len(y), len(x)))

# Calculate 2D intensity pattern
for i in range(len(x)):
    for j in range(len(y)):
        # Calculate theta
        theta = np.arctan(x[i]/L)

        # Calculate alpha
        alpha = ((np.pi * a) / wavelength)*(np.sin(theta))

        # Calculate beta
        beta = ((np.pi * d) / wavelength)*(np.sin(theta))
        
        # Calculate intensity for two slits
        Intensity_two = ( (np.sin(alpha)/alpha) **2 ) * ( (np.cos(beta)) **2 )
        Intensity_2D_two[j, i] = Intensity_two

# Calculate minimum and maximum intensity for levels
Imin = np.min(Intensity_2D_two)   # Minimum Intensity
Imax = np.max(Intensity_2D_two)   # Maximum Intensity

# Set levels for contour plot
levels = np.linspace(Imin, Imax, 400)

# Create 2D Intensity Pattern Plot
plt.contourf(X, Y, Intensity_2D_two, levels=levels, cmap="gist_heat")   # Create contour plot
plt.colorbar(label='I/I0')   # Create colorbar
plt.title("2D Intensity Pattern of Light From Twin Slits (Without fall-off)")   # Title the plot
plt.xlabel('Screen Width (mm)')   # Label x-axis
plt.ylabel('Screen Height (mm)')   # Label y-axis
plt.show()   # Display graph


# ---
# # Task 4:
# ---
# 
# The following program plots the actual 2D intensity pattern of light, taking into consideration the variation of the light observed in the y-direction is that it falls off away from the axis. The final 2-D plot looks like one that would be observed in the actual two-slit experiment.
# 
# ---
# 
# * Input constant
# * Calculate theta, alpha, beta in x and y-directions
# * Calculate intensity
# * Output plot of 2D intensity pattern
# 
# ---

# In[92]:


# Input constants
b = 1.0   # height of slit in milimeters

# Create array for the 2D intensity pattern
Intensity_2D_xy = np.zeros((len(y), len(x)))

# Calculate 2D intensity pattern
for i in range(len(x)):
    for j in range(len(y)):
        # Calculate in x direction:
        # Calculate theta
        theta_x = np.arctan(X[j, i] / L)

        # Calculate alpha
        alpha_x = ((np.pi * a) / wavelength) * (np.sin(theta_x))

        # Calculate beta
        beta_x = ((np.pi * d) / wavelength) * (np.sin(theta_x))
        
        # Calculate intensity for two slits in x-direction
        Intensity_x = ( (np.sin(alpha_x)/alpha_x) ** 2 ) * ( (np.cos(beta_x)) ** 2 )
        
        # Calculate in y direction:
        # Calculate theta
        theta_y = np.arctan(Y[j, i] / L)

        # Calculate alpha
        alpha_y = ((np.pi * b) / wavelength)*(np.sin(theta_y))

        # Calculate intensity for two slits in y-direction
        Intensity_y = ( (np.sin(alpha_y)/alpha_y) **2 )
        
        # Calculate intensity for two slits in x and y-direction
        Intensity_2D_xy[j, i] = Intensity_x * Intensity_y

# Calculate minimum and maximum intensity for levels
Imin = np.min(Intensity_2D_xy)   # Minimum Intensity
Imax = np.max(Intensity_2D_xy)   # Maximum Intensity

# Set levels for contour plot
levels = np.linspace(Imin, Imax, 400)

# Create 2D Intensity Pattern Plot
plt.contourf(X, Y, Intensity_2D_xy, levels=levels, cmap="gist_heat")   # Create contour plot
plt.colorbar(label='I/I0')   # Create colorbar
plt.title("2D Intensity Pattern of Light From Twin Slits (With fall-off)")   # Title the plot
plt.xlabel('Screen Width (mm)')   # Label x-axis
plt.ylabel('Screen Height (mm)')   # Label y-axis
plt.show()   # Display graph

