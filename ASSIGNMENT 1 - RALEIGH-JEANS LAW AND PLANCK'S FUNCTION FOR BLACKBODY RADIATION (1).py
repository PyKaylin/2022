#!/usr/bin/env python
# coding: utf-8

# ## **RAYLEIGH-JEANS LAW AND PLANCK'S FUNCTION FOR BLACKBODY RADIATION:**
# 
# Kaylin Shanahan 2023 
# 
# ---
# This is a python notebook which studies different elements of blackbody radiation. The program first uses the Rayleigh-jeans law to calculate the energy density of a blackbody, then plots the Rayleigh-Jeans law over a range of frequencies which demonstrates the ultraviolet catastrophe. The program then uses Planck's function, plotting Planck's function over a range of temperatures. Finally, the program performs numerical integration of the Planck function, using the Rectangle and Simpson's rule, to verify the area under the curve is proportional to $T^4$.
# 
# ---

# ---
# # Task 1:
# ---
# The Rayleigh-Jeans law for the spectrum of a blackbody that will be used in this programme is:
# $$ U(v, T) = \frac {8 \pi v^2}{c^3} kT $$
# 
# Where $ k = 1.38 \times 10^{-23} J K^{-1} $ and $ c = 3.0 \times 10^8 ms^{-1}. $
# 
# ---
# 
# * Input constants and user values for frequency and temperature
# * Calculate $U(v, T)$ using Rayleigh-Jeans law formula described above
# * Ouput value of $U(v, T)$
# 
# ---

# In[1]:


# TASK 1
# A program to print the value of U(v, T) for a user inputted frequency and temperature.
# Tested with v = 1.5e12 Hz and T = 2500K

# Numpy is needed for calculations
import numpy as np 

# Input Constants
k = 1.38e-23   # Boltzmann constant in J/K
c = 3.0e8      # Speed of light in m/s
pi = np.pi     # Pi constant

# User input for frequency and temperature
v = eval(input("Enter the frequency in Hertz: "))
T = eval(input("Enter the temperature in Kelvin: "))

# Calculate Energy Density value U(v, T)
U = (((8) * (pi) * (v**2)) / (c**3)) * (k * T) # Raleigh-Jeans law

# Output Energy Density value U(v, T)
print("U(v, T) = {0:10.2e}".format(U))


# ---
# # Task 2:
# ---
# 
# The plot of the Rayleigh-Jeans law over a range of frequencies below shows how the Rayleigh-Jeans law leads to the Ultraviolet Catastrophe.
# 
# ---
# 
# - Input values for temperature and frequencies
# - Calculate energy density using the Rayleigh-Jeans Law described above
# - Output plot of The Rayleigh-Jeans Law
# 
# ---

# In[15]:


# TASK 2
# Program to plot the Rayleigh-Jeans law over a range of frequencies 
# to show that it leads to the Ultraviolet Catastrophe

# Matplotlib is needed for plot of the Rayleigh-Jeans law
import matplotlib.pyplot as plt 

# Input values for T, F0, Ff
T = 2500     # Temperature in Kelvin
F0 = 10e11   # Initial Frequency in Hz
Ff = 10e15   # Final Frequency in Hz

# Create range of frequencies
F_Range = np.linspace(F0, Ff)

# Calculate Energy Density
U =  (((8) * (pi) * (F_Range**2)) / (c**3)) * (k * T) # Rayleigh-Jeans Law

# Plot The Rayleigh-Jeans Law
plt.figure(figsize=(10,5))   # Make the plot rectangular
plt.semilogx(F_Range, U, linestyle="--", label="Rayleigh-Jeans Law (T = 2500 K)")   # Make log plot 
plt.title("Rayleigh-Jeans Law")   # Title the plot
plt.xlabel("Frequency ($Hz$)")   # Label x-axis
plt.ylabel("Energy Density ($W/m^2/Hz$)")   # Label y-axis
plt.legend()   # Create legend for graph
plt.grid()   # Create grid on graph
plt.show()   # Display graph


# ---
# 
# The graph above illustrates the ultraviolet catastrophe, as the energy density increases as the square of the frequency, leading to the infinitely increasing graph seen above.
# 
# ---

# ---
# # Task 3:
# ---
# 
# Planck derived a different formula for blackbody radiation, known as Planck's function: 
# $$ u(v, T) = \frac {8 \pi v^2}{c^3} \frac {1}{e^{hv/kT}-1} $$
# 
# This formula contains a new constant, $ h = 6.63 \times 10^{−34} J s $. 
# This can be written in wavelength coordinates as follows:
# $$ u(λ, T) = \frac {8 \pi hc}{λ^5} \frac {1}{e^{hc/λkT}-1} $$
# 
# ---
# 
# Planck's function solved the problem of the Ultaviolet Catastrophe, as it does not diverge from observed intensity a high frequencies like the Rayleigh-Jeans law - The plot below shows the wavelength version of Planck's function for temperatures between $2000 - 10000K$.
# 
# ---
# 
# - Input Plack's constant and wavelength and temperature ranges
# - Calculate energy density using Planck's function described aboce
# - Output plot of Planck's function for range of temperatures
# 
# ---

# In[345]:


# TASK 3
# Plot the wavelength version of Planck’s function for a few temperatures between 2000 and 10000K

# Input new constant
h = 6.63e-34   # Planck's constant in Js

# Create range of wavelengths
W_Range = np.linspace(11e-9, 3000e-9, 1000)   # Wavelength in meters

# Create range of temperatures
T_Range = np.linspace(2000, 10000, 9)   # Temperature in Kelvin

# Create the plot
plt.figure(figsize=(10,5))

# Calculate and plot energy density using Planck's Function for each temperature
for T in T_Range:   # Use a for loop for range of temperatures
    u = ((8 * np.pi * h * c) / ((W_Range) ** 5)) * (1 / (np.exp(h * c / ((W_Range) * k * T)) - 1))   # Calculate energy density using Planck's Function
    plt.semilogx(W_Range, u, linestyle="--", label=f"T = {T} K")   # Make log graph
    
# Plot Planck's Function (wavelength version)
plt.title("Planck's Function")   # Title the function
plt.xlabel("Wavelength ($m$)")   # Label the x-axis
plt.ylabel("Energy Density ($W/m^2/Hz$)")   # Label the y-axis
plt.legend()   # Make a legend for the plot
plt.grid()   # Create grid on graph
plt.show()   # Display graph


# ---
# # Task 4:
# ---
# The peak of the Planck function for each temperature can now be found. From this information, the constant in Wien's displacement law can be evaluated. This gives the peak of the Planck function:
# $$ \lambda_{peak} T = Constant $$
# ---
# 
# - Input list to store Wien's Displacement Law Constants
# - Calculate the peak wavelength of Planck's function
# - Calculate Wien's Displacement Law Constant for each temperature
# - Output the value of Wien's Displacement Law Constant
# 
# ---

# In[343]:


# TASK 4
# Find the peak of the Planck function for each temperature.
# Evaluate the constant in Wien’s displacement law which gives the peak of the Planck function

# Initialise list of Wien's displacement law constants
W_dcs = []

# Calculate the energy density using Planck's function
for T in T_Range:
    u = ((8 * np.pi * h * c) / (W_Range ** 5)) * (1 / (np.exp(h * c / (W_Range * k * T)) - 1))
    
    # Calculate the peak wavelength
    Wavelength_Peak = W_Range[np.argmax(u)]
    
    # Calculate the Wien's displacement law constant
    Constant = Wavelength_Peak * T
    # Store Wien's displacement law constant in list
    W_dcs.append(Constant)

# Calculate average of Wien's displacement constant for most accurate result
W_dcs_avg = np.mean(W_dcs)   # Mean of Wien's displacement constant values

# Output the final result for the Wien's displacement constant
print("Wien's Displacement Law Constant: {0:10.3e} m.k".format(W_dcs_avg))


# ---
# # Task 5:
# ---
# Finally, a numerical integration of the Planck function for a
# few temperatures will be carried out, to verify that the area under the curve is proportional to $T^4$ as in the Stefan-Boltzmann law. 
# 
# The Rectangle and Simpson Rules will be used as the two different integration methods to analyse the Planck function
# 
# ---
# 
# - Input new constants and lists
# - Calculate area under curve and stefan boltzmann constant
# - Output area under curve, stefan boltzmann constant and the proportionality between the two
# 
# ---

# ---
# # Rectangle Rule:
# - The Rectangle rule is the simplest method of integration. The rule approximates the integral of a function, or the area under a curve, by adding the areas of rectangles of width $Δx$ and height of the value of the function $f(x)$ at the midpoint of the length.
# - It can be described with the following formula:
# $$ ∫^{b}_{a} f(x) dx = Δx(f(x_0) + f(x_1) + f(x_2)+ ... + f(x_n-1) + f(x_n)) $$
# where $ Δx = \frac {b - a}{n} $
# ---

# In[346]:


# TASK 5.1
# Numerical integration of the Planck function for a few temperatures using
# the Rectangle rule to verify that the area under the curve is proportional to T4 as in the Stefan-Boltzmann law.

# Input constants
sigma = 5.67e-8 # Stefan Boltzmann constant in W/m^2/K^4

# Initialise list of areas
Rec_A_i = []   # Rectangle areas

# Numerical integration of the Planck function for range of temperatures
for T in T_Range:
    u = ((8 * np.pi * h * c) / (W_Range ** 5)) * (1 / (np.exp(h * c / (W_Range * k * T)) - 1))   # Planck function
    
    # Calculate width
    W_x = W_Range[1] - W_Range[0]
    
    # Calculate area of rectangle
    Rec_A = W_x * (np.sum(u))
    
    # Store area in list
    Rec_A_i.append(Rec_A)
    
    # Sum of rectangle areas
    Rec_A_c = np.sum(Rec_A_i)      # Add all area values in list
    Rec_Area_r = round(Rec_A_c, 2)   # Round area value
    
    # Verify area is proportional to T^4 in Stefan Boltzmann Law
    Rec_SB = sigma * (T**4)
    Rec_SB_r = round(Rec_SB)   # Round Stefan Boltzann value
    
    # Round temperature values
    T_r = round(T)
    
    # Verification proportionality of Area under curve and Stefan Boltzmann law value
    print("Stefan Boltzmann Law Value when T = {} K: {} W/m^2".format(T_r, Rec_SB_r)) 

# Output area under curve
print("Total area under curve of Planck's function: {} m^2".format(Rec_Area_r))


# ---
# 
# From the above values, it can be verified that the area under the curve is proportional to $T^4$ as in the Stefan-Boltzmann law.
# 
# The calculated values could be verified with the expected values calculated by hand. 
# 
# ---

# ---
# # Simpson's Rule:
# - Simpson's Rule is significantly more accurate than other methods of integration, such as the Trapezoidal rule and Rectangle rule. This method uses parabolas at the top of small sections of width $Δx$ instead of straight lines, as used in Rectancle rule above. The parabolas are a much closer fit to the curve providing a more accurate answer (though not exact).
# - It can be described with the following formula:
# $$ ∫^{b}_{a} f(x) dx = \frac {Δx}{3} \times (f(x_0) + 4f(x_1) + 2f(x_2) + ... + 4f(x_{n−1}) + f(x_n)) $$
# ---

# In[340]:


# TASK 5.2
# Numerical integration of the Planck function for a few temperatures using
# Simpson's rule to verify that the area under the curve is proportional to T4 as in the Stefan-Boltzmann law.

# Initialise list of areas
Simp_A_i = []   # Simpson's rule areas

# Numerical integration of the Planck function for range of temperatures
for T in T_Range:
    u = ((8 * np.pi * h * c) / (W_Range ** 5)) * (1 / (np.exp(h * c / (W_Range * k * T)) - 1))
    
    # Calculate area using Simpson's rule
    n = 1000   # Number of steps
    Simp_A = 0
    for i in range(0, n-2, 2):
        Simp_A += (W_x / 3) * (u[i] + 4*u[i+1] + u[i+2])   # sing formula of Simpson's rule
    
    # Store area in list
    Simp_A_i.append(Simp_A)
    
    # Sum of areas
    Simp_A_c = np.sum(Simp_A_i)   # Add all area values in list
    Simp_Area_r = round(Simp_A_c, 2)    # Round area value
    
    # Verify area is proportional to T^4 in Stefan Boltzmann Law
    Simp_SB = sigma * (T**4)
    Simp_SB_r = round(Simp_SB)   # Round Stefan Boltzann value
    
    # Verification proportionality of Area under curve and Stefan Boltzmann law value
    print("Stefan Boltzmann Law Value when T = {} K: {} W/m^2".format(T_r, Simp_SB_r)) 

# Output area under curve
print("Total area under curve of Planck's function: {} m^2".format(Simp_Area_r))


# ---
# 
# From the above values, it can be seen that the area calculated is actually less than the rectangle rule. It could be said, that as the Simpson rule uses the parabola instead of the recatangle, and so can be more accurate, that this value is the more accurate of the two.
# 
# However, it must be noted that the Simpson's rule becomes less accurate for narrow peak-like functions.
# 
# As with the Rectangle rule, from the above values, it can be verified that the area under the curve is proportional to $T^4$ as in the Stefan-Boltzmann law.
# 
# The calculated values could be verified with the expected values calculated by hand.
# 
# ---
