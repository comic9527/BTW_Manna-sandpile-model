# Sandpile model: BTW, Manna and Energy height
Our software provides a general platform, capable of simulating all aspects of SOC dynamics, including driving field, stability criteria, toppling rule, boundary conditions, etc.
Plus, capability of data analysis, such as quantification of critical exponents, exertion of scaling relation, and data interpretation, such as power-law fitting and universality class identification, are all provided. 

Main Programs

BTW_Manna.m : Includes both the Bak-Tang-Wiesenfeld (BTW) and Manna models of self-organized criticality (SOC).This program allows users to set a series of model parameters such as size and number of steps.
Specific features include: simulation execution, data saving to local storage, visualization of results, and simple data fitting.

Continuous_energy.m : Includes the Zhang model, which uses continuous energy variables to simulate SOC phenomena.
Users can set model parameters such as size and number of steps.Similar features include: running simulations, saving data to local storage, visualizing results, and simple data fitting.

Utility Functions

axchange.m : Adjusts the graph's axes, allowing switching between linear and logarithmic displays.

cal_peri.m : Calculates the perimeter of clusters.

fft_bin.m : Filters power spectrum data for easier visualization.

binning.m : Simplifies data for easier visualization.

whisto.m : Calculates the proportional distribution of physical quantities.

line_fit.m & s_linefit.m & w_line_fit.m : Applies suitable linear fitting methods to different datasets.
