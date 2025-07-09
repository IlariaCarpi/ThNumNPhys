# Variational Monte Carlo

## Overview and Project Objectives

In this project, the variational Monte Carlo method is used in order to find the ground states energies of a 1D harmonic oscillator and a He4 nucleus. For each system, the idea is to start from an ansatz for the wave function (trial wave function) and vary the parameters using the Monte Carlo integration to calculate the energy. To vary the parameters we perform a random walk, and each set of parameters is evaluated using the Metropolis algorithm. 

The references for this project are: 

[1] R. Guardiola, "Monte Carlo Methods in Quantum Many Body Theories", 1998, DOI: 10.1007/BFb0104529

[2] R. F. Bishop, E. Buendia, M. F. Flynn and R. Guardiola, "Diffusion Monte Carlo Determination of the Binding Energy of the 4He Nucleus for Model Wigner Potentials", 1992, DOI: 10.1088/0954-3899/18/2/002

[3]  A. Gezerlis, "Numerical Methods in Physics with Python", ch. 7, 2020, DOI: 10.1017/9781108772310


## 1D Harmonic Oscillator: Project Structure

### Libraries

This Python code implements the following libraries: 


| Library                   | Purpose                                                                 |
|--------------------------|-------------------------------------------------------------------------|
| `numpy`                  | Numerical computing with support for multidimensional arrays and mathematical functions                   |
| `math`                   | Implementation of basic mathematical functions |
| `matplotlib`             | Plotting results                 |
| `os`                     | Provides functions for interacting with the operating system, such as file and directory management         |
| `pandas`                 | Organization of output in efficient and readable data structures   |


### Files Initialization

This section initializes the folders and file paths used to store output data and generated figures. Using the os library, it creates a main directory (HOResults) with two subdirectories: one for plots (Figures) and one for saved data (Data). 


### Variables Initialization

In this section, we initialize the variables and arrays that will be used throughout the code. We also set the final goal for the parameter to be \alpha = 1. 

    NAcc = 0 # accepted moves
    NP = 900  # number of moves to find parameters
    NX = 10000  # number of particles' moves
    step = 1.0  # algorithm step for random walk of particles
    step1 = 0.3 # step for parameter random walk

list_x = [] # storing positions
list_alpha = [] # storing alpha values
list_E = [] # storing energies
list_EvalE = [] 
list_variance2 = []

BestPar = 1 # final goal
Par = 0.4 # starting point for parameter
Par_new = 0 # initialize intermediate parameter for random walk

