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

In this section, we initialize the variables and arrays that will be used throughout the code. We also set the final goal for the parameter to be α = 1. 

- Accepted particle moves: NAcc = 0
- Number of moves to find parameter: NP
- Number of particle moves: NX
- Step for particles random walk: step
- Step for parameter random walk: step1

- List for particle positions: list_x = []
- List for accepted α values: list_alpha = []
- List for accepted energies: list_E = []
- Auxiliary list of energies to be used in parameter random walk: list_EvalE = [] 
- List of squared variances of accepted energies: list_variance2 = []

- Best α value: BestPar = 1
- Starting point for parameter: Par
- Initialize intermediate parameter for random walk: Par_new = 0

- Set seed for random numbers, to make the code reproducible: np.random.seed(12231) 

#### IN REALTà SONO QUESTI PER 1DHO, GLI ALTRI SONO PER HE4 MI SONO CONFUSA
NAcc = 0 # accepted moves
NP = 900  # number of moves to find parameters
NX = 10000  # number of particles' moves
step = 1.0  # algorithm step for random walk of particles
step1 = 0.3 # step for parameter random walk

list_x = [] # storing positions
list_alpha = [] # storing alpha values
list_E = [] # storing energies
list_EvalE = [] # NON HO CAPITO A COSA SERVE QUEST'ALTRA LISTA
list_variance2 = []

BestPar = 1 # final goal
Par = 0.4 # starting point for parameter
Par_new = 0 # initialize intermediate parameter for random walk

np.random.seed(12231) 


### Important Functions

We define some functions that will be used throughout the code: 

- Returns a Gaussian wave function where α is a variational parameter that controls the width
  
  `def WF(x,alpha):`
       
  `    return exp(-0.5*alpha*alpha*x*x)`

-  Computes the local energy corresponding to the wavefunction
  
  `def LocalEnergy(x,alpha):`
  
  `    return 0.5*x*x*(1-alpha**4) + 0.5*alpha*alpha`


### Main Steps

1. **Initialize the system and generate initial configuration:**
   
   A random initial position `xOld` is generated and stored in `list_x`

2. **Run a Metropolis algorithm for position sampling:** 

   The Metropolis algorithm is used to build a Markov chain of accepted positions. Each proposed move `xNew` is accepted with probability based on the square of the ratio of the wavefunction values

3. **Perform Variational Monte Carlo (VMC) for α optimization:**

   3.1 A random new value `Par_new` for the variational parameter α is proposed within bounds

   3.2 For each sampled configuration, the local energy is computed using correlated sampling

   3.3 The average energy `E_mean` and its variance `Var2` are calculated for the proposed α

   3.4 Use a Metropolis-like criterion for accepting α: the new parameter is accepted if it lowers the energy or with a probability based on the energy ratio


### Show Results

1. **Print final results**

The code prints the last sampled values of the variational parameter α, the corresponding energy, and variance. It also prints the minimum energy found during the optimization along with the corresponding α and variance

2. **Plot energy and variance vs. α**

Two subplots are created:
* The first shows the energy as a function of α
* The second shows the variance squared as a function of α
These plots help visualize the convergence and behavior of the variational optimization

3. **Save the figure**

The plot is saved as an image file named `VMCHarmonic.png` using the `save_fig` function

4. **Create and print a DataFrame summary**

A Pandas DataFrame is created from the lists of α values, energies, and variances, and printed as a tabular summary of the results


