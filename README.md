# Variational Monte Carlo

## Overview and Project Objectives

In this project, the variational Monte Carlo method is used in order to find the ground states energies of a 1D harmonic oscillator and a He4 nucleus. For each system, the idea is to start from an ansatz for the wave function (trial wave function) and vary the parameters using the Monte Carlo integration to calculate the energy. To vary the parameters we perform a random walk, and each set of parameters is evaluated using the Metropolis algorithm. 

The references for this project are: 

[1] R. Guardiola, "Monte Carlo Methods in Quantum Many Body Theories", 1998, DOI: 10.1007/BFb0104529

[2] R. F. Bishop, E. Buendia, M. F. Flynn and R. Guardiola, "Diffusion Monte Carlo Determination of the Binding Energy of the 4He Nucleus for Model Wigner Potentials", 1992, DOI: 10.1088/0954-3899/18/2/002

[3]  A. Gezerlis, "Numerical Methods in Physics with Python", ch. 7, 2020, DOI: 10.1017/9781108772310


-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

- Best α value: BestPar = 1
- Starting point for parameter: Par

- Set seed for random numbers, to make the code reproducible: np.random.seed(12231) 


### Important Functions

We define some functions that will be used throughout the code: 

- Returns a Gaussian wave function where α is a variational parameter that controls the width
  
  `def WF(x,alpha):
    return exp(-0.5*alpha*alpha*x*x)`

-  Computes the local energy corresponding to the wavefunction
  
  `def LocalEnergy(x,alpha):
     return 0.5*x*x*(1-alpha**4) + 0.5*alpha*alpha`


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



### Results 

edrtfhui



-----------------------------------------------------------------------------------------------------------------------------------------------------------



## Helium 4: Project Structure

### Libraries

This Python code implements the following libraries: 


| Library                   | Purpose                                                                 |
|--------------------------|-------------------------------------------------------------------------|
| `numpy`                  | Numerical computing with support for multidimensional arrays and mathematical functions                   |
| `matplotlib`             | Plotting results                 |
| `os`                     | Provides functions for interacting with the operating system, such as file and directory management         |
| `pandas`                 | Organization of output in efficient and readable data structures   |
| `numba`                 | Acceleration of numerical Python code by translating functions to optimized machine code at runtime   |

### Files Initialization

This section initializes the folders and file paths used to store output data and generated figures. Using the os library, it creates a main directory (4HEResults) with two subdirectories: one for plots (Figures) and one for saved data (Data). 

### Variables Initialization

In this section, we initialize the variables and arrays that will be used throughout the code. We also set the final goal for the parameters as in reference [1]. 

- Accepted particle moves: NA = 0
- Number of moves to find parameter: NP
- Number of particle moves: Nm
- Number of nucleons: A
- Step for particles random walk: step
- Step for γ random walk: step1
- Step for "a" random walk: step2
- Step for β random walk: step3

- Planck constant: hbar = 6.582119*10**(-22)
- Kinetic term factor: C = (197.3269804)**2/2/939.56542052
- Differential increment: h = 0.00001

- Best values for parameters: BestPars = [0.08597, -0.7191, 2.13796] (Guardiola, ref. [1])
- Starting point for parameters: Pars

- Set seed for random numbers, to make the code reproducible: np.random.seed(12231)


### Important Functions

We define some functions that will be used throughout the code: 

- Returns the square distance between two particles i and j
  
  `def d(R,i,j):
    return np.sum((R[:, i] - R[:, j])**2)`

-  Computes the trial wave function as a product of Gaussian two-body terms 
  
  `def WF(R, Par):
    wf = 1
    for i in range(0,A-1):
      for j in range(i+1,A):
        wf = wf * (exp(-Par[0]*d(R,i,j))+Par[1]*exp(-Par[2]*d(R,i,j)))
    return wf`

-  Calculates the total potential energy between nucleon pairs
  
  `def V(R):
    V = 0
    for i in range(0,A-1):
      for j in range(i+1,A):
        V = V + 1000*exp(-3*d(R,i,j))-165.35*exp(-1.05*d(R,i,j))-21.5*exp(-0.6*d(R,i,j))-83*exp(-0.8*d(R,i,j))-11.5*exp(-0.4*d(R,i,j))
    return V`

-  Returns an approximate second derivative of the wavefunction for the kinetic energy
  
  `def Diff (R, i, j, h, Par):
    R_P = R.copy()
    R_M = R.copy()
    R_P[i,j] = R[i,j]+h
    R_M[i,j] = R[i,j]-h
    return (WF(R_P, Par)+WF(R_M, Par)-2*WF(R, Par))/h**2`

-  Computes the kinetic energy by summing second derivatives over all coordinates
  
  `def K(R, Par):
    K = 0 
    for j in range(0,A):
      for i in range(0,3):
        K = K + Diff(R, i, j, h, Par )*C
    return K`



### Main Steps

1. **Initialize the system and generate initial configuration:**
   
   A random initial position `xOld` is generated and stored in `list_x`

2. **Metropolis sampling of coordinates:** 

   Perform a Metropolis random walk of the nucleon coordinates using the reference (target) parameters. Accepted configurations are stored in `list_R`

3. **Precompute fixed values:**

   To speed up the code execution, precompute the wavefunction (WF) with the target parameters and the potential energy (V) for all sampled configurations

4. **Variational Monte Carlo for parameters' optimization:**


   4.1 Random new values for variational parameters are proposed (`Par_new[0]`, `Par_new[1]` and `Par_new[2]`)

   4.2 For each sampled configuration, the local energy is computed using correlated sampling

   4.3 The average energy `E_mean` and its variance `Var2` are calculated for the proposed set of parameters

   4.4 Use a Metropolis-like criterion for accepting `Par_new`: the new parameters are accepted if it lowers the energy or with a probability based on the energy ratio


### Show Results

1. **Print final results**

The code prints the last sampled values of the variational parameters, and the corresponding energy and variance. It also prints the minimum energy `minimumE` found during the optimization along with the corresponding γ (`minimumGamma`), a (`minimumA`) and β (`minimumBeta`), and variance `sqrt(minimumVar2)`


2. **Visualize the energy landscape:**

A 3D scatter plot is generated showing how the energy depends on the parameters γ, a, and β. Each point represents a tested parameter set, and its color indicates the corresponding energy; this helps identify the region where the energy minimum occurs.


3. **Save the figure**

The plot is saved as an image file named `VMC4He.png` using the `save_fig` function
