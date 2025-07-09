# Variational Monte Carlo

## Overview and Project Objectives

In this project, the variational Monte Carlo method is used in order to find the ground states energies of a 1D harmonic oscillator and a He4 nucleus. For each system, the idea is to start from an ansatz for the wave function (trial wave function) and vary the parameters using the Monte Carlo integration to calculate the energy. To vary the parameters we perform a random walk, and each set of parameters is evaluated using the Metropolis algorithm. 

The references for this project are: 

[1] R. Guardiola, "Monte Carlo Methods in Quantum Many Body Theories", 1998, DOI: 10.1007/BFb0104529

[2] R. F. Bishop, E. Buendia, M. F. Flynn and R. Guardiola, "Diffusion Monte Carlo Determination of the Binding Energy of the 4He Nucleus for Model Wigner Potentials", 1992, DOI: 10.1088/0954-3899/18/2/002

[3]  A. Gezerlis, "Numerical Methods in Physics with Python", ch. 7, 2020, DOI: 10.1017/9781108772310


## 1dHO: Project Structure

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

In the first part, we initialize the files to save the data and graphs, and we define the paths to store such files. 

We locate the files in the folder "HOResults", and we create a directory for the figures ("Figures") and one for the data ("Data"). 


PROJECT_ROOT_DIR = "HOResults" 
FIGURE_ID = "HOResults/Figures"
DATA_ID = "HOResults/Data"

if not os.path.exists(PROJECT_ROOT_DIR):
    os.mkdir(PROJECT_ROOT_DIR)
if not os.path.exists(FIGURE_ID):
    os.makedirs(FIGURE_ID)
if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)

def image_path(fig_id):
    return os.path.join(FIGURE_ID, fig_id)
def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

def save_fig(fig_id):
    plt.savefig(image_path(fig_id) + ".png", format='png')

outfile = open(data_path("Data.dat"),'w')

### Variables Initialization

We then initialize the variables and arrays that will be used throughout the code. We set the values for the total number of moves for the parameter variation and for the particles' random walk, also defining the steps used in the respective random walks. Then, we initialize the particles' positions array, the parameter values array, the ...

We also set the final goal to be \alpha = 1

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

