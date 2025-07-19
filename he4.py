import os
import numpy as np
from numpy import exp as exp
from numpy import sqrt as sqrt
from matplotlib import pyplot as plt
from numba import njit



# FILES INITIALIZATION ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

PROJECT_ROOT_DIR = "4HEResults" 
FIGURE_ID = "4HEResults/Figures"
DATA_ID = "4HEResults/Data"

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

# saving figure into designated directory 
def save_fig(fig_id):
    plt.savefig(image_path(fig_id) + ".png", format='png')

# saving data in designated file
outfile = open(data_path("Data4HE.dat"),'w')
 


# DEFINE FUNCTIONS ------------------------------------------------------------------------------------------------------------------------------

# square radius

@njit
def d(R,i,j):
  return np.sum((R[:, i] - R[:, j])**2)

# wavefunction

@njit
def WF(R, Par):
  wf = 1
  for i in range(0,A-1):
    for j in range(i+1,A):
      wf = wf * (exp(-Par[0]*d(R,i,j))+Par[1]*exp(-Par[2]*d(R,i,j)))
  return wf 


# potential

@njit
def V(R):
  V = 0
  for i in range(0,A-1):
    for j in range(i+1,A):
      V = V + 1000*exp(-3*d(R,i,j))-165.35*exp(-1.05*d(R,i,j))-21.5*exp(-0.6*d(R,i,j))-83*exp(-0.8*d(R,i,j))-11.5*exp(-0.4*d(R,i,j))
  return V

# second order differential for kinetic term

@njit
def Diff (R, i, j, h, Par):
  R_P = R.copy()
  R_M = R.copy()
  R_P[i,j] = R[i,j]+h
  R_M[i,j] = R[i,j]-h
  return (WF(R_P, Par)+WF(R_M, Par)-2*WF(R, Par))/h**2 

# kinetic term

@njit
def K(R, Par):
  K = 0 
  for j in range(0,A):
    for i in           range(0,3):
      K = K + Diff(R, i, j, h, Par )*C
  return K



# INITIALIZE VARIABLES -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

NA = 0 # accepted moves
NP = 135  # number of moves to find parameters
NM = 10000  # number of particles' moves
A = 4   # numbers of nucleons
h = 0.00001  # differential increment 
step = 0.55  # algorithm step for random walk of particles
hbar = 6.582119*10**(-22)
C = (197.3269804)**2/2/939.56542052 # constant in front of the kinetic term, the first term is in Mev*fm
step1 = 5/100 #1/100
step2 = 25/100 #5/100
step3 = 50/100 # 10/100
list_gamma = [] # storing gamma values
list_a  = [] # storing a values
list_beta = [] # storing beta values
list_E = [] # storing accepted energies
list_EvalE = [] # storing energies
list_variance2 = [] # store variances^2

# Goal: parameters in Guardiola reference: 
# gamma = 0.08597 
# a = -0.7191 
# beta = 2.13796

# parameters' initialization

BestPars = [0.08597 , -0.7191 , 2.13796 ] # final goal
Bestpars = np.array(BestPars) 

Pars = [0.2, -0.2, 1] # starting point for parameters
Par = np.array(Pars)

Par_new = np.zeros(3) # initialize intermediate parameters for random walk

# random seed for code reproducibility
np.random.seed(12231)


# first coords
R = (np.random.rand(3,A)-0.5)*10*step

list_R = []
list_R.append(R)

list_WFbest = []
list_pot = []





# find the set of steps with Metropolis algorithm --------------------------------------------------------------------------------------------------------

for k in range(0,NM):
  R_new = R.copy() + step*(np.random.rand(3,A)-0.5)
  if ((WF(R_new, Bestpars)/WF(R, Bestpars))**2>1 or (WF(R_new, Bestpars)/WF(R, Bestpars))**2 > np.random.rand() ):
    R = R_new.copy()
    list_R.append(R)
    NA = NA +1
  else:
    R = R.copy()
    list_R.append(R) 

print("Acceptance ratio =", 100*NA/NM)


# calculate in advance the values of WF with the set of goal parameters and of V

for l in range(0,NM): 
    list_WFbest.append(WF(list_R[l], Bestpars))
    list_pot.append(V(list_R[l]))






# VMC TO FIND BEST PARAMETERS ------------------------------------------------------------------------------------------------------------------------------

for l in range(0,NP):

  # random walk
  
  Par_new[0] = Par[0].copy() + (np.random.rand()-0.5)*step1
  Par_new[1] = Par[1].copy() + (np.random.rand()-0.5)*step2
  Par_new[2] = Par[2].copy() + (np.random.rand()-0.5)*step3
  E_loc = 0
  E_loc_2 = 0
  D = 0
  list_tempWF = []
  list_tempK = []
  
  # to speed up calculations
  for p in range(0,NM): 
    list_tempWF.append(WF(list_R[p], Par_new))
    list_tempK.append(K(list_R[p], Par_new))    

    
  # calculate local energy for current set
  for j in range(0,NM):
    E_loc = E_loc + list_pot[j]*(list_tempWF[j]/list_WFbest[j])**2 - list_tempK[j]*list_tempWF[j]/list_WFbest[j]**2 # Modify how the parameters arrive   
    E_loc_2 = E_loc_2 + (list_pot[j]*(list_tempWF[j]/list_WFbest[j])**2 - list_tempK[j]*list_tempWF[j]/list_WFbest[j]**2)**2
    D = D + (list_tempWF[j]/list_WFbest[j])**2
    
  E_mean = E_loc/D
  Var2 = (E_loc_2/D - (E_loc/D)**2)/D
  
  list_tempWF.clear()
  list_tempK.clear()


  # METROPOLIS QUESTION
    
  if (l==0):
    Par = Par_new.copy()  
    list_gamma.append(Par[0])
    list_a.append(Par[1])
    list_beta.append(Par[2])
    list_E.append(E_mean)
    list_EvalE.append(E_mean)
    list_variance2.append(Var2)
    minimumE = E_mean
    minimumVar2 = Var2
    minimumGamma = Par[0]
    minimumA = Par[1]
    minimumBeta = Par[2]
    print("\nFirst: ", "\n", E_mean, "\nWith pars:","\nGamma = ",list_gamma[l],"\na = ",list_a[l],"\nbeta = ",list_beta[l],"\nand variance^2: ",Var2)
  
  else:   
    print("l = ", l)
    if (list_EvalE[l-1]>E_mean or list_EvalE[l-1]/E_mean>np.random.rand()):
    # if (list_EvalE[l-1]>E_mean):
      if (minimumE>E_mean):
        minimumE = E_mean
        minimumVar2 = Var2
        minimumGamma = Par_new[0]
        minimumA = Par_new[1]
        minimumBeta = Par_new[2]
      else:
        minimumE = minimumE
        minimumVar2 = minimumVar2
        minimumGamma = minimumGamma
        minimumA = minimumA
        minimumBeta = minimumBeta
      Par = Par_new.copy()
      list_gamma.append(Par[0])
      list_a.append(Par[1])
      list_beta.append(Par[2])
      list_EvalE.append(E_mean)
      list_E.append(E_mean)
      list_variance2.append(Var2)
      print("ACCEPTED: \nE_mean =", E_mean, "\nWith parameters:","\nGamma = ",list_gamma[-1],"\na = ",list_a[-1],"\nbeta = ",list_beta[-1],"\nand variance^2: ",Var2)
    else:
      Par = Par.copy()
      list_EvalE.append(list_EvalE[l-1])
      
      

# RESULTS OUTPUT ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------      
      
print("last values:", "\nGamma = ",list_gamma[-1],"\na = ",list_a[-1],"\nbeta = ",list_beta[-1], "\nEnergy = ", list_E[-1],"\nand variance: ",sqrt(abs(list_variance2[-1])))
print("minimum values:", "\nGamma = ",minimumGamma,"\na = ",minimumA,"\nbeta = ",minimumBeta, "\nEnergy = ", minimumE,"\nand variance: ",sqrt(abs(minimumVar2)))

outfile.write("       gamma                         a                        beta                        energy                   variance\n")
for i in range(len(list_E)):
    outfile.write(f"{list_gamma[i]}        {list_a[i]}        {list_beta[i]}            {list_E[i]}         {sqrt(abs(list_variance2[i]))}\n")
    
outfile.close()


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
cmap = plt.cm.cool
img = ax.scatter(list_gamma, list_a, list_beta, s = 10, c=list_E, vmin = -27, vmax = 0, cmap = cmap)
cbar = fig.colorbar(img, label= "Energy (MeV)", shrink = 0.5, orientation = 'horizontal')
ax.set_xlabel('gamma')
ax.set_ylabel('a')
ax.set_zlabel('beta')
plt.title('Energy Minimum', fontsize = 30)
save_fig("VMC4HE")


