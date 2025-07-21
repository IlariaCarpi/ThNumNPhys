import os
import numpy as np
from matplotlib import pyplot as plt
from math import exp, sqrt



# FILES INITIALIZATION ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

# saving figure into designated directory 
def save_fig(fig_id):
    plt.savefig(image_path(fig_id) + ".png", format='png')

# saving data in designated file
outfile = open(data_path("DataHO.dat"),'w')




# INITIALIZE VARIABLES ---------------------------------------------------------------------------------------

NAcc = 0 # accepted moves
NP = 200  # number of moves to find parameters
NX = 10000  # number of particles' moves
step = 1.0  # algorithm step for random walk of particles
step1 = 0.8 # step for parameter random walk

list_x = [] # storing positions
list_alpha = [] # storing alpha values
list_E = [] # storing accepted energies
list_Eall = [] # storing energies
list_variance2 = [] # storing square variances


# parameter variables for VMC

BestPar = 1 # final goal
Par = 0.4 # starting point for parameter
Par_new = 0 # initialize intermediate parameter for random walk

# random seed for code reproducibility
np.random.seed(19599)   



# DEFINE FUNCTIONS ------------------------------------------------------------------------------------------------------------------------------

# wave function definition

def WF(x,alpha):     
  return exp(-0.5*alpha*alpha*x*x)

# local energy definition

def LocalEnergy(x,alpha):
  return 0.5*x*x*(1-alpha**4) + 0.5*alpha*alpha


xOld = step * (np.random.rand() - 0.5) # first step
list_x.append(xOld) # add to coordinates list 




# find the set of steps with Metropolis algorithm --------------------------------------------------------------------------------------------------------

for k in range(0,NX):
  xNew = xOld + step*(np.random.rand()-0.5)
  if ((WF(xNew,BestPar)/WF(xOld,BestPar))**2>1 or (WF(xNew,BestPar)/WF(xOld,BestPar))**2 > np.random.rand() ):
    xOld = xNew
    list_x.append(xOld)
    NAcc = NAcc +1
  else:
    xOld = xOld
    list_x.append(xOld) 

print("Acceptance ratio =", 100*NAcc/NX, "%")




# VMC TO FIND BEST PARAMETER ------------------------------------------------------------------------------------------------------------------------------

for l in range(0,NP):

  # random walk
  
  Par_new = Par + (np.random.rand()-0.5)*step1
  while (Par_new<0 or Par_new>1.5):
      Par_new = Par + (np.random.rand()-0.5)*step1

  E_loc = 0 
  E_loc_2 = 0
  D = 0

    
  # calculate local energy for current set
  # technique of correlated sampling
  
  for j in range(0,NX):
    E_loc = E_loc + LocalEnergy(list_x[j],Par_new)*(WF(list_x[j],Par_new)/WF(list_x[j],BestPar))**2
    E_loc_2 = E_loc_2 +(LocalEnergy(list_x[j],Par_new)*(WF(list_x[j],Par_new)/WF(list_x[j],BestPar))**2)**2
    D = D + (WF(list_x[j], Par_new)/WF(list_x[j], BestPar))**2
    
  E_mean = E_loc/D
  Var2 = (E_loc_2/D - (E_loc/D)**2)/D


  # METROPOLIS QUESTION

  if (l==0):
    Par = Par_new  
    list_alpha.append(Par)
    list_E.append(E_mean)
    list_Eall.append(E_mean)
    list_variance2.append(Var2)
    minimumE = E_mean
    minimumAlpha = Par
    minimumVariance2 = Var2
  
  else:   
    print("l = ", l) 
    
    if (list_Eall[l-1]>E_mean or list_Eall[l-1]/E_mean>np.random.rand()):     
    
      if (minimumE>E_mean):
        minimumE = E_mean
        minimumAlpha = Par
        minimumVariance2 = Var2
      else:
        minimumE = minimumE
        minimumAlpha = minimumAlpha
        minimumVariance2 = minimumVariance2

      Par = Par_new
      list_alpha.append(Par)
      list_Eall.append(E_mean)
      list_E.append(E_mean)
      list_variance2.append(Var2)

    else:
      Par = Par
      list_Eall.append(list_Eall[l-1])
      
      


      
# RESULTS OUTPUT ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------      
      
print("\n\n\nlast value:", "\nAlpha = ",list_alpha[-1], "\nEnergy = ", list_E[-1], "\nVariance: ",sqrt(abs(list_variance2[-1])))
print("\n\nENERGY MINIMUM: ",minimumE, "\nALPHA IN MINIMUM: ",minimumAlpha, "\nVARIANCE OF MINIMUM: ",sqrt(abs(minimumVariance2)))


outfile.write("       alpha                    energy                    variance\n")
for i in range(len(list_alpha)):
    outfile.write(f"{list_alpha[i]}        {list_E[i]}        {sqrt(abs(list_variance2[i]))}\n")
    
outfile.close()

plt.subplot(1, 2, 1)
plt.plot(list_alpha, list_E, '.-')
plt.ylabel('Dimensionless energy')
plt.xlabel(r'$\alpha$', fontsize=15)
plt.subplot(1, 2, 2)
plt.plot(list_alpha, list_variance2, '.-')
plt.xlabel(r'$\alpha$', fontsize=15)
plt.ylabel('Variance**2')
plt.suptitle('Energy and variance', fontsize=16)
plt.subplots_adjust(wspace=1)
save_fig("VMCHarmonic")
plt.show()
