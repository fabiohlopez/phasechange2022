# Solving phase change heating from couple of 
#ordinary differential equations
import numpy as np
import math

#N  = int(input("Enter the number of nodes: ")),
#kk = float(input("Enter the number of time steps: "))

#******** Iterative******#
N = 108
kk = 1 # Iterative steps
dt = 4. #step in seg

#####Fluid Properties#######
mf = 0.01 # Kg/seg mass flow fluid
cpf = 1000. # J/Kg K specific heat fluid
Cpf = 1000. # J/K specific heat fluid

#####Solid Properties (Details for Pl16 Wax Paraffin)#######
######from Characterization of Alkanes and Paraffin Waxes for Application as Phase Change Energy Storage Medium SYUKRI HIMRAN  ARYADI SUWONO #####
m = 1. # Kg amount of substance
rho = 800. #Kg/m3 density of substance
lamda = 1.#266 kJ/kg latent heat, release during phase change of substance
k = 0.24 # W/m K thermal conductivity
cs = 2750. # J/kg K specific heat substance
Cs = 1. # J/K specific heat substance
Tm = 308. # K Melting temperature
Tml = 303. # K Liquid melting temperature
Tms = 313. # K Solid melting temperature

#####Storage Recipied#######
At = 2. #m2 cross area of storage
dA = At/N #m2 element area of storage
#Ta = 298. # average ambient temperature in K
length = 2. #length of the storage in meters
U = 1. # losses W/m2 K obtained from UA product
Tref = 1. # Reference temperature in K
#Ta = 298. + 6*sin(j*8640*3.1415/86400.) # ambient temperature in K
P = 3.1415 # wetted perimeter in meters
X = 0. # liquid fraction

theta = (mf * Cpf)/(rho * At * length) # J/Celcius
NTU = U * At * length / (mf * cpf) # Sin Unidades
dtheta = dt * (mf * Cpf)/(rho * At * length) # J/Celcius
omega = (1 - np.exp((-1) * NTU/N)) # Sin Unidades


#float Ts[N], Tss[N], Tf[N],u[N], ai[N], bi[N]

u = np.zeros(N+1)           # array of u[n] values
Ta = np.zeros(N+1)
Ts = np.zeros(N+1)
Tss = np.zeros(N+1)
Tf = np.zeros(N+1)
qac = np.zeros(N+1)

a = np.zeros(N+1)
b = np.zeros(N+1)


for i in range(0, N , 1):
    b [i] = dtheta *0.5 * (rho * At * length * N * omega + U * At /(mf*Cpf))
    Tf [i] = Tss [i] = Ts [i] =  Ta [i] = 298.15
    
with open("Tamb2.txt", "r") as archivo:
    #Tamb = archivo.read()
    Tamb = archivo.readlines()
    Ta = list(Tamb)
Ta = [float(item.rstrip()) for item in Ta]



import matplotlib.pyplot as plt
x = np.linspace(0, length, N+1)

while(kk<=109):
    kk = kk + 8
    Tf [0] = 291.149 - 0.0127* kk * kk + 1.419 * kk
    for i in range(0, N , 1):        
        a [i] = dtheta * (rho * At * length * N * omega * Tf [i] + U * At * Ta [i]/(mf*Cpf)) 
        Tss [i] = Ts [i]
        XX = X
        Ts [i] = ( (1 - b [i]) * Tss [i] + a [i] + (lamda * (X - XX) / Cs) ) / (1 + b [i])
        Ts [N] = Ts [N-1]
        
        Tf [i+1] = (1. - omega) * Tf [i] + omega * 0.5 * (Ts [i] + Tss [i])
        
        u [i] = Cs * (Ts [i] - Tref) + X * lamda
        
        if Ts [i] < Tml:
            X = 0.
        elif (Ts [i] >= Tml) & (Ts [i] <= Tms) :
            X = (Ts [i] - Tms) / (Tml - Tms)
        else:
            X = 1.
        
        #qac [i] = At * rho *cs * (1 - eps) * (Ts [i] - Tss [i])
    
    # Plot Results    
    # Note that even in the OO-style, we use `.pyplot.figure` to create the figure.
    fig, ax = plt.subplots()  # Create a figure and an axes.
    ax.plot(x, Ts, label='Solid Temp')  # Plot some data on the axes.
    #ax.plot(x, u, label='energy')  # Plot more data on the axes...
    ax.plot(x, Tf, label='Fluid Temp')  # ... and some more.
    ax.set_xlabel('axial direction')  # Add an x-label to the axes.
#ax.set_ylabel('y label')  # Add a y-label to the axes.
    ax.set_title("Result")  # Add a title to the axes.
    ax.set_ylim(283., 349.)
    ax.legend()  # Add a legend.
    plt.show()
        #print  "%.2f %.2f" % ( Tf [i], Ts [i])
    #print ("\n")


def plot_history(history, labels):
  """Plots a simulation history."""
  history = np.array(history)
  t = history[:,0]
  n = len(labels) - 1
  plt.figure(figsize=(8,1.95*n))
  
  plt.show ()
