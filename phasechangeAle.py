# Solving phase change heating from couple of 
#ordinary differential equations
import numpy as np
import math
import matplotlib.pyplot as plt

#******** Iterative******#
N = 108
kk = 1 # Iterative steps
dt = 4. #step in seg
#dt critico = 13375 o *1/6 = 2229 seg, segun paper

#####Fluid Properties#######
mf = 0.01 # Kg/seg mass flow fluid
cpf = 1000. # J/Kg K specific heat fluid


#####Solid Properties#######
m = 1. # Kg amount of substance
rho = 2320. # Kg/m3 density of substance
k = 2.15 # W/m K thermal conductivity
cs = 810. # J/kg K specific heat substance


#####Storage Recipied#######
At = 1. # m2 cross area of storage
dA = At/N # m2 element area of storage
length = 2. # length of the storage in meters
dx = length/N
U = 1. # losses W/m2 K obtained from UA product
D = 0.02 # stone's diameter
G0 = mf / At
#hv = 650 * pow((G0 / D), 0.7)
hv = 200. # volumetric heat transfer coefficient
eps = 0.42 # porosity


#theta = (mf * cpf)/(rho * cs * (1 - eps) * At * length) # J/Celcius
NTU = hv * At * length / (mf * cpf) # Sin Unidades
dtheta = dt * (mf * cpf)/ (rho * cs * (1 - eps) * At * length) # J/Celcius
omega = 1. - np.exp((-1) * NTU/N) # Sin Unidades



u = np.zeros(N+1)           # array of u[n] values
Ta = np.zeros(N+1)
Ts = np.zeros(N+1)
Tss = np.zeros(N+1)
Tf = np.zeros(N+1)
qac = np.zeros(N+1)

a = np.zeros(N+1)
b = np.zeros(N+1)


for i in range(0, N + 1 , 1):
    b [i] = dtheta *0.5 * (N * omega + U * dA /(mf*cpf))
    Tf [i] = Tss [i] = Ts [i] = Ta [i] = 288.15

with open("Tamb2.txt", "r") as archivo:
    #Tamb = archivo.read()
    Tamb = archivo.readlines()
    Ta = list(Tamb)
Ta = [float(item.rstrip()) for item in Ta]

#with open("Tcolector.txt", "r") as archivo:
    #Tamb = archivo.readlines()
    #Tf = list(Tamb)
#Tf = [float(item.rstrip()) for item in Tf]



x = np.linspace(0, length, N+1)

while(kk<=109):
    kk = kk + 8
    Tf [0] = 291.149 - 0.0127* kk * kk + 1.419 * kk  
    for i in range(0, N , 1):        
        a [i] = dtheta * ( N * omega * Tf [i]  + U * dA * Ta [int(kk)] * (1. /mf*cpf) )
        Tss [i] = Ts [i]
        Ts [i] =  a [i]/ (1 + b [i]) + ( (1 - b [i]) * Tss [i]) / (1 + b [i])  
        Ts [N] = Ts [N-1]
        Tf [i+1] = (1.0 - omega) * Tf [i] * 0.99 + omega * 0.5 * (Ts [i] + Tss [i])
        
        #qac [i] = At * rho *cs * (1 - eps) * (Ts [i] - Tss [i])
        
    # Plot Results    
    fig, ax = plt.subplots()  # Create a figure and an axes.
    ax.plot(x, Ts, label='Solid Temp')  # Plot some data on the axes.
    ax.plot(x, Tf, label='Fluid Temp')  # ... and some more.
    
    ax.set_xlabel('axial direction')  # Add an x-label to the axes.
    ax.set_ylabel('temperature')  # Add a y-label to the axes.
    ax.set_title("Result")  # Add a title to the axes.
    ax.set_ylim(283., 369.)
    ax.legend()  # Add a legend.
    ax.grid()
    plt.show()
    

