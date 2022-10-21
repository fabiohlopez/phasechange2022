#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy as np
import math
import matplotlib.pyplot as plt
from numpy import linalg as la

''' I switched the WHILE loop for a DO-WHILE loop '''

##### time advancement #####
N = 100;
kk = 1 # Iterative steps
dt = 60. #step in seg


##### Fluid Properties#######
mf = 0.01 # Kg/seg mass flow fluid
cpf = 1007. # J/Kg K specific heat fluid


##### Solid Properties (Details for Pl16 Wax Paraffin)#######
######from Characterization of Alkanes and Paraffin Waxes for Application as Phase Change Energy Storage Medium SYUKRI HIMRAN  ARYADI SUWONO #####
#m = 0.25 # Kg amount of substance
rho = 900. #Kg/m3 density of substance      ###CHECK THIS PARAMETER ###
lamda = 192000. #J/kg latent heat, es 161000 release during phase change of substance
Cs = 2500. # J/kg K average specific heat (between liquid and solid states)
# = 2500. + 28 * 1000. * np.exp( (T[i] - 25.881) / 11.91 ) # J/kg K
x=1./N #liquid fraction


#####Storage Recipied#######
At = 0.0095 #m2 cross area of storage
length = 0.15 #length of the storage in meters
U = 19.5 # losses W/m2 K obtained from UA product
#Tref = 1. # Reference temperature in K
P = 2. * math.pi * math.sqrt (At / math.pi ) # wetted perimeter in meters


#*****Dimensionless parameters*****#
theta = (mf * cpf)/(rho *  At * length * Cs)
dtheta = dt * (mf * cpf)/(rho *  At * length * Cs)
NTU = U * P * length / (mf * cpf)
omega = (1. - np.exp( (-1.) * NTU/N)) 
Y = Cs / dtheta


#*****Vector set ups*****#    
Ta = np.zeros(100*N+1)         #ambient temperature   
T = np.zeros(N+1)          #solid temperature
TT = np.zeros(N+1)         
Tf = np.zeros(N+1)          #fluid temperature
#t = np.linspace(0, S, N+1)  # time mesh

a = np.zeros(N+1)
b = np.zeros(N+1)

#%%
#***** External data*****#
for i in range(0, 100*N +1 , 1): 
    Ta [i] = 293.

##### CI ######
for i in range(0, N+1 , 1): 
    b [i] = 0.5 * ( N * rho * At * length * Cs * omega  + Cs * rho * At * At * length * U  /(mf*cpf)) # J/ºC
    TT [i] = T [i] = 293.
    Tf [i] = 293.

##### BC ######
Tf [0] = Tf [1] = 343. 


###### Print section #####
deltatime = (Y - b[2]) / ( Y + b[2])
dtm = 2. * 0.17 * Cs / (N * omega * mf * cpf + At * U)
formatted_float = "{:.3f}".format(deltatime)
formatted_vector = "{:.3f}".format(dtm)
print(formatted_float)
print (formatted_vector)


#*****Reading ambient temperature file*****#
with open("Tamb2.txt", "r") as archivo:
    #Tamb = archivo.read()
    Tamb = archivo.readlines()
    #Ta = list(Tamb)
#Ta = [float(item.rstrip()) for item in Ta]
#for i in range(0, N , 1):    
    #print  "%.2f" % (Ta [i])
#Ta = np.array(Ta,dtype="float")

#result =  np.la.cond(T)
#print("Condition number of the matrix:")
#print(result)


##### ITERATIVE PROCESS#####
while True:
    xx = x
    x = 0.1 + (kk-5)*(1./N)        
    for i in range(1, N + 1 , 1):     
        if T [i] < 310. or T [i] > 311.:
            a [i] = N * rho * At * length * Cs * omega * Tf [i-1] +  Cs * rho * At * At * length * U * Ta [kk]/(mf*cpf)
        else:
            a [i] = N * rho * At * length * Cs * omega * Tf [i-1] - (lamda/dtheta) * (x - xx) +  Cs * rho * At * At * length * U * Ta [kk]/(mf*cpf)
        TT [i] = T [i]         
        T [i] = (Y - b[2]) * TT [i] / ( Y + b[2]) + a [i] / ( Y + b[2])
        T [0] = T [1] 
        #print ("es el coef positivo?",(Y - b[2]) / ( Y + b[2]))   
        print ("is a variable",a [i])
        if i < N:
            Tf [i] = (1. - omega) * Tf [i-1] + omega * 0.5 * (T [i] + TT [i])
        else:
            Tf [i] = Tf [i-1]
           
       

        # Print section
        if kk % 100 == 0:
            plt.plot(T,'r--o')
            plt.ylabel('PCM temperature')
            #plt.xlabel('Nº of nodes')
            plt.yticks ([300.,320.,340.,360.,380.])
            plt.xticks ([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
            #plt.hold('on')
            plt.savefig('figure'+str(kk)+'.png', format='PNG')
            plt.clf()
            if (la.norm(T)-la.norm(TT))<1.:
                break
    
    kk = kk + 1
    

