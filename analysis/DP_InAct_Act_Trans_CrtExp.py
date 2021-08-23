# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 12:06:36 2021

@author: Koustav
"""

import os
import glob
import matplotlib.pyplot as plt
import seaborn as sea
import numpy as np
import pandas as pan
import math
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit

def leppard():
    #Main Street USA.
    
    t=[50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000]
    p_val= [0.6204698, 0.62027059, 0.62007137, 0.62086823]
    #P values for which we will ascertain the delta and sigma local fit data.
    
    raw_data= np.genfromtxt("../simulations/DP_InAct_Act_Trans_CrtExp/0.62_G_256.csv", delimiter=",", comments='#', skip_header=1)
    #Importing data.
    data_delta=np.zeros((len(p_val)*len(t),3)) # Will store | 1/t ; del(t) ; p |
    data_theta=np.zeros((len(p_val)*len(t),3)) # Will store | 1/t ; sig(t) ; p |
    
    count=0
    for x in range(0, raw_data[:,0].size):
        for p in p_val:
            if(raw_data[x,0] == p):
                print(p)
                #There's a hit!
                for k in range(0, len(t)):
                    data_delta[count*len(t)+k,0] = float(1.0/t[k]) #Storing 1/t values
                    data_delta[count*len(t)+k,1] = raw_data[x,k+1] #Storing del(t) values
                    data_delta[count*len(t)+k,2] = p #Storing appropriate p values
                    
                    data_theta[count*len(t)+k,0] = float(1.0/t[k]) #Storing 1/t values
                    data_theta[count*len(t)+k,1] = raw_data[x,k + 1 +len(t)] #Storing del(t) values
                    data_theta[count*len(t)+k,2] = p #Storing appropriate p values
                count+=1 #To augment our constructed data holders.
    
    '''for x in data_theta:
        print(x)
    print("\n\n")
    for y in data_delta:
        print(y)'''
        
    hurtlocker_del= pan.DataFrame(data_delta, columns= [ "1/t",  r"$\delta(t)$", r"$ p \mbox{ (s.t } p \approx p_{c} \mbox{) }$"])
    hurtlocker_the= pan.DataFrame(data_theta, columns= [ "1/t",  r"$\theta(t)$", r"$ p \mbox{ (s.t } p \approx p_{c} \mbox{) }$"])
    
        
    g= sea.lineplot( data=hurtlocker_del, x='1/t', y=r"$\delta(t)$" , hue=r"$ p \mbox{ (s.t } p \approx p_{c} \mbox{) }$", palette="viridis", style= r"$ p \mbox{ (s.t } p \approx p_{c} \mbox{) }$", marker="s")
    plt.legend()
    g.set_title(r"$ \delta_{p}(t) \quad (where \quad p \approx p_{c})$")
    plt.savefig("../figures/CrtExp/DP/P_C/Best Del(t).png", dpi=400)
    plt.show()
    plt.close()
    
    g= sea.lineplot( data=hurtlocker_the, x='1/t', y=r"$\theta(t)$" , hue=r"$ p \mbox{ (s.t } p \approx p_{c} \mbox{) }$", palette="viridis", style= r"$ p \mbox{ (s.t } p \approx p_{c} \mbox{) }$", marker="s")
    plt.legend()
    g.set_title(r"$ \theta_{p}(t) \quad (where \quad p \approx p_{c})$")
    plt.savefig("../figures/CrtExp/DP/P_C/Best Theta(t).png", dpi=400)
    plt.show()
    plt.close()

leppard()