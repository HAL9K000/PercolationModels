# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 09:23:19 2020

@author: Koustav
"""


import numpy as np
import os
import math
import matplotlib.pyplot as plt
import seaborn as sea
import pandas as pan
from scipy.optimize import curve_fit
import matplotlib.ticker as mtick


def plotter():
    
    # Plots all the necessary plots.
    
    p=0.59274; G=[32]; 
    
    os.chdir("../")
    
    os.chdir("figures")
    if(os.path.isdir("Junkyard_Of_The_Unwanted_Frames")==False):
        os.mkdir("Junkyard_Of_The_Unwanted_Frames")
    os.chdir("Junkyard_Of_The_Unwanted_Frames")
    if(os.path.isdir("p_%6.5f" %(p))==False):
        os.mkdir("p_%6.5f" %(p))
    os.chdir("../../") 
    #Back in the master directory
    
    for i in range(0,60,10):
        for g in G:
            usecol= list(range(0,g))
            data= np.genfromtxt('simulations/dump/SP_p_%6.5f_G_%d_i_%d.csv' %(p, g, i), delimiter=",", comments='#', usecols= usecol)
            print(data.shape)
            
            
            #Plotting starts.
            
            zappa=0; non_spn_clus=0 ; num_non_spn=0;
            
            for j in range(0,g):
                for k in range(0,g):
                    
                    if( data[j,k] != 0):
                        # Site is occupied.
                        
                        if( data[j,k] > 0):
                            #Cluster, but not spanning cluster present.
                            plt.plot(k,g-1-j, marker='o', markerfacecolor='none', markeredgecolor='0.85')
                            plt.plot(k,g-1-j, marker='$%d$' %(data[j,k]), color= 'c', markersize=5)
                            
                            non_spn_clus += 1
                            if ( data[j,k] > num_non_spn):
                                num_non_spn = data[j,k]
                            
                        elif (data[j,k] == -1):
                            #Spanning cluster.
                            zappa+=1
                            plt.plot(k, g-1-j, 'ko')
                            
            print("For frame value %d at grid size %d:" %(i, g))
            th_perc_strngth = float(zappa)/(g*g)
            avg_clus_siz =0.0
            if (zappa == 0):
                avg_clus_siz = float(non_spn_clus)/(num_non_spn)
            else:
                avg_clus_siz = float(non_spn_clus)/(num_non_spn - 1)
            print("Reconstructed P[p]:  %f \t &  S[p]:  %f" %(th_perc_strngth, avg_clus_siz))
            
            
            plt.title("Reconstructed P[p]:  %f   &  S[p]:  %f" %(th_perc_strngth, avg_clus_siz))
            plt.savefig("figures/Junkyard_Of_The_Unwanted_Frames/p_%6.5f/L_%d_Frame_%03d.png" %(p, g, i), dpi=300)
            plt.clf()
            plt.cla()
            plt.close()
            
plotter()