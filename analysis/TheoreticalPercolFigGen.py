# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 16:08:45 2020

@author: Koustav
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sea
import pandas as pan
import os


def plotter():
    
    os.chdir("../")
    # Changing to root directory
    census =12
    r_init =12
    lag=12
    p_start=0.65
    p_end=0.78
    L= [25]
    
    masterbinder={}
    #Empty dictionary to be used as Pandas frame.
    
    masterbinder['p (Birth Prob)']=[]; masterbinder['P (Theoretical Perc Value)']=[]; masterbinder['L (Grid Length)']=[]
    
    for g in L:
        #Iterating over grid sizes.
        data=np.genfromtxt('simulations/Theoretical_Percol/DP_L_%d_p1_%3.2f_p2_%3.2f_Cen_%d_R_%d_Lag_%d.csv' %(g, p_start, p_end, census, r_init, lag), delimiter=",", comments='#', skip_header=1)
        
        for x in range(0, len(data[:,0])):
            #Iterating over all values in a column.
            masterbinder['p (Birth Prob)'].append(data[x,0])
            masterbinder['P (Theoretical Perc Value)'].append(data[x,1])
            masterbinder['L (Grid Length)'].append(g)
            
    heartbreaker= pan.DataFrame(masterbinder)
    
    heartbreaker["L (Grid Length)"] = ["$%s$" % x for x in heartbreaker["L (Grid Length)"]]
    #To fix hue legends
    
    os.chdir("figures") #Changing directory to stores images
    
    g= sea.lineplot(x="p (Birth Prob)", y="P (Theoretical Perc Value)",hue="L (Grid Length)", estimator=None, ci=None ,data=heartbreaker)
    plt.ylim(0,1)
    plt.xlim(0.5, 0.9)
    plt.axvline(0.70548, 0, 1, color='0.65', ls='--') #Theoretical Percolation Point For Infinite 2D Grid.
    plt.axvline(0.728, 0, 1, color='0.9', ls='--') #Ayan's Estimated Percolation Point
    
    census1 = 12    
    plt.savefig("DP Theoretical Percolation Results Mult Ran Trials (Census %d).png" %(census1), dpi=300)
    print("Padua")
    
    os.chdir("../analysis")
    #Returning to original directory we started out in.
    
    
def plotter_acf():
    os.chdir("../")
    # Changing to root directory
    print("Enter details of the CSV file to be read below.")
    g = int(input("Enter Grid Size:\t"))    
    p = float(input("Enter Birth Probability:\t"))
    div= int(input("Enter number of divisions (trials):\t"))
    length= int(input("Enter length:\t"))
    
    masterbinder={}
    #Empty dictionary to be used as Pandas frame.
    
    masterbinder['ΔT (Time Difference)']=[]; masterbinder['ACF(ΔT)']=[]
    
    
    data=np.genfromtxt('simulations/ACF/NP_L_%d_p_%3.2f_Div_%d_Len_%d.csv' %(g, p, div, length), delimiter=",", comments='#', skip_header=2)
    for x in range(0, len(data[:,0])):
        #Iterating over all values in a column.                   
        masterbinder['ΔT (Time Difference)'].append(data[x,0])
        masterbinder['ACF(ΔT)'].append(data[x,1])
        
    heartbreaker= pan.DataFrame(masterbinder)
    
    os.chdir("figures") #Changing directory to stores images
    
    xlength= int(input("Enter length of x-axis (cut-off point):\t"))
    ylength_low =float(input("Enter lower bound on plotted y-axis:\t"))
    ylength_hi =float(input("Enter higher bound on plotted y-axis (select value slightly above 1):\t"))
    
    ge= sea.lineplot(x="ΔT (Time Difference)", y="ACF(ΔT)", estimator='mean', ci='sd' ,data=heartbreaker)
    plt.xlim(0, xlength)
    plt.ylim(ylength_low, ylength_hi)
       
    plt.savefig("ACF NP (L--%d n--%d p--%3.2f).png" %(g, div, p), dpi=300)
    print("Padua")
    
    os.chdir("../analysis")
    #Returning to original directory we started out in.
        
    
        
        


plotter() 
    
        
            
        