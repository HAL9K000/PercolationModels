# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 09:14:28 2021

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


'''Same as Mass Specter.py but for DP'''



def ticks(y, pos):
    return r'$e^{:.0f}$'.format(np.log(y))


def ln_pow_law(x, a, tau):
    return a - tau*x

def pow_law(x, a, expo):
    return a*(np.power(x, expo))


def starter_pack():
    
    
    #crt_exp_gamma_beta()
    crt_exp_nuu()
    # Refer to inline descriptions of these functions to discern the layout and organisation of the CSV data.
    
    base_path = r"DP\12th Jan\Output"
    files = glob.glob(base_path + "/**/*.csv", recursive=True)
        
def crt_exp_nuu():
    ''' The CSV data is grouped into the following columns:
        The first column contains the occupation probability ( p ----> p_c = 0.728 for DP class) at which the simulations were run.
        The second column contains the grid size of the lattices, while the third stores the current simulation number.
        The fourth and fifth column usually store the p and the p^2 values for a given simulation, at which the given system
        has been found to percolate for the first time (Ahorny, Stauffer, Dietrich Pgs 70-75)
        
        |  p  ;  #  ;  t  ;  P(t)  ;  N(t)  |
        
        '''
        
    unspliced_data = []
    #This will store the raw unspliced data without any processing from all the CSV files.

    for i in range(0,19):
        base_path = r"DP\12th Jan\Output" + "\\" + str(i)
        #print(base_path)
        
        files = glob.glob(base_path + "/**/*.csv", recursive=True)
        for file in files:
            if (os.path.getsize(file) > 512):
                #Filtering for file sizes that are greater than 512 Bytes in size.
                print(file)
                
                data_temp = np.genfromtxt('%s' %(file), delimiter=",", comments='#', skip_header=1)
                
                if(len(unspliced_data) == 0):
                    #First chip off the old block.
                    unspliced_data = data_temp
                else:
                    if( unspliced_data[-1,0] == data_temp[-1,0]):
                        ''' The current CSV stores experiments performed using the same relevant parameters (specifically grid size)
                        as the last CSV file that was processed.'''
                        data_temp[:,2] += unspliced_data[-1,2]
                        print("Boris")
                        # We update the trial numbers in the new CSV file to be contiguous with the trial numbers from
                        # the previous CSV, if the experimental data is the same.
                    unspliced_data = np.concatenate((unspliced_data, data_temp), axis=0)
            
    #m_data = exxon_split(unspliced_data)
    
    #m_data = unspliced_data
    
    plt_nuu(unspliced_data)
    
    
def plt_nuu(nudata):
    # Plots average cluster size data using Seaborn and stores it to the relevant directory.
    '''
    
    |  p  ;  #  ;  t  ;  P(t)  ;  N(t)  |
        
    '''
    nudata =np.array(nudata)
    
    p_c = 0.728
    
    g = 256 #Occupation Probability value at which simulations were run.
    p1 = nudata[0,0] #Starting Grid Size
    p2 = nudata[-1,0]    #Ending Grid Size.
    
    print("%d \t p1---- %6.5f \t p2---- %6.5f" %(g,p1,p2))
    
    os.chdir(r"..\..\figures")
    # Changing to relevant directory.
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("DP")==False):
        os.mkdir("DP")
    os.chdir("DP")
    if(os.path.isdir("P_C")==False):
        os.mkdir("P_C")
    os.chdir("P_C")
    
    
    
    p_st=0.615; p_end= 0.63;
    a=0; b=0; gingerman=[]
    
    for x in range(0, nudata[:,0].size):
        
        if(nudata[x,0] < p_st):
            continue
        elif(nudata[x,0] >= p_st and p_st ==0.615 and a==0):
            b=-1; a=x; #First time the above condition is met.
        
        elif(b!=0 and nudata[x,0] > p_st):
            b=x
            print("For size %f, we have b = %d  and  b - a = %d" %(g, b, (b-a)))
            
            if(os.path.isdir("%6.5f" %(nudata[x,0]))==False):
                os.mkdir("%6.5f" %(nudata[x,0]))
            os.chdir("%6.5f" %(nudata[x,0]))
            
            npframe= nudata[a:b,:]
            
            hurtlocker= pan.DataFrame(npframe, columns= [ "p",  "Trial #", r"t", r"P(t)", r"N(t)"])
            
            fm= sea.lineplot(data=hurtlocker, x= r"t" , y=r"P(t)", estimator='mean', ci='sd')
            plt.yscale('log', basey= math.e)
            plt.xscale('log', basex= math.e)
            fm.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
            fm.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
            plt.legend()
            fm.set_title(r'$p = %f \quad ( \xi \longrightarrow \infty ) $' %(p_st))
            plt.savefig("P(t) vs t (g--%d p--%8.7f).png" %(g, p_st), dpi=400)
            plt.close()
            
            os.chdir("../")
            
            a=x; p_st = nudata[x,0]
            
        elif(nudata[x,0] > p_end):
            break;
    
            
            
starter_pack() 
        
        
    