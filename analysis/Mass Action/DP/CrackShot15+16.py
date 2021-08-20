# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 11:41:36 2021

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

def expo(x, a, b, c):
    return a + b*(np.exp(c*x))

def yandu():
    
    binder=[]
    #Saves the values of the critical exponent
    
    for i in range(0,40):
        base_path = r"15+16" + "\\" + str(i)
        
        ''' log = open(base_path + r"\log.txt", "r")
        log_list = log.readlines() 
        g = int(log_list[0]); p = float(log_list[1]); c = int(log_list[2]) '''
        g =128; c= 12500; p = 0
        #print("Grid Size: %d\t p: %f\t Length: %d\t" %(g,p,c))
        
        unspliced_data = []
        #This will store the raw unspliced data without any processing from all the CSV files.
        
        r=0
        files = glob.glob(base_path + "/*.csv")
        for file in files:
            if (os.path.getsize(file) > 512):
                r+=16
                #Filtering for file sizes that are greater than 512 Bytes in size.
                print(file)
            
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1, usecols = (0,1,2,3) )
            
                '''
                data_temp resembles:
                    | t, Correl, k+l, p |
                '''
                p = data_temp[0,3]
                
                
                if(len(unspliced_data) == 0):
                    #First chip off the old block.
                    unspliced_data = data_temp
                else:
                    unspliced_data = np.concatenate((unspliced_data, data_temp), axis=0)
                    #Concatanated.
        print("Grid Size: %d\t p: %f\t Length: %d\t" %(g,p,c))
        a,b = unspliced_data.shape
        print("Unspliced Array Size:\t (%d, %d)" %(a,b))
        #ko=input("Enter to continue:\t")
        
        
        yandex(unspliced_data, g, p, c, r, binder)
    heado = 'p, A, SD(A), B, SD(B), C, SD(C), lag_0.95, lag_0.9, lag_0.8, lag_0.75, lag_0.7, lag_0.6'
    #Header for the CSV file    
    np.savetxt("PissingAbout15+16.csv", binder, delimiter=',', header=heado, comments='#')
        
    find_stat(binder) # Find final stats
        
def yandex(unspliced_data, g, p, c, r, bind):
    
        #Plot the shit out of Compton.
        os.chdir("../../../figures")
            
        if(os.path.isdir("CrossCor")==False):
            os.mkdir("CrossCor")
        os.chdir("CrossCor")
        if(os.path.isdir("DP")==False):
            os.mkdir("DP")
        os.chdir("DP")
        if(os.path.isdir("GandU")==False):
            os.mkdir("GandU")
        os.chdir("GandU")
        
        f=0; hurtlocker=0;
        
        if(g == 128):

            unspliced_data[:,0]/= 4
            unspliced_data = exxon_split(unspliced_data, c)
            unspliced_data[:,0]*= 4
            print('First few lines of final matrix:')
            L=[]
            for x in range(25000,25500):
                if(x%r == 0):
                    print("\n")
                    if(len(L) > 0):
                        print("SD of Max Corr for t = %f is:\t %f\n" %( unspliced_data[x-1,0], np.std(L) ) )
                        L=[]
                print("%f \t %6.5f \t  %f \t %f" %(tuple(unspliced_data[x,:])))
                L.append(unspliced_data[x,1])
                
            hurtlocker= pan.DataFrame(unspliced_data, columns= [r"# Updates (m)", "Max Cross-Correlation",  "i", "j"])

            x1 =np.transpose(unspliced_data[:,0])
            x2= np.transpose(unspliced_data[:,1])
    
            f= sea.lineplot(data=hurtlocker, x=r"# Updates (m)" , y="Max Cross-Correlation", estimator= 'mean', ci='sd', err_style="band")
            
            plt.axvline(x= g*g, color='0.65')
            plt.text(g*(g+1),0.8, r'$N^2$ updates',rotation=90, color ='0.7')
            plt.axvline(x= 2*g*g, color='0.65')
            plt.text(2*g*(g+0.5),0.75,r'$2 \times N^2$ updates',rotation=90, color ='0.7')
            plt.axvline(x= 3*g*g, color='0.65')
            plt.text(3*g*(g+0.33),0.75,r'$3 \times N^2$ updates',rotation=90, color ='0.7')
            
    
            popt, pcov = curve_fit(expo, x1, x2, p0= np.asarray([0, 1, -0.05]))
    
            perr = np.sqrt(np.diag(pcov))
        
            print("SD of exponent:\t" +str(perr[2]) + " for p:\t" +str(p))
    
            tukan= (popt[0], popt[1], popt[2], perr[2])
            plt.plot(x1, expo(x1, *popt), 'm--', label=r'Fit: $ R_{0,0}[p,t] = %3.2f + %3.2f \times e^{(%7.6f \mp %5.5f) \times m} $ ' % tukan )
            plt.legend()
            #plt.xlim(g1,g2+20)
            #plt.yscale('log', basey= math.e)
            #plt.xscale('log', basex= math.e)
            #g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
            #g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    
            f.set_title(r'p = %f, Grid Size (G) = %d, n = %d' %(p,g,r))
            plt.savefig("Cross Correlation --- p_%f - Grid Size (G)_%d - n_%d.png" %(p, g, r), dpi=400)
            #plt.show()
            plt.close()
            
            #Time to get the lag timesteps (normalised by N^2).
            
            lag=[]
            
            iaggo= [0.95, 0.9, 0.8, 0.75, 0.7, 0.6]
            
            for i in iaggo:
                lag.append(float((math.log((i - popt[0])/popt[1])/(popt[2]*g*g))))
                
            #Stores the lag.
            
            bind.append([p, popt[0], perr[0], popt[1], perr[1], popt[2], perr[2], lag[0], lag[1], lag[2], lag[3], lag[4], lag[5]])
            
            


        
        os.chdir(r"..\..\..\..\analysis\Mass Action\DP")
        
def find_stat(bind):
    
    os.chdir("../../../figures")
            
    if(os.path.isdir("CrossCor")==False):
        os.mkdir("CrossCor")
    os.chdir("CrossCor")
    if(os.path.isdir("DP")==False):
        os.mkdir("DP")
    os.chdir("DP")
    if(os.path.isdir("GandU")==False):
        os.mkdir("GandU")
    os.chdir("GandU")
    
    bind =np.array(bind)
    blind = bind[:,0:7]
    g= 128; r=16
    hurtlocker= pan.DataFrame(blind, columns= ["p", "A", "SD(A)", "B", "SD(B)", "Decay Rate", "SD(C)"])
    
    f= sea.lineplot(data=hurtlocker, x="p" , y="Decay Rate", estimator= 'mean')
    f.set_title(r'Decay Rate vs p, Grid Size (G) = %d, n = %d' %(g,r))
    plt.savefig("Decay Rate vs p --- Grid Size (G)_%d - n_%d.png" %( g, r), dpi=400)
    plt.show()
    plt.close()
    
    os.chdir(r"..\..\..\..\analysis\Mass Action\DP")
    
    
    
    
        
def exxon_split(unspliced_data, c):
    # Arranges unspliced_data in strict ascending order of time steps
        
    '''But first, we need to have a comprehensive list of all time steps in ascending order)'''
        
    
    T = [i for i in range(0,c)]
            
            
    '''Now for each grid size, all the revelant data from unspliced_data must be spliced out and concatanated into a 
    new array'''
    
    a=0; b=0 #Stores relevant splices for each grid size.
    
    m_splice =[]
    
    flag=0
        
    for t in T:
        #Iterating over all the grid sizes, in unspliced_data
        print(str(t) +"\n")
        for x in range(t,len(unspliced_data[:,2]),c):
            #Iterating over unspliced_data.
            if(t == unspliced_data[x,0]): 
                # We have a new hit for the given grid size "l".
                if (flag == 0):
                    #First one in the bag.
                    m_splice = unspliced_data[x:x+1,:].tolist()
                    print(m_splice); #a=input("Enter to continue")
                    flag = 1; continue;
                
                #m_splice = np.concatenate((m_splice, unspliced_data[x:x+1,:]), axis=0)
                m_splice.extend(unspliced_data[x:x+1,:].tolist())
                
                
    k=0            
    for x in m_splice:
        print(x)
        k+=1
        if(k == 50):
            break
    m_splice = np.array(m_splice)             
    return m_splice;
        
        
yandu()