# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 11:20:33 2021

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


def yandu():
    
    for i in range(0,7):
        base_path = r"007Apres" + "\\" + str(i)
        
        log = open(base_path + r"\log.txt", "r")
        log_list = log.readlines()
        g = int(log_list[0]); p = float(log_list[1]); c = int(log_list[2])
        print("Grid Size: %d\t p: %f\t Length: %d\t" %(g,p,c))
        
        unspliced_data = []
        #This will store the raw unspliced data without any processing from all the CSV files.
        
        r=0
        files = glob.glob(base_path + "/*.csv")
        for file in files:
            if (os.path.getsize(file) > 512):
                r+=8
                #Filtering for file sizes that are greater than 512 Bytes in size.
                print(file)
            
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
            
                '''
                data_temp resembles:
                    | t, Correl, k, l |
                '''
                
                if(len(unspliced_data) == 0):
                    #First chip off the old block.
                    unspliced_data = data_temp
                else:
                    unspliced_data = np.concatenate((unspliced_data, data_temp), axis=0)
                    #Concatanated.
        a,b = unspliced_data.shape
        print("Unspliced Array Size:\t (%d, %d)" %(a,b))
        #ko=input("Enter to continue:\t")
        
        yandex(unspliced_data, g, p, c, r)
        

        
def yandex(unspliced_data, g, p, c, r):
    
        #Plot the shit out of Compton.
        os.chdir("../../../figures")
            
        if(os.path.isdir("CrossCor")==False):
            os.mkdir("CrossCor")
        os.chdir("CrossCor")
        if(os.path.isdir("DP")==False):
            os.mkdir("DP")
        os.chdir("DP")
        
        f=0; hurtlocker=0;
        
        if(g == 256):

            unspliced_data[:,0]/= (g*g)
            unspliced_data = exxon_split(unspliced_data, c)
            print('First few lines of final matrix:')
            L=[]
            for x in range(2500,3000):
                if(x%r == 0):
                    print("\n")
                    if(len(L) > 0):
                        print("SD of Max Corr for t = %f is:\t %f" %( unspliced_data[x,0], np.std(L) ) )
                        L=[]
                print("%f \t %6.5f \t  %f \t %f" %(tuple(unspliced_data[x,:])))
                L.append(unspliced_data[x,1])
                
            hurtlocker= pan.DataFrame(unspliced_data, columns= [r"t [$N^2$ Updates]", "Max Cross-Correlation",  "i", "j"])

            x1 =np.transpose(unspliced_data[:,0])
            x2= np.transpose(unspliced_data[:,1])
    
            f= sea.lineplot(data=hurtlocker, x=r"t [$N^2$ Updates]" , y="Max Cross-Correlation", estimator= np.nanmean, ci='sd', marker="s") #err_style="band")
    
            '''popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.3, -0.005]))
    
            perr = np.sqrt(np.diag(pcov))
        
            print("SD of C:\t" +str(perr[1]) + " for p:\t" +str(p) + " at q:\t" +str(q))
    
            tukan= (popt[0], -popt[1], perr[1])
            plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ P[p,q] = %5.4f \times L^{-(%5.4f \mp %5.4f)} $ ' % tukan )
            plt.legend()'''
            #plt.xlim(g1,g2+20)
            #plt.yscale('log', basey= math.e)
            #plt.xscale('log', basex= math.e)
            #g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
            #g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    
            f.set_title(r'p = %f, Grid Size (G) = %d, n = %d' %(p,g,r))
            plt.savefig("Cross Correlation --- p_%f - Grid Size (G)_%d - n_%d.png" %(p, g, r), dpi=400)
            plt.show()
            plt.close()
            
        else:
            unspliced_data = exxon_split(unspliced_data, c)
            a,b =unspliced_data.shape
            print("Number of rows & columns in unspliced:\t (%d, %d)" %(a,b))
            L=[]
            for x in range(25000,25500):
                if(x%r == 0):
                    print("\n")
                    if(len(L) > 0):
                        print("SD of Max Corr for t = %f is:\t %f \n" %( unspliced_data[x-1,0], np.std(L) ) )
                        L=[]
                print("%f \t %6.5f \t  %f \t %f" %(tuple(unspliced_data[x,:])))
                L.append(unspliced_data[x,1])
            
            hurtlocker= pan.DataFrame(unspliced_data, columns= [r"No. Of Single Updates", "Max Cross-Correlation",  "i", "j"])
            
            x1 =np.transpose(unspliced_data[:,0])
            x2= np.transpose(unspliced_data[:,1])
    
            f= sea.lineplot(data=hurtlocker, x=r"No. Of Single Updates" , y="Max Cross-Correlation", estimator='mean', ci='sd', err_style="band")
            plt.axvline(x= g*g, color='0.65')
            plt.text(g*(g+1),0.8,r'$N^2$ update',rotation=90, color ='0.65')
            f.set_title(r'p = %f, Grid Size (G) = %d, n = %d' %(p,g,r))
            plt.savefig("Cross Correlation --- p_%f - Grid Size (G)_%d - n_%d.png" %(p, g, r), dpi=400)
            plt.show()
            plt.close()


        
        os.chdir(r"..\..\..\analysis\Mass Action\DP")
        
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
        
                    
            
            
            