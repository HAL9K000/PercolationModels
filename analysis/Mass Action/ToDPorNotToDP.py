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
    crt_exp_nautilus()
    # Refer to inline descriptions of these functions to discern the layout and organisation of the CSV data.
    
    base_path = r"DP\12th Jan\Output"
    files = glob.glob(base_path + "/**/*.csv", recursive=True)
    
def crt_exp_nautilus():
    '''This function will extract data from CSVs ( each corresponding to a specific p-point).
       After plotting data from that CSV, it will move onto other CSVs. Data stored in CSVs as:
       
       |  p  ;  #  ;  t  ;  P(t)  ;  N(t)  |
       
       '''
       
    b= int(input("Enter a value you would like to use for delta(t) and theta(t) determination [25-100 may be a half-decent place to start off]:\t"))
       
    g=256 #Occupation Probability value at which simulations were run.
    counter =0
    delt_thet_t = [] #2D array to store values of delta(t) and theta(t) for all encountered values of p.
    
    T=[b, 2*b, 4*b, 6*b, 8*b, 10*b, 12*b, 14*b, 16*b, 18*b, 20*b, 30*b, 40*b]
    #Creating Header Line for delt_thet_t
    header="# p,"
    for t in T:
        header+= "del(%d)," %(t)
    for t in T:
        header+= "th(%d)," %(t)
    header+= "<delta>, <theta>, sigma(delta), sigma(theta)"
    
    base_path = r"DP\18thJan\0.62"
    files = glob.glob(base_path + "/*.csv", recursive=True)
    for file in files:
        if (os.path.getsize(file) > 1024):
            #Filtering for file sizes that are greater than 1024 Bytes in size.
            print(file)
            
            counter +=1
            data_temp = np.genfromtxt('%s' %(file), delimiter=",", comments='#', skip_header=1)   
            
            plt_nautilus_data(data_temp, g, counter, delt_thet_t, T, b)
            
            del data_temp #Flushing Data Temp from Memory.
    
    delt_thet_t = np.array(delt_thet_t)
    np.savetxt("%s_G_%d.csv" %(base_path, g), delt_thet_t, delimiter=",", header=header, comments='# ')
    
    min_d =[0, 0, 0, 0, 0] #Stores indices of min arguments of sig_del
    min_t =[0, 0, 0, 0, 0] #Stores indices of min arguments of sig_th
    
    for i in range(0, delt_thet_t[:,-1].size):
        if( delt_thet_t[i,-2] < delt_thet_t[min_d[-1],-2]):
            for x in range(0,5):
                if(delt_thet_t[i,-2] < delt_thet_t[min_d[x],-2]):
                   for y in range(4, x, -1):
                      min_d[y]= min_d[y-1]
                   min_d[x] = i
                   break
        if( delt_thet_t[i,-1] < delt_thet_t[min_t[-1],-1]):
            for x in range(0,5):
                if(delt_thet_t[i,-1] < delt_thet_t[min_t[x],-1]):
                   for y in range(4, x, -1):
                      min_t[y]= min_t[y-1]
                   min_t[x] = i
                   break
    print("5 most accurate values of delta:")
    for i in range(0,5):
        print("p:\t %f \t <del>:\t %f sig_del:\t %f" %(delt_thet_t[min_d[i],0],delt_thet_t[min_d[i],-4], delt_thet_t[min_d[i],-2]))
    print("5 most accurate values of theta:")    
    for i in range(0,5):
        print("p:\t %f \t <theta>:\t %f sig_th:\t %f" %(delt_thet_t[min_t[i],0],delt_thet_t[min_t[i],-3], delt_thet_t[min_t[i],-1]))
            
def plt_nautilus_data(data, g, count, delt_thet_t, T, b):
    
    '''Remember:
       |  p  ;  #  ;  t  ;  P(t)  ;  N(t)  |
       
       '''
    
    p_c = 0.728
    p = data[0,0] #Occ Prob at Which Simulations were Run
    
    print("%d \t count---- %d \t p1---- %6.5f" %(g,count, p))
    
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
    if(os.path.isdir("%4.3f" %(p))==False):
        os.mkdir("%4.3f" %(p))
    os.chdir("%4.3f" %(p))
    
    
    #We are creating a dictionary to sort (and average) the raw data by time stamp.
    
    half_step_Pt={}; half_step_Nt={}
    
    for i in range(0, 5000):
        half_step_Pt[i] =[]; half_step_Nt[i] =[]
        #Each dictionary key will store all the data values associated with that time stamp.
    
    for x in range(0, data[:,0].size):
        half_step_Pt[data[x,2]].append(data[x,3]) #Generating list of all P(t) values across all random trials for a given t.
        half_step_Nt[data[x,2]].append(data[x,4]) #Generating list of all N(t) values across all random trials for a given t.
     
    CrtData=[] #Array that stores processed data.
    '''Processed data has form: |  t  ;  <P(t)>  ;  <N(t)>  ;  sigma_P(t)  ;  sigma_N(t)  | '''
    for key in half_step_Pt:
        #Iterating through all the time steps.
        Pt_mean = np.mean(half_step_Pt[key]); Nt_mean = np.mean(half_step_Nt[key])
        Pt_sd = np.std(half_step_Pt[key]); Nt_sd = np.std(half_step_Nt[key])
        CrtData.append([key, Pt_mean, Nt_mean, Pt_sd, Nt_sd])
    
    CrtData=np.array(CrtData) #Converting it into Numpy Format
    
    x1= np.transpose(CrtData[3:4915,0]) #Taking a slice corresponding to e^1 to e^8
    x2 =np.transpose(CrtData[3:4915,1])
    
    hurtlocker= pan.DataFrame(data, columns= [ "p",  "Trial #", r"t", r"P(t)", r"N(t)"])
    
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([1, -0.413]))
    #Taking delta= -0.413 as a rough guess
    
    perr = np.sqrt(np.diag(pcov)); tukan= (popt[0], -popt[1], perr[1]) #Storing coefficients of delta
            
    fm= sea.lineplot(data=hurtlocker, x= r"t" , y=r"P(t)", estimator='mean', ci='sd')
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ P(t) = %5.4f \times t^{-(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    fm.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    fm.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    plt.legend()
    fm.set_title(r'$p = %f$' %(p))
    plt.savefig("P(t) vs t (g--%d p--%f).png" %(g,p), dpi=400)
    if(count%33 ==1):
        plt.show()
    plt.close()
    
    x2 =np.transpose(CrtData[3:4915,2])
    popt1, pcov1 = curve_fit(pow_law, x1, x2) #p0= np.asarray([1, 0.16])
    #Taking delta= -0.413 as a rough guess
    
    perr1 = np.sqrt(np.diag(pcov1)); tukan= (popt1[0], popt1[1], perr1[1]) #Storing coefficients of delta
            
    fm= sea.lineplot(data=hurtlocker, x= r"t" , y=r"N(t)", estimator='mean', ci='sd')
    plt.plot(x1, pow_law(x1, *popt1), 'm--', label=r'Th Fit: $ N(t) = %5.4f \times t^{(%5.4f \pm %5.4f)} $ ' % tukan )
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    fm.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    fm.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    plt.legend()
    fm.set_title(r'$p = %f$' %(p))
    plt.savefig("N(t) vs t (g--%d p--%f).png" %(g,p), dpi=400)
    if(count%33 ==1):
        plt.show()
    plt.close()
    
    L=[] #Will store row that will be appended to delt_thet_t
    
    L = delt_thet_determination(p, CrtData, T, L, b)
    #Stores delta and theta(t) values for given p and b.
    
    L.extend([-popt[1], popt1[1], perr[1], perr1[1] ]) 
    #Appending [<del>, <theta>, sig_del, sig_theta]
    
    delt_thet_t.append(L) #Appending completed row.
    
    del L; del CrtData; del x1; del x2; del half_step_Pt; del half_step_Nt
    
            
    os.chdir("../../../../../analysis/Mass Action")
    #Reverting to root directory
        
def delt_thet_determination(p, CrtData, T, L, b):
    # Homes in a new row to delt_thet_t corresponding to a particular p.

    print(T)
    L=[p]
    for t in T:
        t_b = int(t/b)
        if(CrtData[t,1] == 0):
            #If <P(t)> == 0 for some t.
            delt_t= -1
        else:
            delt_t= -math.log10((CrtData[t,1]/CrtData[t_b,1]))/math.log10(b)
        L.append(delt_t)
    for t in T:
        #If <N(t)> == 0 for some t.
        t_b = int(t/b)
        if(CrtData[t,2] == 0):
            thet_t= -1
        else:
            thet_t= math.log10((CrtData[t,2]/CrtData[t_b,2]))/math.log10(b)
        L.append(thet_t)
    
    return L
        
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
        base_path = r"DP\12thJan\0.63"# + "\\" + str(i)
        #print(base_path)
        
        files = glob.glob(base_path + "/*.csv", recursive=True)
        for file in files:
            if (os.path.getsize(file) > 1024):
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
    
    m_data = unspliced_data
    
    plt_nuu(m_data)
    
    
def plt_nuu(nudata):
    # Plots average cluster size data using Seaborn and stores it to the relevant directory.
    '''
    
    |  p  ;  #  ;  t  ;  P(t)  ;  N(t)  |
        
    '''
    
    p_c = 0.728
    
    g = 256 #Occupation Probability value at which simulations were run.
    p1 = int(nudata[0,0]) #Starting Grid Size
    p2 = int(nudata[-1,0])    #Ending Grid Size.
    
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
        
        
    