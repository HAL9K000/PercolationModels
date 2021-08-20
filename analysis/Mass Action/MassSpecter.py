# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 17:41:47 2020

@author: Koustav
"""
import os
import glob
import math
import matplotlib.pyplot as plt
import seaborn as sea
import numpy as np
import pandas as pan
from scipy.optimize import curve_fit
import matplotlib.ticker as mtick

''' This is a script that is specifically designed to ascertain, read and subsequently operate on data present
    in an arbitrary number of CSV files scattered across various folders.'''
    
    
def ticks(y, pos):
    return r'$e^{:.0f}$'.format(np.log(y))


def ln_pow_law(x, a, tau):
    return a - tau*x

def pow_law(x, a, expo):
    return a*(np.power(x, expo))


def starter_pack():
    
    
    crt_exp_gamma_beta()
    #crt_exp_nuu()
    # Refer to inline descriptions of these functions to discern the layout and organisation of the CSV data.
    
    base_path = r"13th Nov\Output"
    files = glob.glob(base_path + "/**/*.csv", recursive=True)
    m=0; L=[]
    for file in files:
        print(file)
        if (os.path.getsize(file) > 1024):
            print(True)
        else:
            print(False)
        abba = (m, int(''.join(filter(lambda i: i.isdigit(), file)))) 
        #Tuple that stores index number of files list alongside a concatanation of all the digits present inside "file".
        
        L.append(abba)
        m+=1
    
    for x in range(0, len(files)):
        print(L[x], end="\t")
        print(files[x])
        
            
    #sort_files = sorted(files, key=int(''.join(filter(lambda i: i.isdigit(), files))))
    
    
    
def crt_exp_gamma_beta():
    
    ''' The CSV data is grouped into the following columns:
        The first column contains the occupation probability ( p ----> p_c = 0.728 for DP class) at which the simulations were run.
        The second column contains the grid size of the lattices, while the third stores the current simulation number.
        The fourth and fifth column usually store the cluster size and the corresponding n_s(p) value for a given simulation,
        with the exception of the last entry into each of these columns for each trial number, which stores -10 (indicative of a
        spanning cluster) and P[p] respectively.
        '''
        
    unspliced_data = []
    #This will store the raw unspliced data without any processing from all the CSV files.

    for i in range(0,9):
        base_path = r"13th Nov\Output" + "\\" + str(i)
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
                    if( unspliced_data[-1,1] == data_temp[-1,1]):
                        ''' The current CSV stores experiments performed using the same relevant parameters (specifically grid size)
                        as the last CSV file that was processed.'''
                        data_temp[:,2] += unspliced_data[-1,2]
                        print("Boris")
                        # We update the trial numbers in the new CSV file to be contiguous with the trial numbers from
                        # the previous CSV, if the experimental data is the same.
                    unspliced_data = np.concatenate((unspliced_data, data_temp), axis=0)
            
    m_data = exxon_split(unspliced_data)
    
    m_data = np.array(m_data)
    
    #s= input("Enter any key to continue.")
    s=0
    '''for  x in m_data:
        if(s % 25 == 0):
            print("%d \t %d \t %f \t %f" %(int(x[1]), int(x[2]), x[3], x[4]))
        s+=1
    print("Total:\t %d" %(s))'''
    
    ''' Now to thresh out the percolation strength data from m_data.
        Remember, P[p] data is characterised by a -10 entry in the 4th column.'''
        
    split_data = m_data[:,3] == -10
    
    perc_data = m_data[split_data]              #Final data form for percolation strength.
    post_transcript = m_data[~split_data]       #Final data form for average cluster size calculations
    
    '''print("Collated Percolation Data:")
    for x in perc_data:
        print("%d \t %d \t %f \t %f" %(int(x[1]), int(x[2]), x[3], x[4]))
        s+=1
    print("Total:\t %d" %(s))
    s= input("Enter any key to continue.")
    s=0
    for  x in post_transcript:
        if(s % 25 == 0):
            print("%d \t %d \t %f \t %f" %(int(x[1]), int(x[2]), x[3], x[4]))
        s+=1
    print("Total:\t %d" %(s))'''
    
    plt_beta(perc_data)
    #Makes Finite Sized Scaling Plots For The Beta Critical Exponent.
    
    plt_gamma(post_transcript)
    #Makes Finite Sized Scaling Plots For The Gamma Critical Exponent.
    
    
def crt_exp_nuu():
    ''' The CSV data is grouped into the following columns:
        The first column contains the occupation probability ( p ----> p_c = 0.728 for DP class) at which the simulations were run.
        The second column contains the grid size of the lattices, while the third stores the current simulation number.
        The fourth and fifth column usually store the p and the p^2 values for a given simulation, at which the given system
        has been found to percolate for the first time (Ahorny, Stauffer, Dietrich Pgs 70-75)
        
        |  p_c  ;  L  ;  #  ;  p'  ;  (p')^2  |
        
        '''
        
    unspliced_data = []
    #This will store the raw unspliced data without any processing from all the CSV files.

    for i in range(0,19):
        base_path = r"20th Nov\Output" + "\\" + str(i)
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
                    if( unspliced_data[-1,1] == data_temp[-1,1]):
                        ''' The current CSV stores experiments performed using the same relevant parameters (specifically grid size)
                        as the last CSV file that was processed.'''
                        data_temp[:,2] += unspliced_data[-1,2]
                        print("Boris")
                        # We update the trial numbers in the new CSV file to be contiguous with the trial numbers from
                        # the previous CSV, if the experimental data is the same.
                    unspliced_data = np.concatenate((unspliced_data, data_temp), axis=0)
            
    m_data = exxon_split(unspliced_data)
    
    m_data = np.array(m_data)
        
        
    plt_nuu(m_data)
    
    
    
def plt_beta(perc_data):
    # Plots percolation data using Seaborn and stores it to the relevant directory.
    
    p = perc_data[0,0] #Occupation Probability value at which simulations were run.
    g1 = int(perc_data[0,1]) #Starting Grid Size
    g2 = int(perc_data[-1,1])    #Ending Grid Size.
    
    print("%f \t G1---- %d \t G2---- %d" %(p,g1,g2))
    
    os.chdir(r"..\..\figures")
    # Changing to relevant directory.
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("DP")==False):
        os.mkdir("DP")
    os.chdir("DP")
    if(os.path.isdir("Finite Scaling")==False):
        os.mkdir("Finite Scaling")
    os.chdir("Finite Scaling")
    if(os.path.isdir("Beta")==False):
        os.mkdir("Beta")
    os.chdir("Beta")
    
    hurtlocker= pan.DataFrame(perc_data, columns= ["p", "L", "Trial Number",  "-10", r"P[p]"])
    
    x1 =np.transpose(perc_data[:,1])
    x2= np.transpose(perc_data[:,4])
    
    g= sea.lineplot(data=hurtlocker, x="L" , y="P[p]", estimator='mean', ci='sd', marker="s", err_style="band")
    
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.3, -0.005]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of C:\t" +str(perr[1]) + " for p:\t" +str(p))
    
    tukan= (popt[0], -popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ P[p] = %5.4f \times L^{-(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.xlim(math.exp(3),g2+20)
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    plt.legend()
    
    g.set_title(r'$p = %5.4f, ( \xi \longrightarrow \infty ) $' %(p))
    plt.savefig("Log Line Beta P(p) vs L (p--%8.7f) (Range-- %d-%d).png" %(p, g1, g2), dpi=400)
    plt.show()
    plt.close()
    
    
    
    os.chdir(r"..\..\..\..\..\analysis\Mass Action")
    #Returning to our home directory.
    
    
def plt_gamma(post_transcript):
    # Plots average cluster size data using Seaborn and stores it to the relevant directory.
    
    p = post_transcript[0,0] #Occupation Probability value at which simulations were run.
    g1 = int(post_transcript[0,1]) #Starting Grid Size
    g2 = int(post_transcript[-1,1])    #Ending Grid Size.
    
    print("%f \t G1---- %d \t G2---- %d" %(p,g1,g2))
    
    os.chdir(r"..\..\figures")
    # Changing to relevant directory.
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("DP")==False):
        os.mkdir("DP")
    os.chdir("DP")
    if(os.path.isdir("Finite Scaling")==False):
        os.mkdir("Finite Scaling")
    os.chdir("Finite Scaling")
    if(os.path.isdir("Gamma")==False):
        os.mkdir("Gamma")
    os.chdir("Gamma")
    
    g=g1; trl_no=1;a=0; b=0; nu_data=[]
    
    for x in range(0, post_transcript[:,1].size):
        
        if(g != post_transcript[x,1]):
            b=x;
            denom = float(np.sum(post_transcript[a:b,4])); denom2= 2 - (1/post_transcript[x,0]);
            s_nsp= np.multiply(post_transcript[a:b,3], post_transcript[a:b,4])
            denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
            s2_nsp = np.multiply(post_transcript[a:b,3], s_nsp)
            numer = float(np.sum(s2_nsp))
            nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
            
            g = post_transcript[x,1]
            trl_no=1; a=x;
            
        elif(trl_no != post_transcript[x,2]):
            b=x;
            denom = float(np.sum(post_transcript[a:b,4])); denom2= 2 - (1/post_transcript[x,0]);
            s_nsp= np.multiply(post_transcript[a:b,3], post_transcript[a:b,4])
            denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
            s2_nsp = np.multiply(post_transcript[a:b,3], s_nsp)
            numer = float(np.sum(s2_nsp))
            nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
            
            trl_no= post_transcript[x,2]; a=x;
            
        elif(x == post_transcript[:,1].size -1):
            #last entry in series
            b=x+1
            denom = float(np.sum(post_transcript[a:b,4])); denom2= 2 - (1/post_transcript[x,0]);
            s_nsp= np.multiply(post_transcript[a:b,3], post_transcript[a:b,4])
            denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
            s2_nsp = np.multiply(post_transcript[a:b,3], s_nsp)
            numer = float(np.sum(s2_nsp))
            nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
            break
        
    print(" L ,  # ,  <S[p]\t")
    
    new_data= np.array(nu_data)
    
    zerodark30= pan.DataFrame(new_data, columns= ["L", "$Trial Number$", r"$\langle S[p] \rangle$", r"$\langle S'[p] \rangle$"])
    
    x1= np.transpose(new_data[:,0])
    x2= np.transpose(new_data[:,2])
    
    g= sea.lineplot(data=zerodark30, x="L" , y=r"$\langle S[p] \rangle$", estimator='mean', ci='sd', marker="s", err_style="band")
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([7.5, 0.15]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of p_avg - p_c:\t" +str(perr[1]))
    
    tukan= (popt[0], popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ \langle S[p] \rangle = %5.4f \times L^{(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.xlim(math.exp(3),g2+20)
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    plt.legend()
    g.set_title(r'$p = %f \quad ( \xi \longrightarrow \infty ) $' %(post_transcript[0,0]))
    plt.savefig("Log Band Gamma S(p) vs L (p--%8.7f) (Range-- %d-%d).png" %(post_transcript[0,0], g1, g2), dpi=400)
    plt.show()
    plt.close()
    
    
    os.chdir(r"..\..\..\..\..\analysis\Mass Action")
    #Returning to our home directory.
    
def plt_nuu(nudata):
    # Plots average cluster size data using Seaborn and stores it to the relevant directory.
    
    p_c = 0.728
    
    p = nudata[0,0] #Occupation Probability value at which simulations were run.
    g1 = int(nudata[0,1]) #Starting Grid Size
    g2 = int(nudata[-1,1])    #Ending Grid Size.
    
    print("%f \t G1---- %d \t G2---- %d" %(p,g1,g2))
    
    os.chdir(r"..\..\figures")
    # Changing to relevant directory.
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("DP")==False):
        os.mkdir("DP")
    os.chdir("DP")
    if(os.path.isdir("Finite Scaling")==False):
        os.mkdir("Finite Scaling")
    os.chdir("Finite Scaling")
    if(os.path.isdir("Nuu")==False):
        os.mkdir("Nuu")
    os.chdir("Nuu")
    
    
    g=g1; a=0; b=0; gingerman=[]
    
    for x in range(0, nudata[:,1].size):
        
        if(g != nudata[x,1]):
            b=x;
            print("For size %f, we have b = %d  and  b - a = %d" %(g, b, (b-a)))
            mean_p= np.mean(nudata[a:b,3])
            mean_p2= np.mean(nudata[a:b,4])
            sd_p = np.std(nudata[a:b,3])
            gingerman.append([g, mean_p, mean_p2, sd_p])
            
            a=x; g = nudata[x,1]
        elif(x == nudata[:,1].size -1):
            #last entry in series
            b=x+1
            print("For size %f, we have b = %d  and  b - a = %d" %(nudata[x,1], b, (b-a)))
            mean_p= np.mean(nudata[a:b,3])
            mean_p2= np.mean(nudata[a:b,4])
            sd_p = np.std(nudata[a:b,3])
            gingerman.append([g, mean_p, mean_p2, sd_p])
            break
    
    npframe = np.array(gingerman)
    print(npframe.shape)
    npframe[:,1] -= p_c
    npframe[:,1] = np.fabs(npframe[:,1])
    
    hurtlocker= pan.DataFrame(npframe, columns= [ "L",  r"$ | \langle p \rangle - p_c | $", r"$ \langle p^{2} \rangle $", r"$ \sigma_{p} $"])
    
    
    x1 =np.transpose(npframe[:,0])
    x2= np.transpose(npframe[:,3])
    
    g= sea.scatterplot(data=hurtlocker, x="L" , y= r"$ \sigma_{p} $")
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.1, -0.75]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of Sigma_p:\t" +str(perr[1]))
    
    tukan= (popt[0], -popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ \sigma_{p} = %5.4f \times L^{-(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.xlim(g1- 10, g2 + 10)
    
    plt.legend()
    g.set_title(r'$ \sigma_{p} \quad vs \quad L$')
    plt.savefig("Nu Sigma_p vs L (G1--%d G2-- %d).png" %(g1, g2), dpi=400)
    plt.show()
    plt.close()
    
    x2= np.transpose(npframe[:,1])
    
    g= sea.scatterplot(data=hurtlocker, x="L" , y= r"$ | \langle p \rangle - p_c | $")
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.05, -0.75]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of p_avg - p_c:\t" +str(perr[1]))
    
    tukan= (popt[0], -popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ | \langle p \rangle - p_c | = %5.4f \times L^{-(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.xlim(g1- 10, g2 + 10)
    
    plt.legend()
    g.set_title(r'$ | \langle p \rangle - p_c | \quad vs \quad L $')
    plt.savefig("Nu p_avg - p_c vs L (G1--%d G2-- %d).png" %(g1, g2), dpi=400)
    plt.show()
    plt.close()
            
            
            
        
        
                
def exxon_split(unspliced_data):
    # Arranges unspliced_data in strict ascending order of grid sizes.
        
    '''But first, we need to have a comprehensive list of all grid sizes in ascending order)'''
        
    a=unspliced_data[0,1]; L=[a] #Initialsing variables so as to detect all grid sizes.
        
    for x in range(0,len(unspliced_data[:,1])):
        #Iterating over all possible grid values to create a list of grid sizes in ascending  order.
        b = unspliced_data[x,1]
        if( b > a and (b not in L)):
            print("Bonobo:\t %3.0f" %(b))
            # A new grid size has been detected.
            L.append(b)
            
            
    '''Now for each grid size, all the revelant data from unspliced_data must be spliced out and concatanated into a 
    new array'''
    
    a=0; b=0 #Stores relevant splices for each grid size.
    
    m_splice =[]        
    for l in L:
        #Iterating over all the grid sizes, in unspliced_data
        flag =0; a=0; b=0
        for x in range(0,len(unspliced_data[:,1])):
            #Iterating over unspliced_data.
            if(l == unspliced_data[x,1] and flag==0): 
                # We have a new hit for the given grid size "l".
                a=x
                flag=1
            elif(unspliced_data[x,1] != l and flag==1):
                # The splice for the given grid size "l" just ended and we must extract the relevant slice.
                
                b=x; flag=0
                print("Slice for grid of size %d is:\t [ %d , %d ]" %(int(l), int(a), int(b)))
                if (len(m_splice) == 0):
                    #First one in the bag.
                    m_splice = unspliced_data[a:b,:]
                else:
                    m_splice = np.concatenate((m_splice, unspliced_data[a:b,:]), axis=0)
                    
            if( x == len(unspliced_data[:,1])-1 and flag == 1):
                #Special case that only applies to very last row of unspliced_data.
                b= x+1; flag=0
                print("Slice for grid of size %d is:\t [ %d , %d ]" %(int(l), int(a), int(b)))
                m_splice = np.concatenate((m_splice, unspliced_data[a:b,:]), axis=0)
                
                
                
            
            
    return m_splice;       
        
    
        
        
    
starter_pack()