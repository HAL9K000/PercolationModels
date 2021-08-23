# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 00:14:24 2021

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

''' This is a script that is specifically designed to ascertain, read and subsequently operate on data present
    in an arbitrary number of CSV files scattered across various folders.'''
    
    
def ticks(y, pos):
    return r'$e^{:.0f}$'.format(np.log(y))


def ln_pow_law(x, a, tau):
    return a - tau*x

def pow_law(x, a, expo):
    return a*(np.power(x, expo))

def crt_exp_gamma_beta_1stMar():
    
    ''' The CSV data is POORLY FORMATED & ERRONEOUSLY grouped into the following columns:
        The first column contains the occupation probability ( p ----> p_c(q=0) = 0.728 for DP class) at which the simulations were run.
        The second column contains BOTH the q-value as well as the grid size of the lattices  ("0.03 64" for example)
        while the third stores the current simulation number.
        The fourth and fifth column usually store the cluster size and the corresponding n_s(p) value for a given simulation,
        with the exception of the last entry into each of these columns for each trial number, which stores -10 (indicative of a
        spanning cluster) and P[p] respectively.
        '''
        
    unspliced_data = []
    #This will store the raw unspliced data without any processing from all the CSV files.
    
    q_L_unspl=[]; tr_unspl=[]
    #ThIS will store the last (q,L, #) value-pair present in each CSV file.

    #for i in range(0,9):
    base_path = r"1stMar" #+ "\\" + str(i)
    #print(base_path)
        
    files = glob.glob(base_path + "/*.csv", recursive=True)
    for file in files:
        if (os.path.getsize(file) > 512):
            #Filtering for file sizes that are greater than 512 Bytes in size.
            print(file)
            
            raw_colnames = ['p', 'q_L' , 'Tr No', 's', 'n_s(p,q)']
            raw_df= pan.read_csv('%s' %(file), skiprows=1, names=raw_colnames)
            
            print(raw_df)
            
            err_list= raw_df.q_L.to_list() #Storing column borne under "q_L" as a list.
            
            q=[]; L=[]
            for x in err_list:
                #print(x);
                xsplit = str(x).split('  ') #Split string by whitespace.
                #print(xsplit[0]); print(xsplit[1])
                q.append(float(xsplit[0]))
                L.append(float(xsplit[1]))
            
            roG_col= ['p', 'Tr No', 's', 'n_s(p,q)']
            data_temp = raw_df[roG_col].to_numpy()
            
            #print(data_temp)
            
            
            #Pruning the errors in the data.
            data_temp= np.insert(data_temp,1,q, axis=1); data_temp= np.insert(data_temp,2,L, axis=1);
            
            '''
            data_temp resembles:
                | p, q, L, #, s, n_s(p,q) |
            '''
            
            print(data_temp)
            #ko=input("Enter to continue:\t")
                
            if(len(unspliced_data) == 0):
                #First chip off the old block.
                unspliced_data = data_temp
                q_L_unspl.append((data_temp[-1,1],data_temp[-1,2],data_temp[-1,3])) #Appending end q, L, # values of the CSV

            else:
                for (a,b,c) in q_L_unspl:
                    if( (data_temp[-1,1],data_temp[-1,2]) ==(a,b)):
                        ''' The current CSV stores experiments performed using the same relevant parameters (specifically q & grid size)
                        as the previous CSV file that were processed.'''
                        data_temp[:,3] += c
                        print("Boris")
                        # We update the trial numbers in the new CSV file to be contiguous with the trial numbers from
                        # the previous CSVs, if the experimental data is the same.
                unspliced_data = np.concatenate((unspliced_data, data_temp), axis=0)
                q_L_unspl.append((unspliced_data[-1,1],unspliced_data[-1,2],unspliced_data[-1,3])) #Appending end q, L, # values of the CSV
            
    m_data = exxon_split(unspliced_data)
    # Arranges unspliced_data in strict ascending order of grid sizes.
    m_data = np.array(m_data)
    
    for i in range(0, 10000, 200):
        print(m_data[i,:])
        
    ko=input("Enter to continue:\t")
    
    #s= input("Enter any key to continue.")
    s=0
    '''for  x in m_data:
        if(s % 25 == 0):
            print("%d \t %d \t %f \t %f" %(int(x[1]), int(x[2]), x[3], x[4]))
        s+=1
    print("Total:\t %d" %(s))'''
    
    ''' Now to thresh out the percolation strength data from m_data.
        Remember, P[p] data is characterised by a -10 entry in the 4th column.'''
        
    split_data = m_data[:,4] == -10
    
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
    
def crt_exp_gamma_beta_17Mar():
    
    ''' The CSV data is PROPERLY FORMATED grouped into the following columns:
        The first column contains the occupation probability ( p ----> p_c(q=0) = 0.728 for DP class) at which the simulations were run.
        The second column contains the q-value, the third column represents the grid size of the lattices
        while the fourth stores the current simulation number.
        The fifth and sixth column usually store the cluster size and the corresponding n_s(p) value for a given simulation,
        with the exception of the last entry into each of these columns for each trial number, which stores -10 (indicative of a
        spanning cluster) and P[p] respectively.
        '''
        
    unspliced_data = []
    #This will store the raw unspliced data without any processing from all the CSV files.
    
    q_L_unspl=[]; tr_unspl=[]
    #ThIS will store the last (q,L, #) value-pair present in each CSV file.

    #for i in range(0,9):
    base_path = r"1stMar\17thMar" #+ "\\" + str(i)
    #print(base_path)
        
    files = glob.glob(base_path + "/*.csv", recursive=True)
    for file in files:
        if (os.path.getsize(file) > 512):
            #Filtering for file sizes that are greater than 512 Bytes in size.
            print(file)
            
            data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
            
            '''
            data_temp resembles:
                | p, q, L, #, s, n_s(p,q) |
            '''
            
            print(data_temp)
            #ko=input("Enter to continue:\t")
                
            if(len(unspliced_data) == 0):
                #First chip off the old block.
                unspliced_data = data_temp
                q_L_unspl.append((data_temp[-1,1],data_temp[-1,2],data_temp[-1,3])) #Appending end q, L, # values of the CSV

            else:
                for (a,b,c) in q_L_unspl:
                    if( (data_temp[-1,1],data_temp[-1,2]) ==(a,b)):
                        ''' The current CSV stores experiments performed using the same relevant parameters (specifically q & grid size)
                        as the previous CSV file that were processed.'''
                        data_temp[:,3] += c
                        print("Boris")
                        # We update the trial numbers in the new CSV file to be contiguous with the trial numbers from
                        # the previous CSVs, if the experimental data is the same.
                unspliced_data = np.concatenate((unspliced_data, data_temp), axis=0)
                q_L_unspl.append((unspliced_data[-1,1],unspliced_data[-1,2],unspliced_data[-1,3])) #Appending end q, L, # values of the CSV
            
    m_data = exxon_split(unspliced_data)
    # Arranges unspliced_data in strict ascending order of grid sizes.
    m_data = np.array(m_data)
    
    for i in range(0, 10000, 200):
        print(m_data[i,:])
        
    ko=input("Enter to continue:\t")
    
    #s= input("Enter any key to continue.")
    s=0
    '''for  x in m_data:
        if(s % 25 == 0):
            print("%d \t %d \t %f \t %f" %(int(x[1]), int(x[2]), x[3], x[4]))
        s+=1
    print("Total:\t %d" %(s))'''
    
    ''' Now to thresh out the percolation strength data from m_data.
        Remember, P[p] data is characterised by a -10 entry in the 4th column.'''
        
    split_data = m_data[:,4] == -10
    
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
    
def plt_beta(perc_data):
    # Plots percolation data using Seaborn and stores it to the relevant directory.
    
    p = perc_data[0,0] #Occupation Probability value at which simulations were run.
    g1 = int(perc_data[0,2]) #Starting Grid Size
    g2 = int(perc_data[-1,2])    #Ending Grid Size.
    q1 = int(perc_data[0,1]) #Starting q
    q2 = int(perc_data[-1,1])    #Ending q.
    
    print("%f \t G1---- %d \t G2---- %d \t q1---- %d \t q2---- %d" %(p,g1,g2, q1, q2))
    
    os.chdir(r"..\..\..\figures")
    # Changing to relevant directory.
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("TCP")==False):
        os.mkdir("TCP")
    os.chdir("TCP")
    if(os.path.isdir("Finite Scaling")==False):
        os.mkdir("Finite Scaling")
    os.chdir("Finite Scaling")
    if(os.path.isdir("Beta")==False):
        os.mkdir("Beta")
    os.chdir("Beta")
    
    b=perc_data[0,1]; Q=[b] #Initialsing variables so as to detect all q values.
        
    for x in range(0,len(perc_data[:,2])):
        #Iterating over all possible grid values to create a list of grid sizes in ascending  order.
        b= perc_data[x,1]
        if( (b not in Q)):
            print("Gibbon:\t %3.3f" %(b))
            # A new q value has been detected.
            Q.append(b)
    
    Q.sort()
    
    for q in Q:
        
        split_data = perc_data[:,1] == q
    
        perc_data_q = perc_data[split_data]     #Final percolation data for q.
        p= perc_data_q[0,0]
    
        hurtlocker= pan.DataFrame(perc_data_q, columns= ["p", "q", "L", "Trial Number",  "-10", r"P[p,q]"])
    
        x1 =np.transpose(perc_data_q[:,2])
        x2= np.transpose(perc_data_q[:,5])
    
        g= sea.lineplot(data=hurtlocker, x="L" , y="P[p,q]", estimator='mean', ci='sd', marker="s", err_style="band")
    
        popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.3, -0.005]))
    
        perr = np.sqrt(np.diag(pcov))
        
        print("SD of C:\t" +str(perr[1]) + " for p:\t" +str(p) + " at q:\t" +str(q))
    
        tukan= (popt[0], -popt[1], perr[1])
        plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ P[p,q] = %5.4f \times L^{-(%5.4f \mp %5.4f)} $ ' % tukan )
        plt.legend()
        plt.xlim(g1,g2+20)
        plt.yscale('log', basey= math.e)
        plt.xscale('log', basex= math.e)
        g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
        g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    
        g.set_title(r'$p = %8.7f ( \xi \longrightarrow \infty ), q = %3.2f $' %(p,q))
        plt.savefig("TCP Exp Beta P(p,q) vs L (p--%8.7f) (q--%3.2f) (Range-- %d-%d).png" %(p,q, g1, g2), dpi=400)
        plt.show()
        plt.close()
    
    
    
    os.chdir(r"..\..\..\..\..\analysis\Mass Action\TCP")
    #Returning to our home directory.
    
     
def plt_gamma(post_transcript):
    # Plots average cluster size data using Seaborn and stores it to the relevant directory.
    
    p = post_transcript[0,0] #Occupation Probability value at which simulations were run.
    g1 = int(post_transcript[0,2]) #Starting Grid Size
    g2 = int(post_transcript[-1,2])    #Ending Grid Size.
    q1 = int(post_transcript[0,1]) #Starting q
    q2 = int(post_transcript[-1,1])    #Ending q.
    
    print("%f \t G1---- %d \t G2---- %d \t q1---- %d \t q2---- %d" %(p,g1,g2, q1, q2))
    
    os.chdir(r"..\..\..\figures")
    # Changing to relevant directory.
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("TCP")==False):
        os.mkdir("TCP")
    os.chdir("TCP")
    if(os.path.isdir("Finite Scaling")==False):
        os.mkdir("Finite Scaling")
    os.chdir("Finite Scaling")
    if(os.path.isdir("Gamma")==False):
        os.mkdir("Gamma")
    os.chdir("Gamma")
    
    b=post_transcript[0,1]; Q=[b] #Initialsing variables so as to detect all q values.
        
    for x in range(0,len(post_transcript[:,2])):
        #Iterating over all possible grid values to create a list of grid sizes in ascending  order.
        b= post_transcript[x,1]
        if( (b not in Q)):
            print("Gibbon:\t %3.3f" %(b))
            # A new q value has been detected.
            Q.append(b)
    
    Q.sort()
    
    for q in Q:
        
        split_data = post_transcript[:,1] == q
    
        post_transcript_q = post_transcript[split_data]     #Final ns_p for q.
        
        p= post_transcript_q[0,0]
    
        g=g1; trl_no=1;a=0; b=0; nu_data=[]
    
        for x in range(0, post_transcript_q[:,1].size):
        
            if(g != post_transcript_q[x,2]):
                b=x;
                denom = float(np.sum(post_transcript_q[a:b,5])); denom2= 2 - (1/post_transcript_q[x,0]);
                s_nsp= np.multiply(post_transcript_q[a:b,4], post_transcript_q[a:b,5])
                denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
                print("Denom:\t" +str(denom)+ " Denom2:\t" +str(denom2))
                s2_nsp = np.multiply(post_transcript_q[a:b,4], s_nsp)
                numer = float(np.sum(s2_nsp))
                nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
            
                g = post_transcript_q[x,2]
                trl_no=1; a=x;
            
            elif(trl_no != post_transcript_q[x,3]):
                b=x;
                denom = float(np.sum(post_transcript_q[a:b,5])); denom2= 2 - (1/post_transcript_q[x,0]);
                s_nsp= np.multiply(post_transcript_q[a:b,4], post_transcript_q[a:b,5])
                denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
                print("Denom:\t" +str(denom)+ " Denom2:\t" +str(denom2))
                s2_nsp = np.multiply(post_transcript_q[a:b,4], s_nsp)
                numer = float(np.sum(s2_nsp))
                nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
            
                trl_no= post_transcript_q[x,3]; a=x;
            
            elif(x == post_transcript_q[:,1].size -1):
                #last entry in series
                b=x+1
                denom = float(np.sum(post_transcript_q[a:b,5])); denom2= 2 - (1/post_transcript_q[x,0]);
                s_nsp= np.multiply(post_transcript_q[a:b,4], post_transcript_q[a:b,5])
                denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
                s2_nsp = np.multiply(post_transcript_q[a:b,4], s_nsp)
                numer = float(np.sum(s2_nsp))
                nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
                break
        
        print(" L ,  # ,  <S[p]\t")
    
        new_data= np.array(nu_data)
    
        zerodark30= pan.DataFrame(new_data, columns= ["L", "$Trial Number$", r"$\langle S[p,q] \rangle$", r"$\langle S'[p,q] \rangle$"])
    
        x1= np.transpose(new_data[:,0])
        x2= np.transpose(new_data[:,2])
        
        g= sea.lineplot(data=zerodark30, x="L" , y=r"$\langle S[p,q] \rangle$", estimator='mean', ci='sd', marker="s", err_style="band")
        popt=0; pcov=0;
        if(q != 0.06):
            popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([7.5, 0.15]))
        else:
            popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.1, 1.78]), bounds=((0,1.66), (0.6, 1.78)))
    
        perr = np.sqrt(np.diag(pcov))
        
        print("SD of p_avg - p_c:\t" +str(perr[1]))
    
        tukan= (popt[0], popt[1], perr[1])
        plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ \langle S[p,q] \rangle = %5.4f \times L^{(%5.4f \mp %5.4f)} $ ' % tukan )
        plt.xlim(g1-5, g2 + 10)
        plt.yscale('log', basey= math.e)
        plt.xscale('log', basex= math.e)
        g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
        g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
        plt.legend()
        g.set_title(r'$p = %f \quad ( \xi \longrightarrow \infty ), q = %3.2f $' %(post_transcript_q[0,0],q))
        if(q== 0.06):
            plt.savefig("TCP %7.5f Exp Gamma S(p,q) vs L (p--%8.7f) (q--%3.2f) (Range-- %d-%d).png" %(popt[1], post_transcript_q[0,0],q, g1, g2), dpi=400)
        else:
            plt.savefig("TCP Exp Gamma S(p,q) vs L (p--%8.7f) (q--%3.2f) (Range-- %d-%d).png" %(post_transcript_q[0,0],q, g1, g2), dpi=400)
        plt.show()
        plt.close()
    
    
    os.chdir(r"..\..\..\..\..\analysis\Mass Action")
    #Returning to our home directory.
    
def exxon_split(unspliced_data):
    # Arranges unspliced_data in strict ascending order of grid sizes.
        
    '''But first, we need to have a comprehensive list of all grid sizes in ascending order)'''
        
    a=unspliced_data[0,2]; L=[a] #Initialsing variables so as to detect all grid sizes.
    b=unspliced_data[0,1]; Q=[b] #Initialsing variables so as to detect all q values.
        
    for x in range(0,len(unspliced_data[:,2])):
        #Iterating over all possible grid values to create a list of grid sizes in ascending  order.
        a = unspliced_data[x,2]; b= unspliced_data[x,1]
        if( (a not in L)):
            print("Bonobo:\t %3.0f" %(a))
            # A new grid size has been detected.
            L.append(a)
        if( (b not in Q)):
            print("Gibbon:\t %3.0f" %(b))
            # A new q value has been detected.
            Q.append(b)
    
    L.sort(); Q.sort()
            
            
    '''Now for each grid size, all the revelant data from unspliced_data must be spliced out and concatanated into a 
    new array'''
    
    a=0; b=0 #Stores relevant splices for each grid size.
    
    m_splice =[]
    for q in Q:        
      for l in L:
        #Iterating over all the grid sizes, in unspliced_data
        flag =0; a=0; b=0
        for x in range(0,len(unspliced_data[:,2])):
            #Iterating over unspliced_data.
            if(l == unspliced_data[x,2] and q == unspliced_data[x,1] and flag==0): 
                # We have a new hit for the given grid size "l".
                a=x
                flag=1
            elif(unspliced_data[x,2] != l and q == unspliced_data[x-1,1] and flag==1):
                # The splice for the given grid size "l" just ended and we must extract the relevant slice.
                
                b=x; flag=0
                print("Slice at q value %f for grid of size %d is:\t [ %d , %d ]" %(float(q), int(l), int(a), int(b)))
                if (len(m_splice) == 0):
                    #First one in the bag.
                    m_splice = unspliced_data[a:b,:]
                else:
                    m_splice = np.concatenate((m_splice, unspliced_data[a:b,:]), axis=0)
                    
            if( x == len(unspliced_data[:,1])-1 and q == unspliced_data[x,1] and flag == 1):
                #Special case that only applies to very last row of unspliced_data.
                b= x+1; flag=0
                print("Slice at q value %f for grid of size %d is:\t [ %d , %d ]" %(float(q), int(l), int(a), int(b)))
                m_splice = np.concatenate((m_splice, unspliced_data[a:b,:]), axis=0)
                
                
                
            
            
    return m_splice;

crt_exp_gamma_beta_1stMar()