# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 09:48:40 2020

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


def ticks(y, pos):
    return r'$e^{:.0f}$'.format(np.log(y))


def ln_pow_law(x, a, tau):
    return a - tau*x

def pow_law(x, a, expo):
    return a*(np.power(x, expo))

def plotter():
    
    s = input("Enter Critical Exponent You Would Like To Plot (Choice of Beta, Gamma, Sigma, Tau, Nuu (Nu) or Fin Sc (Note Spelling):\t")
    
    if ( s== "Gamma" or s == "gamma"):
        crtexpt_gamma()
    elif( s== "Tau" or s == "tau"):
        crtexpt_tau()
    elif ( s== "Beta" or s == "beta"):
        crtexpt_beta()
    elif ( s=="Sigma" or s== "sigma"):
        crtexpt_sigma()
    elif (s=="Nuu" or s== "nuu" or s== "nu" or s== "Nu"):
        crtexpt_nuu()
    elif (s=="Fin Sc" or s== "fin sc"):
        crtexpt_finsc()
    else:
        print("Your input didn't match any known exponent. Learn how to spell properly, dirtwad.")
    
        

def crtexpt_gamma():
    os.chdir("..\simulations\CrtExp")
    # Changing to relevant directory.
    
    #s = input("The default file to be analysed is file: 'Gamma_NP_L_256_p1_0.578_p2_0.608_R_5_Cen_25.csv'. To modify some other Gamma file, type 'Y':\t")
    
    s = input("The default file to be analysed is file: 'FinScGam_SP_p_0.5927_Div_241_G1_25_G2_401.csv'. To modify some other Gamma file, type 'Y':\t")
    
    p =0.5927; g1 = 25; g2 = 401; R= 5; Div = 241;
    
    if (s == "Y" or s == "y"):
        p= float(input("Enter Grid Size:\t"))
        g1= int(input("Enter p1 (starting value of p):\t"))
        g2= int(input("Enter p2 (ending value of p):\t"))
        R= int(input("Enter number of random trials:\t"))
        Div= int(input("Enter number of divisions:\t"))
        
    data= np.genfromtxt('FinScGam_SP_p_%5.4f_Div_%d_G1_%d_G2_%d.csv' %(p, Div, g1, g2), delimiter=",", comments='#', skip_header=1)
                        
    os.chdir(r"..\..\figures")
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("Finite Scaling")==False):
        os.mkdir("Finite Scaling")
    os.chdir("Finite Scaling")
    if(os.path.isdir("Gamma")==False):
        os.mkdir("Gamma")
    os.chdir("Gamma")

    # p values are such that p ----> p_c
    p_c =0.592746; rtot = Div*R;
    
    '''data[:,1] -= p_c        #From p, creating p -p_c.
    
    split_data = data[:,1] > 0
    
    datahigh = data[split_data]
    datalow = data[~split_data]
    
    datalow[:,1] =np.abs(datalow[:,1]) # Generating |p - p_c|
    
    #data[:,1] = np.abs(data[:,1]) # Generating |p - p_c|'''
    
    g=g1; trl_no=1;a=0; b=0; nu_data=[]
    
    for x in range(0, data[:,1].size):
        
        if(g != data[x,1]):
            b=x;
            denom = float(np.sum(data[a:b,4])); denom2= data[x,0];
            s_nsp= np.multiply(data[a:b,3], data[a:b,4])
            denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
            s2_nsp = np.multiply(data[a:b,3], s_nsp)
            numer = float(np.sum(s2_nsp))
            nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
            
            g = data[x,1]
            trl_no=1; a=x;
            
        elif(trl_no != data[x,2]):
            b=x;
            denom = float(np.sum(data[a:b,4])); denom2= data[x,0];
            s_nsp= np.multiply(data[a:b,3], data[a:b,4])
            denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
            s2_nsp = np.multiply(data[a:b,3], s_nsp)
            numer = float(np.sum(s2_nsp))
            nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
            
            trl_no= data[x,2]; a=x;
            
        elif(x == data[:,1].size -1):
            #last entry in series
            b=x+1
            denom = float(np.sum(data[a:b,4])); denom2= data[x,0];
            s_nsp= np.multiply(data[a:b,3], data[a:b,4])
            denom= float(np.sum(s_nsp)) #Calculating S[p] for a given trial number as per defn.
            s2_nsp = np.multiply(data[a:b,3], s_nsp)
            numer = float(np.sum(s2_nsp))
            nu_data.append([g, trl_no, (numer/denom), (numer/denom2)])
            break
        
    print(" L ,  # ,  <S[p]\t")
    for x in nu_data:
        print(x)
    
    new_data= np.array(nu_data)
    
    zerodark30= pan.DataFrame(new_data, columns= ["L", "$Trial Number$", r"$\langle S[p] \rangle$", r"$\langle S'[p] \rangle$"])
    
    x1= np.transpose(new_data[:,0])
    x2= np.transpose(new_data[:,2])
    
    g= sea.lineplot(data=zerodark30, x="L" , y=r"$\langle S[p] \rangle$", estimator='mean', ci='sd')
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([7.5, 0.15]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of p_avg - p_c:\t" +str(perr[1]))
    
    tukan= (popt[0], popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ \langle S[p] \rangle = %5.4f \times L^{(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.xlim(math.exp(3), g2 + 10)
    #plt.ylim(math.exp(0), math.exp(9))
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    plt.legend()
    g.set_title(r'$p = %f \quad ( \xi \longrightarrow \infty ) $' %(data[0,0]))
    plt.savefig("Exp Fig V S(p) vs L (p--%8.7f).png" %(data[0,0]), dpi=400)
    plt.show()
    plt.close()
    
    
    '''zerodark30= pan.DataFrame(datahigh, columns= ["Trial Number", "$p - p_{c}$", r"$\langle S(p) \rangle$"])
    
    g=sea.scatterplot(data=zerodark30, x="$p - p_{c}$" , y=r"$\langle S(p) \rangle$", hue="Trial Number", hue_order= ["1","25","50", "75", "100", "125"], palette= "viridis")
    g.set_title('$p > p_{c}$')
    plt.savefig("Scatter High Hue (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    #plt.flush()
    plt.close()
    
    g= sea.lineplot(data=zerodark30, x="$p - p_{c}$" , y=r"$\langle S(p) \rangle$", estimator='mean', ci='sd')
    g.set_title('$p > p_{c}$')
    plt.savefig("Line High Hue (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    plt.close()
    
    zerodark30= pan.DataFrame(datalow, columns= ["Trial Number", "$|p - p_{c}|$", r"$\langle S(p) \rangle$"])
    
    g=sea.scatterplot(data=zerodark30, x="$|p - p_{c}|$" , y=r"$\langle S(p) \rangle$", hue="Trial Number", hue_order= ["1","25","50", "75", "100", "125"], palette= "viridis")
    g.set_title('$p < p_{c}$')
    plt.savefig("Scatter Low Hue (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    #plt.flush()
    plt.close()
    
    g= sea.lineplot(data=zerodark30, x="$|p - p_{c}|$" , y=r"$\langle S(p) \rangle$", estimator='mean', ci='sd')
    g.set_title('$p < p_{c}$')
    plt.savefig("Line Low Hue (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    plt.close()
    
    
    
    # On log scale below.
    
    ln_p_pc = np.log(datahigh[:,1])
    ln_Sp = np.log(datahigh[:,2])
    
    lnPframe= np.zeros((datahigh[:,0].size, 3))
    
    lnPframe[:,0]= datahigh[:,0]; lnPframe[:,1]= ln_p_pc ; lnPframe[:,2]= ln_Sp; 
    
    hurtlocker= pan.DataFrame(lnPframe, columns= ["Trial Number", "$ln(p - p_{c})$", r"$\langle ln(S(p)) \rangle$"])'''
    #Creating Pandas Dataframe.
    
    '''g=sea.scatterplot(data=hurtlocker, x="$ln|p - p_{c}|$" , y=r"$\langle ln(S(p)) \rangle$", hue="Trial Number", hue_order= ["1","25","50", "75", "100", "125"], palette= "viridis")
    
    plt.savefig("Scatter Hue Ln (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.close()'''
    
    '''g= sea.lineplot(data=hurtlocker, x="$ln(p - p_{c})$" , y=r"$\langle ln(S(p)) \rangle$", estimator='mean', ci='sd')
    g.set_title(r"$p > p_{c}$")
    plt.savefig("Line Hue Ln High (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    plt.close()
    
    ln_p_pc = np.log(datalow[:,1])
    ln_Sp = np.log(datalow[:,2])
    
    lnPframe= np.zeros((datalow[:,0].size, 3))
    
    lnPframe[:,0]= datalow[:,0]; lnPframe[:,1]= ln_p_pc ; lnPframe[:,2]= ln_Sp; 
    
    hurtlocker= pan.DataFrame(lnPframe, columns= ["Trial Number", "$ln|p - p_{c}|$", r"$\langle ln(S(p)) \rangle$"])
    #Creating Pandas Dataframe.
    
    g= sea.lineplot(data=hurtlocker, x="$ln|p - p_{c}|$" , y=r"$\langle ln(S(p)) \rangle$", estimator='mean', ci='sd')
    g.set_title(r"$p < p_{c}$")
    plt.savefig("Line Hue Ln Low (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    plt.close()'''
    
    
    
    
    
    
    
    

def crtexpt_beta():
    os.chdir("..\simulations\CrtExp")
    # Changing to relevant directory.
    
    s = input("The default file to be analysed is file: 'Beta_NP_L_256_p1_0.593_p2_0.613_Cen_20_R_5.csv'. To modify some other Beta file, type 'Y':\t")
    
    GrS =256; p1 = 0.593; p2 = 0.613; Cen = 20; R= 5;
    
    if (s == "Y" or s == "y"):
        GrS= int(input("Enter Grid Size:\t"))
        p1= float(input("Enter p1 (starting value of p):\t"))
        p2= float(input("Enter p2 (ending value of p):\t"))
        R= int(input("Enter number of random trials:\t"))
        Cen= int(input("Enter number of censuses:\t"))
        
    data= np.genfromtxt('Beta_SP_L_%d_p1_%4.3f_p2_%4.3f_Cen_%d_R_%d.csv' %(GrS, p1, p2, Cen, R), delimiter=",", comments='#', skip_header=2)
    
    # p values are such that p ----> p_c+
    p_c =0.592746
    
    data[:,0] -= p_c        #From p, creating p -p_c = |p - p_c| as p> p_c
    
    ln_p_pc = np.log(data[:,0])
    ln_P = np.log(data[:,2])
    
    lnPframe= np.zeros((data[:,0].size, 2))
    
    lnPframe[:,0]= ln_p_pc ; lnPframe[:,1]= ln_P; 
    
    hurtlocker= pan.DataFrame(lnPframe, columns= ["$ln|p - p_{c}|$", "$ln(P(p))$"])
    #Creating Pandas Dataframe.
    
    os.chdir(r"..\..\figures")
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("Beta")==False):
        os.mkdir("Beta")
    os.chdir("Beta")
    
    
    
    g= sea.scatterplot(data=hurtlocker, x="$ln|p - p_{c}|$" , y="$ln(P(p))$")
    
    rtot = Cen*R
    
    plt.savefig("Scatter P(p) (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    #plt.flush()
    plt.close()
    
    g= sea.lineplot(data=hurtlocker, x="$ln|p - p_{c}|$" , y="$ln(P(p))$")
    '''plt.xlim(2,6)
    plt.xticks(np.arange(2, 6.0, step=0.20))
    plt.ylim(0,4)
    plt.yticks(np.arange(0, 4.2, step=0.20))'''
    
    plt.savefig("Line P(p) (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    #plt.flush()
    plt.close()
    
    
    '''split_data = lnPframe[:,0] >= -7
    filter_data = lnPframe[split_data] #Only considering ln|p - p_{c}| > 10^(-7)'''
    
    
    g= sea.scatterplot(data=hurtlocker, x="$ln|p - p_{c}|$" , y="$ln(P(p))$")
    
    x1 = np.transpose(ln_p_pc)
    x2 = np.transpose(ln_P)
    
    #print(ln_P.shape)
    
    popt, pcov = curve_fit(ln_pow_law, x1, x2, p0= np.asarray([0.3, -0.5]))
    
    perr = np.sqrt(np.diag(pcov))
    tukan= (popt[0], -popt[1], perr[1])
    
    plt.plot(x1, ln_pow_law(x1, *popt), 'm--', label=r'Th Fit: $ ln(P(p)) = %5.4f + (%5.4f \mp %4.3f) \times ln|p - p_{c}|$' % tukan)
    plt.legend()
    '''plt.xlim(-7.0,-3.0)
    plt.xticks(np.arange(-7.0, -3.0, step=0.5))
    plt.ylim(-2.0,-0.4)
    plt.yticks(np.arange(-2.0,-0.3, step=0.10))'''
    plt.savefig("Line Beta Fit (G--%d  N--%d).png" %(GrS, rtot), dpi=400)
    plt.show()
    plt.close()
    
    
    os.chdir(r"..\..\..\analysis")
    #Back to our home base.
    
    


def crtexpt_tau():
    os.chdir("..\simulations\CrtExp")
    # Changing to relevant directory.
    
    s = input("The default file to be analysed is file: 'Tau_NP_L_256_p1_0.552_p2_0.592_R_5_Cen_25.csv'. To modify some other Tau file, type 'Y':\t")
    
    GrS =256; p1 = 0.552; p2 = 0.592; R= 5; Cen = 25;
    
    if (s == "Y" or s == "y"):
        GrS= int(input("Enter Grid Size:\t"))
        p1= float(input("Enter p1 (starting value of p):\t"))
        p2= float(input("Enter p2 (ending value of p):\t"))
        R= int(input("Enter number of random trials:\t"))
        Cen= int(input("Enter number of censuses:\t"))
    
    data= np.genfromtxt('Tau_SP_L_%d_p1_%4.3f_p2_%4.3f_R_%d_Cen_%d.csv' %(GrS, p1, p2, R, Cen), delimiter=",", comments='#', skip_header=1)
                        
    p= float(input("Choose a p-value to analyse (b/w p1 & p2):\t"))
    
    a=0; b=0 #Slices for p.
    
    os.chdir(r"..\..\figures")
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("Tau")==False):
        os.mkdir("Tau")
    os.chdir("Tau")
    
    p_c = 0.592746
    lsquare = float(GrS*GrS)
    
    
    
    #The following loop finds relevant slices for p
    flag=0
    for x in range(0, data[:,1].size):
        
        if (data[x,1] == p and flag == 0):
            a = x
            flag=1
        elif (data[x,1] != p and flag == 1):
            #p slice just ended.
            b = x
            flag =2
        elif( x == data[:,1].size -1):
            #Last element.
            b = x+1
            flag =2
        elif (flag == 2):
            break
    
    s= data[a:b, 2]
    nsp = data[a:b, 3]
    #tr = data[a:b, 0]  #Trial number.
    
    npframe=data[a:b, :]
    
    npframe[:,3] /= lsquare
    hurtlocker= pan.DataFrame(npframe, columns= ["Trial Number", "p", "s", r"$\langle n_{s}(p) \rangle$"])
    #Creating Pandas Dataframe.
    
    g= sea.scatterplot(data=hurtlocker, x="s" , y=r"$\langle n_{s}(p) \rangle$")
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    
    plt.savefig("Scatter (p--%6.5f).png" %(p), dpi=400)
    plt.show()
    #plt.flush()
    plt.close()
    
    '''g = sea.scatterplot(data=hurtlocker, x="ln(s)" , y=r"$\langle ln(n_{s}(p)) \rangle$", hue="Trial Number", hue_order= ["1","25","50", "75", "100", "125"], palette= "viridis")
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig("Scatter Hue (p--%4.3f).png" %(p), dpi=400)
    plt.show()
    #plt.flush()
    plt.close()'''
    
    g= sea.lineplot(data=hurtlocker, x="s" , y=r"$\langle n_{s}(p) \rangle$", estimator='mean', ci='sd')
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    plt.xlim(math.exp(1),math.exp(6))
    plt.ylim(math.exp(-12),math.exp(-4))
    '''plt.xticks(np.arange(2, 6.1, step=0.20))
    plt.ylim(0,4)
    plt.yticks(np.arange(0, 4.2, step=0.20)) '''
    
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    
    plt.savefig("Line (p--%6.5f).png" %(p), dpi=400)
    
    a=[];b=[];c=[] #Used to store select trial numbers, ln(s) entries and ln(nsp) values from npframe.
    
    for x in range(0,s.size):
        if (math.exp(2) <= npframe[x,2] and npframe[x,2] <=math.exp(4)):
            a.append(npframe[x,0]) #Filtering relevant trial numbers.
            b.append(npframe[x,2]) #Filtering relevant (s) numbers.
            c.append(npframe[x,3]) #Filtering relevant (nsp) numbers.
            
    x1= np.asarray(b)
    x2= np.asarray(c)
    
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([1, -2]))
    
    perr = np.sqrt(np.diag(pcov))
    
    tukan= (popt[0], popt[1], perr[1])
    
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ \langle n_{s}(p) \rangle = %5.4f \times s^{(%5.4f \mp %4.3f)}$' % tukan )
    plt.legend()
    plt.xlim(math.exp(2),math.exp(5))
    plt.ylim(math.exp(-12),math.exp(-7))
    g.set_title('p = %6.5f' %(p))
    if p == p_c :
       g.set_title('$p = p_{c}$') 
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    plt.savefig("Line Tau Expo Fit (p--%6.5f).png" %(p), dpi=400)
    plt.show()
    plt.close()
    print("SD in tau: \t %f" %(perr[1]))
    
    
    
    
    #LOG PLOTS BELOW
    
    ln_s = np.log(s)
    nsp /= lsquare
    ln_nsp =np.log(nsp)
    
    npframe[:,2] = ln_s; npframe[:,3] = ln_nsp;
    
    hurtlocker= pan.DataFrame(npframe, columns= ["Trial Number", "p", "ln(s)", r"$\langle ln(n_{s}(p)) \rangle$"])
    #Creating Pandas Dataframe.
    
    g= sea.scatterplot(data=hurtlocker, x="ln(s)" , y=r"$\langle ln(n_{s}(p)) \rangle$")
    
    plt.savefig("Scatter Ln (p--%6.5f).png" %(p), dpi=400)
    plt.show()
    #plt.flush()
    plt.close()
    
    g = sea.scatterplot(data=hurtlocker, x="ln(s)" , y=r"$\langle ln(n_{s}(p)) \rangle$", hue="Trial Number", hue_order= ["1","25","50", "75", "100", "125"], palette= "viridis")
    
    plt.savefig("Scatter Hue Ln (p--%6.5f).png" %(p), dpi=400)
    plt.show()
    #plt.flush()
    plt.close()
    
    p_c= 0.592746; s_e = 0; sig = 36.0/91.0;
    if( p != p_c):
        s_e =math.pow(math.fabs((p - p_c)), -(1.0/sig))
        s_e = math.log(s_e)
    
    
    g= sea.lineplot(data=hurtlocker, x="ln(s)" , y=r"$\langle ln(n_{s}(p)) \rangle$", estimator='mean', ci='sd')
    plt.xlim(1,6)
    plt.xticks(np.arange(1, 6.2, step=0.5))
    plt.ylim(-23,-16)
    plt.yticks(np.arange(-23, -16, step=0.5)) 
    
    plt.savefig("Line Ln (p--%6.5f).png" %(p), dpi=400)
    
    
    #Now to plot the best fit power law curve.
    
    # We are choosing the range 2 <= ln(s) <= 3 as our plotting schema.
    
    a=[];b=[];c=[] #Used to store select trial numbers, ln(s) entries and ln(nsp) values from npframe.
    
    for x in range(0,ln_s.size):
        if (2 <= npframe[x,2] and npframe[x,2] <=3.8):
            a.append(npframe[x,0]) #Filtering relevant trial numbers.
            b.append(npframe[x,2]) #Filtering relevant ln(s) numbers.
            c.append(npframe[x,3]) #Filtering relevant ln(nsp) numbers.
            
    '''split_data = npframe[:,2] >= 2
    
    print("First twenty elements of split_data")
    for i in range(0,15):
        print("%f \t" %(split_data[i])),
    print("\n")'''
    
    
    '''plot_filt= np.zeros((len(b),3)) #NP array to store a,b,c as it's three columns.
    
    plot_filt[:,0] = np.reshape(a, (-1,1))
    plot_filt[:,1] = np.reshape(b, (-1,1))
    plot_filt[:,2] = np.reshape(c, (-1,1))'''
    
    x1= np.asarray(b)
    x2= np.asarray(c)
    
    popt, pcov = curve_fit(ln_pow_law, x1, x2, p0= np.asarray([6, 2]))
    
    perr = np.sqrt(np.diag(pcov))
    
    tukan= (popt[0], popt[1], perr[1])
    
    plt.plot(x1, ln_pow_law(x1, *popt), 'm--', label=r'Th Fit: $\langle ln_{s}(p) \rangle = %5.4f - (%5.4f \pm %4.3f) \times ln(s)$' % tukan )
    plt.legend()
    g.set_title('p = %6.5f' %(p))
    if p == p_c :
       g.set_title('$p = p_{c}$')
    plt.xlim(2,5)
    plt.xticks(np.arange(2, 5.1, step=0.20))
    plt.ylim(-23.5,-17)
    plt.yticks(np.arange(-23.5, -17, step=0.5)) 
    plt.savefig("Line Tau Fit Ln (p--%6.5f).png" %(p), dpi=400)
    plt.show()
    plt.close()
    print("SD in log plot tau:\t %f" %(perr[1]))
    
    
    os.chdir(r"..\..\..\analysis")
    #Back to our home base.
    
    
def crtexpt_sigma():
    
    os.chdir("..\simulations\CrtExp")
    # Changing to relevant directory.
    
    s = input("The default file to be analysed is file: 'Tau_NP_L_256_p1_0.552_p2_0.592_R_5_Cen_25.csv'. To modify some other Tau file, type 'Y':\t")
    
    GrS =256; p1 = 0.552; p2 = 0.592; R= 5; Cen = 25;
    
    if (s == "Y" or s == "y"):
        GrS= int(input("Enter Grid Size:\t"))
        p1= float(input("Enter p1 (starting value of p):\t"))
        p2= float(input("Enter p2 (ending value of p):\t"))
        R= int(input("Enter number of random trials:\t"))
        Cen= int(input("Enter number of censuses:\t"))
    
    data= np.genfromtxt('Tau_NP_L_%d_p1_%4.3f_p2_%4.3f_R_%d_Cen_%d.csv' %(GrS, p1, p2, R, Cen), delimiter=",", comments='#', skip_header=1)
                        
    #p= float(input("Choose a p-value to analyse (b/w p1 & p2):\t"))
    
    lsquare= float(256*256)
    data[:,3] /= lsquare
    
    os.chdir(r"..\..\figures")
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("Sigma")==False):
        os.mkdir("Sigma")
    os.chdir("Sigma")
    
    pend= float(input("Choose ending point of p (needs to be less than p2):\t"))
    
    p_st = p1
    
    a=0; b=0 #Slices for p.
    
    L= [p1, p1 +0.006, round(p1 + 0.012,3), round(p1 + 0.018,3)] 
    #Points at which to obtain plots
    
    C_fits =[]
    
    while (p_st <= pend):
        
        
        split_data = data[:,1] == p_st
        filter_data = data[split_data] #Taking the slice where p = p_st
        
        ln_nsp = np.log(filter_data[:,3])
        
        npframe = filter_data
        
        npframe[:,3] = ln_nsp 
        
        a=[];b=[];c=[] #Used to store select trial numbers, ln(s) entries and ln(nsp) values from npframe.
    
        for x in range(0,ln_nsp.size):
            if (math.exp(3.5) <= npframe[x,2] and npframe[x,2] <=math.exp(4.5)):
                a.append(npframe[x,0]) #Filtering relevant trial numbers.
                b.append(npframe[x,2]) #Filtering relevant (s) numbers.
                c.append(npframe[x,3]) #Filtering relevant (nsp) numbers.
            
        x1= np.asarray(b)
        x2= np.asarray(c)
        
        '''if(p_st == 0.563):
            print("X1:")
            for x in x1:
                print(str(x)+ "  ", end='')
            print("\n\nX2:")
            for x in x2:
                print(str(x)+ "  ", end='')
            hurtlocker= pan.DataFrame(npframe, columns= ["Trial Number", "p", "s", r"$\langle ln(n_{s}(p)) \rangle$"])
            g= sea.lineplot(data=hurtlocker, x="s" , y=r"$\langle ln(n_{s}(p)) \rangle$", estimator='mean', ci='sd')
            plt.xlim(math.exp(3.5),math.exp(4.5))
            plt.ylim(0,4.4)
            plt.yticks(np.arange(0, 4.5, step=0.20))
            plt.savefig("Test Fit For C (p--%4.3f).png" %(p_st), dpi=400)
            plt.show()
            plt.close()'''
    
        popt, pcov = curve_fit(ln_pow_law, x1, x2, p0= np.asarray([2, 0.015]))
    
        perr = np.sqrt(np.diag(pcov))
        
        print("SD of C:\t" +str(perr[1]) + " for p:\t" +str(p_st))
        
        C_fits.append([p_st, popt[1], perr[1]])
    
        tukan= (popt[0], popt[1], perr[1])
        
        if (p_st in L):
            #Plot only if in L
            hurtlocker= pan.DataFrame(npframe, columns= ["Trial Number", "p", "s", r"$\langle ln(n_{s}(p)) \rangle$"])
            g= sea.lineplot(data=hurtlocker, x="s" , y=r"$\langle ln(n_{s}(p)) \rangle$", estimator='mean', ci='sd')
            x1=np.sort(x1, kind='mergesort')
            plt.plot(x1, ln_pow_law(x1, *popt), 'm--', label=r'Th Fit: $ \langle ln(n_{s}(p)) \rangle \approx %5.4f - (%5.4f \pm %5.4f) \times s $ ' % tukan )
            plt.legend()
            #plt.xscale('log', basex= math.e)
            plt.xlim(math.exp(2),math.exp(5))
            plt.ylim(-11.5,-7)
            plt.yticks(np.arange(-11.5,-7, step=0.25))
            g.set_title('p = %4.3f' %(p_st))
            plt.savefig("1 Test Fit For C (p--%4.3f).png" %(p_st), dpi=400)
            plt.show()
            plt.close() 
        
        #g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
        #g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    
        
        
        p_st +=0.001
        p_st= round(p_st, 3)
        print(p_st)

    print(C_fits)
    p_c =0.592746
    C_fit=np.array(C_fits)
    C_fit[:,0] -= p_c
    p_pc_abs = np.abs(C_fit[:,0])
    C_fit[:,0] = np.log(p_pc_abs)
    C_fit[:,1] = np.log(C_fit[:,1])
    
    x1= np.transpose(C_fit[:,0])
    x2= np.transpose(C_fit[:,1])
    
    final_locker= pan.DataFrame(C_fit, columns= ["$ln|p - p_{c}|$", r"$\langle ln(C) \rangle$", "error in C"])
    
    '''g= sea.scatterplot(data=final_locker, x="$ln|p - p_{c}|$" , y=r"$\langle ln(C) \rangle$")
    plt.savefig("lnC Scatter (p_st--%4.3f, p_end--%4.3f).png" %(p_st, pend), dpi=400)
    plt.show()
    plt.close()'''
    
    g= sea.scatterplot(data=final_locker, x="$ln|p - p_{c}|$" , y=r"$\langle ln(C) \rangle$")
    popt, pcov = curve_fit(ln_pow_law, x1, x2, p0= np.asarray([-2.75, -0.42]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of lnC:\t" +str(perr[1]))
    
    tukan= (popt[0], -popt[1], perr[1])
    plt.plot(x1, ln_pow_law(x1, *popt), 'm--', label=r'Th Fit: $ \langle ln(C) \rangle = %5.4f + (%5.4f \pm %5.4f) \times ln|p - p_{c}| $ ' % tukan )
    plt.legend()
    g.set_title(r'$ %4.3f \leq p \leq %4.3f $' %(p1, pend))
    plt.savefig("1 lnC Scatter Fit (p_st--%4.3f, p_end--%4.3f).png" %(p1, pend), dpi=400)
    plt.show()
    plt.close()
      
    
    
def crtexpt_finsc():
    # Makes plots based on finite scaling.
    
    os.chdir("..\simulations\CrtExp")
    # Changing to relevant directory.
    
    s = input("The default file to be analysed is file: 'FinSc_SP_p_0.5927_Div_162_G1_25_G2_401.csv'. To modify some other Finite Scaling file, type 'Y' or 'y':\t")
    
    p =0.5927; g1 = 25; g2 = 401; R= 50; Div = 162;
    
    if (s == "Y" or s == "y"):
        g1= float(input("Enter Starting Grid Size:\t"))
        g2= float(input("Enter Ending Grid Size:\t"))
        p= float(input("Enter p (value at which data was collected):\t"))
        R= int(input("Enter number of random trials:\t"))
        Div= int(input("Enter number of divisions:\t"))
    
    data= np.genfromtxt('FinSc_SP_p_%5.4f_Div_%d_G1_%d_G2_%d.csv' %(p, Div, g1, g2), delimiter=",", comments='#', skip_header=1)
    
    p= data[0,0] #Getting exact p value.
    p_c= 0.592746
    
    data[:,0] -= p_c
    
    r=0; GrS_List=[]
    for i in data[:,1]:
        if (r != i):
            r=i
            GrS_List.append(r)
    #Collecting list of all grid sizes available.
    
    npframe = data
    
    hurtlocker= pan.DataFrame(npframe, columns= ["$|p - p_{c}|$", "L", "Trial Number",  "P[p]", r"S[p]"])
    
    
    x1 =np.transpose(npframe[:,1])
    x2= np.transpose(npframe[:,3])
    
    os.chdir(r"..\..\figures")
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("Finite Scaling")==False):
        os.mkdir("Finite Scaling")
    os.chdir("Finite Scaling")
    
    g= sea.lineplot(data=hurtlocker, x="L" , y="P[p]", estimator='mean', ci='sd')
    
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.3, -0.005]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of C:\t" +str(perr[1]) + " for p:\t" +str(p))
    
    tukan= (popt[0], -popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ P[p] = %5.4f \times L^{-(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.legend()
    
    plt.xlim(g1,g2)
    plt.ylim(math.exp(-5), math.exp(0))
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.set_title(r'$p = %8.7f, ( \xi \longrightarrow \infty ) $' %(p))
    plt.savefig("Fig III Beta P(p) vs L (p--%8.7f) Exp.png" %(p), dpi=400)
    plt.show()
    plt.close()
    
    x2= np.transpose(npframe[:,4])
    
    g= sea.lineplot(data=hurtlocker, x="L" , y=r"S[p]", estimator='mean', ci='sd')
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([4, 0.15]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of C:\t" +str(perr[1]) + " for p:\t" +str(p))
    
    tukan= (popt[0], popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ S[p] = %5.4f \times L^{(%5.4f \pm %5.4f)} $ ' % tukan )
    plt.legend()
    plt.xlim(g1,g2)
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.set_title(r'$p = %8.7f, ( \xi \longrightarrow \infty ) $' %(p))
    plt.savefig(" Gamma S(p) vs L (p--%8.7f) Exp.png" %(p), dpi=400)
    plt.show()
    plt.close()
    
    
def crtexpt_nuu():
    # Makes plots based on finite scaling (nuu).
    
    os.chdir("..\simulations\CrtExp")
    # Changing to relevant directory.
    
    p_c=0.592746
    
    s = input("The default file to be analysed is file: 'FinScNu_SP_p_0.5947_Div_123_G1_25_G2_207.csv'. To modify some other Finite Scaling file, type 'Y' or 'y':\t")
    
    p =0.5947; g1 = 25; g2 = 207; R= 100; Div = 123;
    
    if (s == "Y" or s == "y"):
        g1= float(input("Enter Starting Grid Size:\t"))
        g2= float(input("Enter Ending Grid Size:\t"))
        p= float(input("Enter p (value at which data was collected):\t"))
        R= int(input("Enter number of random trials:\t"))
        Div= int(input("Enter number of divisions:\t"))
    
    data= np.genfromtxt('FinScNu_SP_p_%5.4f_Div_%d_G1_%d_G2_%d.csv' %(p, Div, g1, g2), delimiter=",", comments='#', skip_header=1)
                        
    r=0; GrS_List=[]
    for i in data[:,1]:
        if (r != i):
            r=i
            GrS_List.append(r)
    #Collecting list of all grid sizes available.
    
    gingerbread =[] #Stores data in final format
    
    g = GrS_List[0] ; flag =0; a =0; b=0;
    for x in range(0, data[:,1].size):
        
        if (data[x,1] != g and flag == 1 or x == data[:,1].size -1):
            #p slice just ended.
            b = x
            print("For size %f, we have b = %d  and  b - a = %d" %(data[x,1], b, (b-a)))
            mean_p= np.mean(data[a:b,3])
            mean_p2= np.mean(data[a:b,4])
            sd_p = np.std(data[a:b,3])
            gingerbread.append([g, mean_p, mean_p2, sd_p])
            
            g= data[x,1]
            flag =0
        
        if (data[x,1] == g and flag == 0):
            a = x
            flag=1
    
    npframe = np.array(gingerbread)
    print(npframe.shape)
    npframe[:,1] -= p_c
    npframe[:,1] = np.fabs(npframe[:,1])
    
    hurtlocker= pan.DataFrame(npframe, columns= [ "L",  r"$ | \langle p \rangle - p_c | $", r"$ \langle p^{2} \rangle $", r"$ \sigma_{p} $"])
    
    
    x1 =np.transpose(npframe[:,0])
    x2= np.transpose(npframe[:,3])
    
    os.chdir(r"..\..\figures")
    
    if(os.path.isdir("CrtExp")==False):
        os.mkdir("CrtExp")
    os.chdir("CrtExp")
    if(os.path.isdir("Finite Scaling")==False):
        os.mkdir("Finite Scaling")
    os.chdir("Finite Scaling")
    if(os.path.isdir("Nuu")==False):
        os.mkdir("Nuu")
    os.chdir("Nuu")
    
    g= sea.scatterplot(data=hurtlocker, x="L" , y= r"$ \sigma_{p} $")
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.1, -0.75]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of Sigma_p:\t" +str(perr[1]))
    
    tukan= (popt[0], -popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ \sigma_{p} = %5.4f \times L^{-(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.xlim(GrS_List[0]- 2, GrS_List[-1] + 10)
    plt.xlim(math.exp(3), math.exp(5.5))
    plt.ylim(math.exp(-5.25), math.exp(-3))
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    plt.legend()
    g.set_title(r'$ \sigma_{p} \quad vs \quad L$')
    plt.savefig("Nu Sigma_p vs L (G1--%3.0f G2-- %3.0f) Exp.png" %(GrS_List[0], GrS_List[-1]), dpi=400)
    plt.show()
    plt.close()
    
    x2= np.transpose(npframe[:,1])
    
    g= sea.scatterplot(data=hurtlocker, x="L" , y= r"$ | \langle p \rangle - p_c | $")
    popt, pcov = curve_fit(pow_law, x1, x2, p0= np.asarray([0.05, -0.75]))
    
    perr = np.sqrt(np.diag(pcov))
        
    print("SD of p_avg - p_c:\t" +str(perr[1]))
    
    tukan= (popt[0], -popt[1], perr[1])
    plt.plot(x1, pow_law(x1, *popt), 'm--', label=r'Th Fit: $ | \langle p \rangle - p_c | = %5.4f \times L^{-(%5.4f \mp %5.4f)} $ ' % tukan )
    plt.xlim(GrS_List[0]- 2, GrS_List[-1] + 10)
    plt.xlim(math.exp(3), math.exp(5.5))
    plt.yscale('log', basey= math.e)
    plt.xscale('log', basex= math.e)
    g.xaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    g.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
    plt.legend()
    g.set_title(r'$ | \langle p \rangle - p_c | \quad vs \quad L $')
    plt.savefig("Nu p_avg - p_c vs L (G1--%3.0f G2-- %3.0f) Exp.png" %(GrS_List[0], GrS_List[-1]), dpi=400)
    plt.show()
    plt.close()
    

    
plotter()
    
    
            
        
        