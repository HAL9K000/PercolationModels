# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 17:18:39 2021

@author: Koustav
"""

import os
import glob
import matplotlib.pyplot as plt
import seaborn as sea
import numpy as np
import pandas as pan
import math
import collections
import matplotlib.ticker as mtick
from mpl_toolkits import mplot3d
from matplotlib.collections import LineCollection
from scipy.optimize import curve_fit
import powerlaw


def pow_law(x, a, expo):
    return a*(np.power(x, expo))

def trunc_pow_law(x, a, expo, trunc_expo):       #Truncated Power Law
    return a*(np.power(x, expo))*np.exp(trunc_expo*x)


def main_ind():
    
    fandango = np.genfromtxt("PissingAbout15+16.csv", delimiter=",", comments='#', skip_header=1)
    #Stores decay data of cross-correlation between frames as a function of p.
    gaol={} #Stores truncated power law fit data.
    gaol[0.60] =[]; gaol[0.70] =[]; gaol[0.75] =[]; 
    gaol[0.80] =[]; gaol[0.90] =[]; gaol[0.95] =[];
    L=0
    for i in range(6,7):
        base_path = r"22Apret\Apres 256+512\256" + "\\" + str(i)
        files = glob.glob(base_path + "**/**/*.csv", recursive=True)
        for file in files:
            if (file == base_path + r"\dump\15_16_KungF---U.csv"):
                continue
            if (os.path.getsize(file) > 4096):
                #Keeping unwanted files out.
                print(file)
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
                
                '''
                data_temp resembles:
                    | L, p, lag, #, s, s + del(s) |
                    Hai
                '''
                
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                '''if(p == 0.728):
                    print("Skipped")
                    continue'''
                
                data_temp[:,5] -= data_temp[:,4]
                data_temp[:,5] = np.abs(data_temp[:,5])
                temp_freqs = dict(collections.Counter(data_temp[:,5]))
                a,b = data_temp.shape 
                DP_freqs = {k: v / (a) for k, v in temp_freqs.items()}
                
                DP_freqs = np.array(list(DP_freqs.items())) #Converting dictionary to numpy array.
                
                #Sorting array in increasing order of del(s).
                
                #DP_freqs = DP_freqs[DP_freqs[:,0].argsort()]
                
                #Next, to convert PDF into 1 - CDF (P(S >= (DEL(S))))
                
                print("Sorted del(s) PDF:")
                print(DP_freqs)
                
                '''DP_freqs[-2,1] += DP_freqs[-1,1]; #DP_freqs[-1,1] = 0 
                k= len(DP_freqs[:,1]) #Finding total number of del(s) elements
                
                print("Total distinct del(s) samples:\t" +str(k))
                
                for j in range(k-3, -1, -1):
                    #Iterate over the PDF function in reverse.
                    DP_freqs[j,1] += DP_freqs[j+1,1]
                
                print("Sorted del(s) 1-CDF:")
                print(DP_freqs)'''
                
                os.chdir("../../../figures")
            
                if(os.path.isdir("del_S")==False):
                    os.mkdir("del_S")
                os.chdir("del_S")
                if(os.path.isdir("DP")==False):
                    os.mkdir("DP")
                os.chdir("DP")
                if(os.path.isdir("Individual")==False):
                    os.mkdir("Individual")
                os.chdir("Individual")
                '''if(os.path.isdir("1-CDF")==False):
                    os.mkdir("1-CDF")
                os.chdir("1-CDF")'''
                if(os.path.isdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))==False):
                    os.mkdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                os.chdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                
                
                print("p:\t" +str(p) + " L:\t"+ str(L) + " CC:\t" +str(CC))
                #hurtlocker= pan.DataFrame(DP_freqs, columns= [r"$|\Delta s|$", r"$P (S \geq \Delta s)$"])
                hurtlocker= pan.DataFrame(DP_freqs, columns= [r"$|\Delta s|$", r"$P (S = \Delta s)$"])
                fig = plt.figure(figsize=(6.4,4.8))
                f = sea.scatterplot(data=hurtlocker, x=r"$|\Delta s|$" , y=r"$P (S = \Delta s)$")
                f.set_title('p = %f, Grid Size (G) = %d, Cross-Correlation = %3.2f' %(p, L, CC))
                
                
                #Overlaying two seaborn plots.
                #ax = fig.add_subplot(111)
                #sea.scatterplot(data=hurtlocker, x=r"$|\Delta s|$" , y=r"$P (S \geq \Delta s)$", alpha=0.5, s=2, ax= ax)
                #sea.lineplot(data=hurtlocker, x=r"$|\Delta s|$" , y=r"$P (S \geq \Delta s)$", alpha=0.2, ax= ax) #, s=1)
                #ax.set_title('p = %f, Grid Size (G) = %d, Cross-Correlation = %3.2f' %(p, L, CC))
                plt.yscale('log'); plt.xscale('log')
                plt.xlim(1, 10**5)
                plt.ylim(10**(-6.4), 10**(0.1))
                
                plt.savefig("0P(del(s)) vs del(s)  --- p_%f - Grid Size (G)_%d - CC_%3.2f.png" %(p,L,CC), dpi=400)
                #plt.show()
                plt.close()
                
                '''x1 = np.transpose(DP_freqs[:,0])
                x2 = np.transpose(DP_freqs[:,1])
                
                popt, pcov = curve_fit(trunc_pow_law, x1, x2, p0= np.asarray([1, -0.75, -0.0005]), maxfev=5000 )
    
                perr = np.sqrt(np.diag(pcov))
        
                print("SD of exponent:\t" +str(perr[1]) + " for p:\t" +str(p))
    
                tukan= (popt[0], popt[1], perr[1], popt[2], perr[2])
                
                plt.plot(x1, trunc_pow_law(x1, *popt), 'm--', label=r'Fit: $ P (S \geq \Delta s) =  %3.2f \times \Delta s^{(%4.3f \mp %4.3f)}\times e^{(%4.3f \mp %4.3f)\times \Delta s}$ ' % tukan )
                plt.ylim(10**(-6.4), 10**(0.1)); plt.xlim(1, 10**5)
                plt.legend()
                
                plt.savefig("Fit 1- CDF(del(s)) vs del(s)  --- p_%f - Grid Size (G)_%d - CC_%3.2f.png" %(p,L,CC), dpi=400)
                #plt.show()
                plt.close()
                
                #Saving best fit data.
                
                gaol[float(round(CC,2))].append([L, p, -popt[1], perr[1], -popt[2], perr[2]])'''
                
                os.chdir(r"..\..\..\..\..\analysis\Mass Action\DP")
                #break;
    
    #Saving as CSVs.
    '''if(os.path.isdir("del_S")==False):
        os.mkdir("del_S")
    os.chdir("del_S")
    if(os.path.isdir("%d" %(L))==False):
        os.mkdir("%d" %(L))
    os.chdir("%d" %(L))
    K= [0.6, 0.7, 0.75, 0.8, 0.9, 0.95]
    heado = 'L, p, alpha, SD(alpha), lambda, SD(lambda)'
    for k in K:
        np.savetxt("BestFitCDF_CC_%3.2F.csv" %(k), gaol[k], delimiter=',', header=heado, comments='#')
    os.chdir(r"../../")'''
    
    
def main_ccdf_fit():
    
    fandango = np.genfromtxt("PissingAbout15+16.csv", delimiter=",", comments='#', skip_header=1)
    #Stores decay data of cross-correlation between frames as a function of p.
    gaol={} #Stores truncated power law fit data.
    gaol[0.60] =[]; gaol[0.70] =[]; gaol[0.75] =[]; 
    gaol[0.80] =[]; gaol[0.90] =[]; gaol[0.95] =[];
    L=0; crosc= 0.7
    for i in range(0,10):
        base_path = r"22Apret\Apres 256+512\512" + "\\" + str(i)
        files = glob.glob(base_path + "**/**/*.csv", recursive=True)
        for file in files:
            if (file == base_path + r"\dump\15_16_KungF---U.csv"):
                continue
            if (os.path.getsize(file) > 4096):
                #Keeping unwanted files out.
                print(file)
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1, max_rows=3)
                
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                
                if( p == 0.678):
                    print(str(CC) + "  " + str(p) + " shall be skipped.")
                    continue
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
                
                '''
                data_temp resembles:
                    | L, p, lag, #, s, s + del(s) |
                '''
                
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                
                data_temp[:,5] -= data_temp[:,4]
                data_temp[:,5] = np.abs(data_temp[:,5])
                
                fit = powerlaw.Fit(data_temp[:,5],discrete=True,estimate_discrete = False) #If you already know xmin pass it as an argument (xmin=value) for speed
                print("p:\t" +str(p) + " L:\t"+ str(L) + " CC:\t" +str(CC))
                print('x_min: ',fit.xmin)
                print('alpha: ',fit.truncated_power_law.parameter1)
                print('1/lambda: ',1/fit.truncated_power_law.parameter2)
                tukan = (-fit.truncated_power_law.parameter1, -fit.truncated_power_law.parameter2)
                
                fig = fit.plot_ccdf(color ='cornflowerblue', ls='-', linewidth=1.1, alpha=0.2)
                fit.plot_ccdf(color='darkcyan',marker='o', linestyle='', ms=1.2, alpha=0.35, ax=fig)
                #ax = fig.add_subplot(111)
                fit.truncated_power_law.plot_ccdf(color='darkslateblue', linestyle='--', label=r'Fit: $ P (S \geq \Delta s) \propto \Delta s^{(%4.3f)}\times e^{(%6.5f)\times \Delta s}$ ' % tukan, ax=fig)
                fig.set_title('p = %f, Grid Size (G) = %d, Cross-Correlation = %3.2f' %(p, L, CC))
                #x = fit.xmins
                #y = fit.Ds
                #plt.ylim(10**(-6.4), 10**(0.1)); 
                plt.xlim(1, 10**5.3)
                plt.xlabel(r"$|\Delta s|$")
                plt.ylabel(r"$P (S \geq \Delta s)$")
                plt.legend()
                
                os.chdir("../../../figures")
            
                if(os.path.isdir("del_S")==False):
                    os.mkdir("del_S")
                os.chdir("del_S")
                if(os.path.isdir("DP")==False):
                    os.mkdir("DP")
                os.chdir("DP")
                if(os.path.isdir("Individual")==False):
                    os.mkdir("Individual")
                os.chdir("Individual")
                if(os.path.isdir("1-CDF")==False):
                    os.mkdir("1-CDF")
                os.chdir("1-CDF")
                if(os.path.isdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))==False):
                    os.mkdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                os.chdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                
                plt.savefig("Better Fit 1- CDF(del(s)) vs del(s)  --- p_%f - Grid Size (G)_%d - CC_%3.2f.png" %(p,L,CC), dpi=400)
                #plt.show()
                plt.close()
                
                os.chdir("../../")
                if(os.path.isdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))==False):
                    os.mkdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                os.chdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                
                print("Done with CDF Plots And Fits. Moving On To PDF Plots...")
                
                fig = fit.plot_pdf(color='darkcyan',marker='o', linestyle='', ms=1.5, alpha=0.4)
                #fit.plot_pdf(color='darkcyan',marker='o', linestyle='', ms=1.2, alpha=0.35, ax=fig)
                #ax = fig.add_subplot(111)
                fit.truncated_power_law.plot_pdf(color='darkslateblue', linestyle='--', label=r'Fit: $ P (S = \Delta s) \propto \Delta s^{(%4.3f)}\times e^{(%6.5f)\times \Delta s}$ ' % tukan, ax=fig)
                fig.set_title('p = %f, Grid Size (G) = %d, Cross-Correlation = %3.2f' %(p, L, CC))
                #x = fit.xmins
                #y = fit.Ds
                #plt.ylim(10**(-6.4), 10**(0.1)); 
                plt.xlim(1, 10**5.3)
                plt.xlabel(r"$|\Delta s|$")
                plt.ylabel(r"$P (S = \Delta s)$")
                plt.legend()

                plt.savefig("Better Fit PDF(del(s)) vs del(s)  --- p_%f - Grid Size (G)_%d - CC_%3.2f.png" %(p,L,CC), dpi=400)
                #plt.show()
                plt.close()
                
                

                comparison_tpl_exp = fit.distribution_compare('truncated_power_law','exponential',normalized_ratio=True)
                comparison_tpl_streched_exp = fit.distribution_compare('truncated_power_law','stretched_exponential',normalized_ratio=True)
                comparison_tpl_log_normal = fit.distribution_compare('truncated_power_law','lognormal',normalized_ratio=True)
                comparison_tpl_pl = fit.distribution_compare('truncated_power_law','power_law',normalized_ratio=True)
                
                f = open("Taupe.txt", "w+")
                f.write("LR (Power Law): " + str(comparison_tpl_pl[0]) +" p-value: "+ str(comparison_tpl_pl[1]) +"\n")
                f.write("LR (Exponential): " + str(comparison_tpl_exp[0]) +" p-value: "+ str(comparison_tpl_exp[1]) +"\n")
                f.write("LR (Log-Normal): " + str(comparison_tpl_log_normal[0]) +" p-value: "+ str(comparison_tpl_log_normal[1]) +"\n")
                f.write("LR (Stretched-Exponential): " + str(comparison_tpl_streched_exp[0]) +" p-value: "+ str(comparison_tpl_streched_exp[1]) +"\n")
                f.close()
                
                
                print("LR (Power Law): ",comparison_tpl_pl[0]," p-value: ",comparison_tpl_pl[1])
                print("LR (Exponential): ",comparison_tpl_exp[0]," p-value: ",comparison_tpl_exp[1])
                print("LR (Log-Normal): ",comparison_tpl_log_normal[0]," p-value: ",comparison_tpl_log_normal[1])
                print("LR (Stretched-Exponential): ",comparison_tpl_streched_exp[0]," p-value: ",comparison_tpl_streched_exp[1])
                
                gaol[float(round(CC,2))].append([L, p, fit.xmin, fit.truncated_power_law.parameter1, 1/fit.truncated_power_law.parameter2])
                
                os.chdir(r"..\..\..\..\..\analysis\Mass Action\DP")
                
    if(os.path.isdir("del_S")==False):
        os.mkdir("del_S")
    os.chdir("del_S")
    if(os.path.isdir("%d" %(L))==False):
        os.mkdir("%d" %(L))
    os.chdir("%d" %(L))
    K= [0.6, 0.7, 0.75, 0.8, 0.9, 0.95]
    heado = 'L, p, x_min, alpha,  1/lambda'
    for k in K:
        np.savetxt("Nu_Pow_0_6_BestFitCDF_CC_%3.2F.csv" %(k), gaol[k], delimiter=',', header=heado, comments='#')
    os.chdir(r"../../")
                
def main_cumulative():
    
    p_c = 0.725194
    
    crosc = float(input("Enter a Cross-Correlation Value To Be Analysed (Choose Between 0.95, 0.9, 0.8, 0.75, 0.7 & 0.6):\t"))
    
    fandango = np.genfromtxt("PissingAbout15+16.csv", delimiter=",", comments='#', skip_header=1)
    #Stores decay data of cross-correlation between frames as a function of p.
    binder=[]; L=0;
    for i in range(0,10):
        base_path = r"22Apret\Apres 256+512\512" + "\\" + str(i)
        files = glob.glob(base_path + "**/**/*.csv", recursive=True)
        for file in files:
            if (file == base_path + r"\dump\15_16_KungF---U.csv"):
                print('Gandu')
                continue
            if (os.path.getsize(file) > 4096):
                #Keeping unwanted files out.
                print(file)
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1, max_rows=3)
                
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                
                if( CC <= crosc - 0.01 or CC  >= crosc + 0.01):
                    print(str(CC) + " shall be skipped.")
                    continue
                if( p == 0.678):
                    print("Fuck You")
                    continue
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
                
                '''
                data_temp resembles:
                    | L, p, lag, #, s, s + del(s) |
                '''
                
                data_temp[:,5] -= data_temp[:,4]
                data_temp[:,5] = np.abs(data_temp[:,5])
                temp_freqs = dict(collections.Counter(data_temp[:,5]))
                a,b = data_temp.shape 
                DP_freqs = {k: v / a for k, v in temp_freqs.items()}
                
                DP_freqs = np.array(list(DP_freqs.items())) #Converting dictionary to numpy array.
                a,b =DP_freqs.shape
                #col_P= np.zeros((a,1)); col_P = p 
                DP_freqs = np.insert(DP_freqs, 0, p, axis=1)
                
                '''DP_freqs looks like: 
                    | p, del(s), P(del(s))|
                '''
                
                '''DP_freqs = list(DP_freqs.items()) #Converting dictionary to list.
                for j in range(0,len(DP_freqs)):
                    DP_freqs[j].append(p)'''
                    
                print(DP_freqs)
                
                if(len(binder)==0):
                    #First one in the bag.
                    binder = DP_freqs.tolist()
                else:
                    binder.extend(DP_freqs.tolist())
                    
    os.chdir("../../../figures")
    if(os.path.isdir("del_S")==False):
        os.mkdir("del_S")
    os.chdir("del_S")
    if(os.path.isdir("DP")==False):
        os.mkdir("DP")
    os.chdir("DP")
    if(os.path.isdir("3D")==False):
        os.mkdir("3D")
    os.chdir("3D")
    if(os.path.isdir("%d" %(L))==False):
        os.mkdir("%d" %(L))
    os.chdir("%d" %(L))
                
    binder= np.array(binder)
    
    fig=plt.figure()
    ax = plt.axes(projection='3d')
    #surf1 =ax.plot_trisurf(np.log10(binder[:,1]), binder[:,0], np.log10(binder[:,2]), cmap='viridis', edgecolor='none')
    '''for k in range(0,len(self.x1)):
        #Plotting SD bars
        ax.plot([self.x1[k], self.x1[k]], [self.y1[k], self.y1[k]], [self.z1[k] + self.sd_z1[k], self.z1[k] - self.sd_z1[k]], marker="_", markerfacecolor='k', color='k')
    '''
    surf1 =ax.scatter(np.log10(binder[:,1]), binder[:,0], np.log10(binder[:,2]), c= np.log10(binder[:,2]), cmap='viridis', linewidth=0.5)
    cbar1=fig.colorbar(surf1, shrink=0.75)
    cbar1.ax.get_yaxis().labelpad = 12
    cbar1.ax.set_ylabel(r"$P (S=\Delta s)$", rotation=270)
    
    ax.set_xlabel(r"$log_{10}|\Delta s|$")
    ax.set_zlabel(r"$log_{10}|P (S=\Delta s)|$")
    ax.set_ylabel("Occupancy rate (p)")
    #plt.zscale('log'); plt.xscale('log')
                    
    ax.view_init(elev=36.0, azim=-52.0)
    ax.legend()
    ax.set_title(r"$P (S=\Delta s)$ vs $|\Delta s|$, L = %d, $R_{0,0}$ = %3.2f" %(L,crosc))
    plt.savefig("Cumulative Scatter P(del(s)) vs del(s)  --- Grid Size (G)_%d - CC_%3.2f.png" %(L,crosc), dpi=550)
    plt.show()
    plt.close()
    
    
    '''Now for scatter plot'''
    
    fig=plt.figure(figsize=(6.4,4.8))
    #ax = plt.axes(projection='3d')
    ax = fig.add_subplot(111,projection='3d')
    surf1 =ax.scatter(np.log10(binder[:,1]), binder[:,0], np.log10(binder[:,2]), c= np.log10(binder[:,2]), cmap='viridis', linewidth=0.5)
    '''for k in range(0,len(self.x1)):
        #Plotting SD bars
        ax.plot([self.x1[k], self.x1[k]], [self.y1[k], self.y1[k]], [self.z1[k] + self.sd_z1[k], self.z1[k] - self.sd_z1[k]], marker="_", markerfacecolor='k', color='k')
    '''
    cbar1=fig.colorbar(surf1, shrink=0.75)
    cbar1.ax.get_yaxis().labelpad = 12
    cbar1.ax.set_ylabel(r"$log|P (S=\Delta s)|$", rotation=270)
    
    ax.set_xlabel(r"$log_{10}|\Delta s|$")
    ax.set_xlim(-0.1, 5)
    ax.set_zlabel(r"$log_{10}|P (S=\Delta s)|$")
    ax.set_zlim(-6.1, 0)
    ax.set_ylabel("Occupancy rate (p)")
    
    #plt.zscale('log'); plt.xscale('log')
    
    #Plotting p_c plane.
    x = np.linspace(-1,5.5,10)
    z = np.linspace(-7,1,10)
    X,Z = np.meshgrid(x,z)
    Y= 0*X +0*Z + p_c
    
    #ax.hold(True) #Preserve pre-plotted elements.
    
    ax.plot_surface(X,Y,Z, alpha= 0.3, color='k', antialiased=True)
    ax.text(5, p_c, -1, "$p_{c}(q)$", color='0.5')
    '''p_clin = np.array([[0,p_c], [5,p_c]])
    lines = LineCollection([p_clin],zorder=1000,color='0.65',lw=2)
    ax.add_collection3d(lines, zs=-90)'''
    
    ax.view_init(elev=36.0, azim=-52.0)
    ax.legend()
    ax.set_title(r"$log|P (S=\Delta s)|$ vs $log|\Delta s|$, L = %d, $R_{0,0}$ = %3.2f" %(L,crosc))
    plt.savefig("Cumulative Scatter Plane P(del(s)) vs del(s)  --- Grid Size (G)_%d - CC_%3.2f.png" %(L,crosc), dpi=550)
    
    ax.view_init(elev=62.0, azim=-3.0)
    plt.savefig("Cumulative Scatter Plane P(del(s)) vs del(s) Top Down --- Grid Size (G)_%d - CC_%3.2f.png" %(L,crosc), dpi=550)
    plt.show()
    plt.close()
    
    os.chdir(r"..\..\..\..\..\analysis\Mass Action\DP")
    
    
def main_del_s_count():
    
    p_c = 0.725194
    
    crosc = float(input("Enter a Cross-Correlation Value To Be Analysed (Choose Between 0.95, 0.9, 0.8, 0.75, 0.7 & 0.6):\t"))
    
    fandango = np.genfromtxt("PissingAbout15+16.csv", delimiter=",", comments='#', skip_header=1)
    #Stores decay data of cross-correlation between frames as a function of p.
    binder=[]; L=0;
    for i in range(0,10):
        base_path = r"22Apret\Apres 256+512\256" + "\\" + str(i)
        files = glob.glob(base_path + "**/**/*.csv", recursive=True)
        for file in files:
            if (file == base_path + r"\dump\15_16_KungF---U.csv"):
                print('Gandu')
                continue
            if (os.path.getsize(file) > 4096):
                #Keeping unwanted files out.
                print(file)
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1, max_rows=3)
                
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                
                if( CC <= crosc - 0.01 or CC  >= crosc + 0.01):
                    print(str(CC) + " shall be skipped.")
                    continue
                if( p == 0.678):
                    print("Fuck You")
                    continue
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
                
                '''
                data_temp resembles:
                    | L, p, lag, #, s, s + del(s) |
                '''
                
                data_temp[:,5] -= data_temp[:,4]
                data_temp[:,5] = np.abs(data_temp[:,5])
                temp_freqs = dict(collections.Counter(data_temp[:,5]))
                a,b = data_temp.shape 
                DP_freqs = {k: v / a for k, v in temp_freqs.items()}
                
                DP_freqs = np.array(list(DP_freqs.items())) #Converting dictionary to numpy array.
                a,b =DP_freqs.shape
                #col_P= np.zeros((a,1)); col_P = p 
                DP_freqs = np.insert(DP_freqs, 0, p, axis=1)
                
                '''DP_freqs looks like: 
                    | p, del(s), P(del(s))|
                '''
                
                '''DP_freqs = list(DP_freqs.items()) #Converting dictionary to list.
                for j in range(0,len(DP_freqs)):
                    DP_freqs[j].append(p)'''
                    
                print(DP_freqs)
                print("Number of del s counts: " + str(a))
                binder.append([p, a])
                    
    os.chdir("../../../figures")
    if(os.path.isdir("del_S")==False):
        os.mkdir("del_S")
    os.chdir("del_S")
    if(os.path.isdir("DP")==False):
        os.mkdir("DP")
    os.chdir("DP")
    if(os.path.isdir("Bifurcation")==False):
        os.mkdir("Bifurcation")
    os.chdir("Bifurcation")
    if(os.path.isdir("S Count")==False):
        os.mkdir("S Count")
    os.chdir("S Count")
                
    binder= np.array(binder)
    
    hurtlocker= pan.DataFrame(binder, columns= ["p", r"Number of unique $|\Delta s|$ observations"])
    f = sea.scatterplot(data=hurtlocker, x="p" , y=r"Number of unique $|\Delta s|$ observations")#, marker="+")
    #sea.lineplot(data=hurtlocker, x=r"$|\Delta s|$" , y=r"$P (S \geq \Delta s)$", alpha=0.2, ax= ax) #, s=1)
    f.set_title('Unique $|\Delta s|$ observations, Grid Size (G) = %d, Cross-Correlation = %3.2f' %( L, crosc))
    #plt.yscale('log'); #plt.xscale('log')
    #plt.ylim(1, 10**5)
    plt.axvline(x= p_c, color='0.65')
    plt.text(p_c+ 0.003,10**2,r'$p_{c}$',rotation=90, color ='0.65')
    
    plt.savefig("S Count, Grid Size (G) = %d, CC = %3.2f.png" %(L, crosc), dpi=400)
    plt.show()
    plt.close()
    
    os.chdir(r"..\..\..\..\..\analysis\Mass Action\DP")
    
    
def main_del_s_symmetry():
    
    p_mask=[0.658, 0.666, 0.678, 0.689, 0.701, 0.728, 0.739, 0.743, 0.755, 0.773 ]
    
    fandango = np.genfromtxt("PissingAbout15+16.csv", delimiter=",", comments='#', skip_header=1)
    #Stores decay data of cross-correlation between frames as a function of p.
    for i in range(0,10):
        base_path = r"22Apret\Apres 256+512\256" + "\\" + str(i)
        files = glob.glob(base_path + "**/**/*.csv", recursive=True)
        MastBind=[]; L=0
        for file in files:
            if (file == base_path + r"\dump\15_16_KungF---U.csv"):
                continue
            if (os.path.getsize(file) > 4096):
                #Keeping unwanted files out.
                print(file)
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1, max_rows=3)
                
                
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                
                if( p not in p_mask):
                    continue
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
                
                '''
                data_temp resembles:
                    | L, p, lag, #, s, s + del(s) |
                '''
                
                data_temp[:,5] -= data_temp[:,4]
                #data_temp[:,5] = np.abs(data_temp[:,5])
                temp_freqs = dict(collections.Counter(data_temp[:,5]))
                a,b = data_temp.shape 
                DP_freqs = {k: v / (a) for k, v in temp_freqs.items()}
                
                DP_freqs = np.array(list(DP_freqs.items())) #Converting dictionary to numpy array.
                
                #Sorting array in increasing order of del(s).
                
                DP_freqs = DP_freqs[DP_freqs[:,0].argsort()]
                
                #Next, to convert PDF into 1 - CDF (P(S >= (DEL(S))))
                
                print("Sorted del(s) PDF:")
                print(DP_freqs)
                
                #DP_freqs[-2,1] += DP_freqs[-1,1]; #DP_freqs[-1,1] = 0 
                k= len(DP_freqs[:,1]) #Finding total number of del(s) elements
                
                print("Total distinct del(s) samples:\t" +str(k))
                
                '''Performing a log-mod transform
                https://blogs.sas.com/content/iml/2014/07/14/log-transformation-of-pos-neg.html
                https://juluribk.com/dealing-with-plotting-negative-zero-and-positive-values-in-log-scale.html
                '''
                
                DP_freqs[:,0] = np.sign(DP_freqs[:,0])*(np.log10(np.abs(DP_freqs[:,0])+1))
                DP_freqs = np.insert(DP_freqs, 2, float(round(CC,2)), axis=1)
                DP_freqs = np.insert(DP_freqs, 3, p, axis=1)
                
                '''DP_freqs looks like: 
                    |del(s), P(del(s)), CC, p|
                '''
                
                print("Final del(s) PDF:")
                print(DP_freqs)
                
                if(len(MastBind)== 0):
                    #Empty
                    MastBind = DP_freqs
                else:
                    MastBind = np.concatenate((MastBind, DP_freqs), axis=0)
                
                
                '''for j in range(k-3, -1, -1):
                    #Iterate over the PDF function in reverse.
                    DP_freqs[j,1] += DP_freqs[j+1,1]
                
                print("Sorted del(s) 1-CDF:")
                print(DP_freqs)'''
                
                os.chdir("../../../figures")
            
                if(os.path.isdir("del_S")==False):
                    os.mkdir("del_S")
                os.chdir("del_S")
                if(os.path.isdir("DP")==False):
                    os.mkdir("DP")
                os.chdir("DP")
                if(os.path.isdir("Individual")==False):
                    os.mkdir("Individual")
                os.chdir("Individual")
                if(os.path.isdir("Symmetry")==False):
                    os.mkdir("Symmetry")
                os.chdir("Symmetry")
                if(os.path.isdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))==False):
                    os.mkdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                os.chdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                
                
                print("p:\t" +str(p) + " L:\t"+ str(L) + " CC:\t" +str(CC))
                hurtlocker= pan.DataFrame(DP_freqs, columns= [r"$\Delta s$", r"$P (S = \Delta s)$", "Cross-Correlation", "p"])
                fig = plt.figure(figsize=(6.4,4.8))
                
                #Overlaying two seaborn plots.
                #ax = fig.add_subplot(111)
                f= sea.scatterplot(data=hurtlocker, x=r"$\Delta s$" , y=r"$P (S = \Delta s)$")#, alpha=0.5, s=2, ax= ax)
                #sea.lineplot(data=hurtlocker, x=r"$\Delta s$" , y=r"$P (S = \Delta s)$", alpha=0.2, ax= ax) #, s=1)
                f.set_title('p = %f, Grid Size (G) = %d, Cross-Correlation = %3.2f' %(p, L, CC))
                plt.yscale('log'); #plt.xscale('log')
                #plt.xlim(1, 10**5)
                plt.ylim(10**(-6.4), 10**(0.1))
                #plt.xlim(-5, 5)
                
                '''x1 = np.transpose(DP_freqs[:,0])
                x2 = np.transpose(DP_freqs[:,1])
                
                popt, pcov = curve_fit(trunc_pow_law, x1, x2, p0= np.asarray([1, -0.75, -0.0005]), maxfev=5000 )
    
                perr = np.sqrt(np.diag(pcov))
        
                print("SD of exponent:\t" +str(perr[1]) + " for p:\t" +str(p))
    
                tukan= (popt[0], popt[1], perr[1], popt[2], perr[2])
                plt.plot(x1, trunc_pow_law(x1, *popt), 'm--', label=r'Fit: $ P (S \geq \Delta s) =  %3.2f \times \Delta s^{(%4.3f \mp %4.3f)}\times e^{(%4.3f \mp %4.3f)\times \Delta s}$ ' % tukan )
                plt.legend()'''
                
                plt.savefig("Symmetry PDF(del(s)) vs del(s)  --- p_%f - Grid Size (G)_%d - CC_%3.2f.png" %(p,L,CC), dpi=400)
                #plt.show()
                plt.close()
                os.chdir(r"..\..\..\..\..\..\analysis\Mass Action\DP")
                #break;
        #Plotting cumulative results.
        
        os.chdir("../../../figures")
        
        if(os.path.isdir("del_S")==False):
            os.mkdir("del_S")
        os.chdir("del_S")
        if(os.path.isdir("DP")==False):
            os.mkdir("DP")
        os.chdir("DP")
        if(os.path.isdir("Individual")==False):
            os.mkdir("Individual")
        os.chdir("Individual")
        if(os.path.isdir("Symmetry")==False):
            os.mkdir("Symmetry")
        os.chdir("Symmetry")
        if(os.path.isdir("Cum")==False):
            os.mkdir("Cum")
        os.chdir("Cum")
        
        hurtlocker= pan.DataFrame(MastBind, columns= [r"$\Delta s$", r"$P (S = \Delta s)$", "Cross-Correlation", "p"])
        fig = plt.figure(figsize=(6.4,4.8))
                
        #Overlaying two seaborn plots.
        #ax = fig.add_subplot(111)
        f= sea.scatterplot(data=hurtlocker, x=r"$\Delta s$" , y=r"$P (S = \Delta s)$", hue="Cross-Correlation")#, alpha=0.5, s=2, ax= ax)
        #sea.lineplot(data=hurtlocker, x=r"$\Delta s$" , y=r"$P (S = \Delta s)$", alpha=0.2, ax= ax) #, s=1)
        f.set_title('p = %f, Grid Size (G) = %d' %(MastBind[0,3], L))
        plt.yscale('log'); #plt.xscale('log')
        #plt.xlim(1, 10**5)
        plt.ylim(10**(-6.4), 10**(0.1))
        plt.xlim(-5, 5)
        
        plt.savefig("Alt Cum Symmetry PDF(del(s)) vs del(s)  --- p_%f - Grid Size (G)_%d.png" %(MastBind[0,3], L), dpi=400)
        #plt.show()
        plt.close()
        os.chdir(r"..\..\..\..\..\..\analysis\Mass Action\DP")
        
    
    
def main_bifurcation():
    
    p_c = 0.725194
    
    crosc = float(input("Enter a Cross-Correlation Value To Be Analysed (Choose Between 0.95, 0.9, 0.8, 0.75, 0.7 & 0.6):\t"))
    #crosc =0.8
    fandango = np.genfromtxt("PissingAbout15+16.csv", delimiter=",", comments='#', skip_header=1)
    #Stores decay data of cross-correlation between frames as a function of p.
    binder=[]; L=0;
    for i in range(0,10):
        base_path = r"22Apret\Apres 256+512\256" + "\\" + str(i)
        files = glob.glob(base_path + "**/**/*.csv", recursive=True)
        for file in files:
            if (file == base_path + r"\dump\15_16_KungF---U.csv"):
                print('Gandu')
                continue
            if (os.path.getsize(file) > 4096):
                #Keeping unwanted files out.
                print(file)
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1, max_rows=3)
                
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                
                if( CC <= crosc - 0.01 or CC  >= crosc + 0.01):
                    print(str(CC) + " shall be skipped.")
                    continue
                if( p == 0.678):
                    print("Fuck You")
                    continue
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
                
                '''
                data_temp resembles:
                    | L, p, lag, #, s, s + del(s) |
                '''
                
                data_temp[:,5] -= data_temp[:,4]
                data_temp[:,5] = np.abs(data_temp[:,5])
                temp_freqs = dict(collections.Counter(data_temp[:,5]))
                a,b = data_temp.shape 
                DP_freqs = {k: v / a for k, v in temp_freqs.items()}
                
                DP_freqs = np.array(list(DP_freqs.items())) #Converting dictionary to numpy array.
                a,b =DP_freqs.shape
                
                split_data = DP_freqs[:,1] < 10**(-5.6)
                DP_freqs = DP_freqs[split_data]
                print("Half Done:")
                print(DP_freqs)
                split_data = DP_freqs[:,1] > 10**(-6)
                DP_freqs_band = DP_freqs[split_data]   #Stores the band of del(s) values whose probability lie between 10^(-5.85) and 10^(-5.85)
                
                #col_P= np.zeros((a,1)); col_P = p 
                DP_freqs_band = np.insert(DP_freqs_band, 0, p, axis=1)
                
                DP_freqs_band = DP_freqs_band[DP_freqs_band[:,1].argsort()]
                #Sorting in increasing values of del(s)
                
                print("Total number of points in given gap for p:\t"+str(p) +" is: \t" +str(len(DP_freqs_band[:,2])) +"\n")
                print(DP_freqs_band)
                
                '''DP_freqs looks like: 
                    | p, del(s), P(del(s))|
                '''
                flag=0
                for j in range(1, len(DP_freqs_band[:,2])-1):
                    if(abs(DP_freqs_band[j,1] -DP_freqs_band[j-1,2]) > 411 or abs(DP_freqs_band[j,1] -DP_freqs_band[j+1,2]) > 411):
                        # 10^(3.3) - 10^(3.2) = 410.369
                        binder.append([p,DP_freqs_band[j,1]])
                        flag=1
                
                if(flag==0):
                    #No del(s) value satisfied the bandwidth demand.
                    #if()
                    binder.append([p,DP_freqs_band[-1,1]])  
                    #Append the very last value
                    
    os.chdir("../../../figures")
    if(os.path.isdir("del_S")==False):
        os.mkdir("del_S")
    os.chdir("del_S")
    if(os.path.isdir("DP")==False):
        os.mkdir("DP")
    os.chdir("DP")
    if(os.path.isdir("Bifurcation")==False):
        os.mkdir("Bifurcation")
    os.chdir("Bifurcation")
    if(os.path.isdir("%d" %(L))==False):
        os.mkdir("%d" %(L))
    os.chdir("%d" %(L))
                
    binder= np.array(binder)
    
    hurtlocker= pan.DataFrame(binder, columns= ["p", r"$|\Delta s|$ s.t. $P (\Delta s \geq 10^{-6})$"])
    f = sea.scatterplot(data=hurtlocker, x="p" , y=r"$|\Delta s|$ s.t. $P (\Delta s \geq 10^{-6})$", marker="+")
    #sea.lineplot(data=hurtlocker, x=r"$|\Delta s|$" , y=r"$P (S \geq \Delta s)$", alpha=0.2, ax= ax) #, s=1)
    f.set_title('Bifurcation Map, Grid Size (G) = %d, Cross-Correlation = %3.2f' %( L, crosc))
    plt.yscale('log'); #plt.xscale('log')
    plt.ylim(1, 10**5)
    plt.axvline(x= p_c, color='0.65')
    plt.text(p_c+ 0.003,10**1,r'$p_{c}$',rotation=90, color ='0.65')
    
    plt.savefig("Bifurcation Map, Grid Size (G) = %d, CC = %3.2f.png" %(L, crosc), dpi=400)
    plt.show()
    plt.close()
    
    os.chdir(r"..\..\..\..\..\analysis\Mass Action\DP")
    
    
    
def plot_fit_pdf():
    
    twist =(-1.2912647288993737, -(1/37.72480211483688))
    fandango = np.genfromtxt("PissingAbout15+16.csv", delimiter=",", comments='#', skip_header=1)
    #Stores decay data of cross-correlation between frames as a function of p.
    gaol={} #Stores truncated power law fit data.
    gaol[0.60] =[]; gaol[0.70] =[]; gaol[0.75] =[]; 
    gaol[0.80] =[]; gaol[0.90] =[]; gaol[0.95] =[];
    L=0; crosc= 0.7
    for i in range(0,1):
        base_path = r"22Apret\Apres 256+512\256" + "\\" + str(i)
        files = glob.glob(base_path + "**/**/*.csv", recursive=True)
        for file in files:
            if (file == base_path + r"\dump\15_16_KungF---U.csv"):
                continue
            if (os.path.getsize(file) > 4096):
                #Keeping unwanted files out.
                print(file)
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1, max_rows=3)
                
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                
                if( CC <= crosc - 0.01 or CC  >= crosc + 0.01 or p != 0.66):
                    print(str(CC) + "  " + str(p) + " shall be skipped.")
                    continue
                
                data_temp= np.genfromtxt(file, delimiter=",", comments='#', skip_header=1)
                
                '''
                data_temp resembles:
                    | L, p, lag, #, s, s + del(s) |
                '''
                '''
                p= data_temp[0,1]; L= int(data_temp[0,0]); CC= cross_cor(fandango, data_temp[0,2], L, p)
                
                data_temp[:,5] -= data_temp[:,4]
                data_temp[:,5] = np.abs(data_temp[:,5])
                
                temp_freqs = dict(collections.Counter(data_temp[:,5]))
                a,b = data_temp.shape 
                DP_freqs = {k: v / (a) for k, v in temp_freqs.items()}
                
                DP_freqs = np.array(list(DP_freqs.items())) #Converting dictionary to numpy array.
                
                #Sorting array in increasing order of del(s).
                
                DP_freqs = DP_freqs[DP_freqs[:,0].argsort()]
                
                
                
                print("Sorted del(s) PDF:")
                print(DP_freqs)
                
                os.chdir("../../../figures")
            
                if(os.path.isdir("del_S")==False):
                    os.mkdir("del_S")
                os.chdir("del_S")
                if(os.path.isdir("DP")==False):
                    os.mkdir("DP")
                os.chdir("DP")
                if(os.path.isdir("Individual")==False):
                    os.mkdir("Individual")
                os.chdir("Individual")
                if(os.path.isdir("1-CDF")==False):
                    os.mkdir("1-CDF")
                os.chdir("1-CDF")
                if(os.path.isdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))==False):
                    os.mkdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                os.chdir("L_%d_p_%4.3f" %(int(data_temp[0,0]), data_temp[0,1]))
                
                print("p:\t" +str(p) + " L:\t"+ str(L) + " CC:\t" +str(CC))
                hurtlocker= pan.DataFrame(DP_freqs, columns= [r"$|\Delta s|$", r"$P (S = \Delta s)$"])
                fig = plt.figure(figsize=(6.4,4.8))
                
                #Overlaying two seaborn plots.
                ax = fig.add_subplot(111)
                sea.scatterplot(data=hurtlocker, x=r"$|\Delta s|$" , y=r"$P (S = \Delta s)$", ax= ax)#, alpha=0.5, s=2, ax= ax)
                #sea.lineplot(data=hurtlocker, x=r"$|\Delta s|$" , y=r"$P (S = \Delta s)$", alpha=0.2, ax= ax) #, s=1)
                ax.set_title('p = %f, Grid Size (G) = %d, Cross-Correlation = %3.2f' %(p, L, CC))
                plt.yscale('log'); plt.xscale('log')
                plt.xlim(1, 10**5)
                plt.ylim(10**(-6.3), 10**(0.1))
                
                x1 = np.transpose(DP_freqs[:,0])
                x2 = np.transpose(DP_freqs[:,1])
                
                #popt, pcov = curve_fit(trunc_pow_law, x1, x2, p0= np.asarray([1, -0.75, -0.0005]), maxfev=5000 )
    
                #perr = np.sqrt(np.diag(pcov))
        
                #print("SD of exponent:\t" +str(perr[1]) + " for p:\t" +str(p))
    
                #tukan= (popt[0], popt[1], perr[1], popt[2], perr[2])
                
                plt.plot(x1, trunc_pow_law(x1, *twist), color='darkslateblue', linestyle='--', label=r'Fit: $ P (S = \Delta s) =  %3.2f \times \Delta s^{(%4.3f)}\times e^{(%6.5f)\times \Delta s}$ ' % tukan )
                plt.ylim(10**(-6.4), 10**(0.1)); plt.xlim(1, 10**5)
                plt.legend()
                
                plt.savefig("Fit 1- CDF(del(s)) vs del(s)  --- p_%f - Grid Size (G)_%d - CC_%3.2f.png" %(p,L,CC), dpi=400)
                #plt.show()
                plt.close()                
                
                
                
                #Next, to convert PDF into 1 - CDF (P(S >= (DEL(S))))
                
                DP_freqs[-2,1] += DP_freqs[-1,1]; #DP_freqs[-1,1] = 0 
                k= len(DP_freqs[:,1]) #Finding total number of del(s) elements
                
                print("Total distinct del(s) samples:\t" +str(k))
                
                for j in range(k-3, -1, -1):
                    #Iterate over the PDF function in reverse.
                    DP_freqs[j,1] += DP_freqs[j+1,1]
                
                print("Sorted del(s) 1-CDF:")
                print(DP_freqs)
                
                
                plt.savefig("Even Better Fit 1- CDF(del(s)) vs del(s)  --- p_%f - Grid Size (G)_%d - CC_%3.2f.png" %(p,L,CC), dpi=400)
                #plt.show()
                plt.close()
                
                

                comparison_tpl_exp = fit.distribution_compare('truncated_power_law','exponential',normalized_ratio=True)
                comparison_tpl_streched_exp = fit.distribution_compare('truncated_power_law','stretched_exponential',normalized_ratio=True)
                comparison_tpl_log_normal = fit.distribution_compare('truncated_power_law','lognormal',normalized_ratio=True)
                comparison_tpl_pl = fit.distribution_compare('truncated_power_law','power_law',normalized_ratio=True)
                
                print("LR (Power Law): ",comparison_tpl_pl[0]," p-value: ",comparison_tpl_pl[1])
                print("LR (Exponential): ",comparison_tpl_exp[0]," p-value: ",comparison_tpl_exp[1])
                print("LR (Log-Normal): ",comparison_tpl_log_normal[0]," p-value: ",comparison_tpl_log_normal[1])
                print("LR (Stretched-Exponential): ",comparison_tpl_streched_exp[0]," p-value: ",comparison_tpl_streched_exp[1])
                
                gaol[float(round(CC,2))].append([L, p, fit.xmin, fit.truncated_power_law.parameter1, 1/fit.truncated_power_law.parameter2])
                
                os.chdir(r"..\..\..\..\..\..\analysis\Mass Action\DP")
                '''

                
                
                
def cross_cor(grim_fandango, lag, L, p):
    CC=0; k= 128/L
    for t in  range(0, len(grim_fandango[:,0])):
        if grim_fandango[t,0] == p:
            CC = grim_fandango[t,1]+ grim_fandango[t,3]*(math.exp(lag*grim_fandango[t,5]*k*k)); break;
    
    #Calculating cross-correlation b/w frames.
    
    print("CC:\t"+ str(CC))
    return CC;
                
main_ind()            