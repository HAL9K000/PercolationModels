# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 01:44:35 2021

@author: Koustav
"""
import matplotlib.pyplot as plt
import seaborn as sea
import numpy as np
import os

def main():
    #Plot for the p_c(q) thresholds for TCP models.
    
    p_c= [0.7284, 0.725484212, 0.72450002, 0.723381678, 0.722279867, 0.720339457, 
          0.71917661, 0.718261719, 0.716491699, 0.715192159, 0.71309789]
    
    q= [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    
    p_c_sd = [0.0024, 0.001929008, 0.001599747, 0.002012191, 0.001688647, 0.002135829,
              0.001709018, 0.001660967, 0.001987551, 0.00216563, 0.001992864]
    
    bet_nu= [0.264, 0.2941, 0.35, 0.3154, 0.4026, 0.4663, 0.4783, 0.5706, 0.596, 0.6411, 0.7982]
    bet_nu_sd = [0.0111, 0.0287, 0.0289, 0.0282, 0.0291, 0.031, 0.0355, 0.0316, 0.0333, 0.0338, 0.0362]
    
    gam_nu=[1.65522, 1.7507, 1.7946, 1.8484, 1.7683, 1.7376, 1.80 ,1.7809, 1.7207, 1.6845, 1.6616]
    gam_nu_sd = [0.0977, 0.0421, 0.0429, 0.0438, 0.0415, 0.0392, 0.04, 0.0373, 0.0356, 0.0356, 0.0336]
    
    low_bound= list(np.array(p_c) - np.array(p_c_sd))
    upp_bound= list(np.array(p_c) + np.array(p_c_sd))
    
    fig, ax = plt.subplots(figsize=(6,4))
    plt.plot(p_c, q, label="$p_{c}(q)$", marker="s")
    ax.plot(low_bound, q, color='tab:blue', alpha=0.1)
    ax.plot(upp_bound, q, color='tab:blue', alpha=0.1)
    ax.fill_betweenx(q, low_bound, upp_bound, alpha=0.2)
    ax.set_xlabel("$p_{c}(q)$ for TCP Process")
    ax.set_ylabel('q for TCP Process')
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    plt.xlim(0.6, 0.8)
    plt.ylim(0, 0.15)
    ax.set_title("$p_c(q)$ Percolation Thresholds For TCP Process")
    plt.legend()
    
    os.chdir("../../../figures/CrtExp/TCP")
    plt.savefig("Alt Percolation Threshold Phase Space TCP.png", dpi=400, transparent=True)
    plt.show()
    plt.close()
    
    
    #Now plotting critical exponents as a function of q.
    
    os.chdir("Finite Scaling")
    
    low_bound= list(np.array(bet_nu) - np.array(bet_nu_sd))
    upp_bound= list(np.array(bet_nu) + np.array(bet_nu_sd))
    
    fig, ax = plt.subplots(figsize=(6,4))
    plt.plot(q, bet_nu, label=r"$(\beta/\nu)(q)$", color='tab:green', marker="s")
    ax.plot(q, low_bound, color='tab:green', alpha=0.1)
    ax.plot(q, upp_bound, color='tab:green', alpha=0.1)
    ax.fill_between(q, low_bound, upp_bound, color='tab:green', alpha=0.2)
    ax.set_ylabel(r"$(\beta/\nu)(q)$ for TCP Process")
    ax.set_xlabel('q for TCP Process')
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #plt.xlim(0.1, 0.9)
    plt.xlim(0, 0.12)
    ax.set_title(r"Finite Scaling Critical Exponent: $(\beta/\nu)(q)$ For TCP Process")
    plt.legend()
    
    plt.savefig("Alt Finite Scaling Critical Exponent Beta vs Q TCP.png", dpi=400, transparent=True)
    plt.show()
    plt.close()
    
    low_bound= list(np.array(gam_nu) - np.array(gam_nu_sd))
    upp_bound= list(np.array(gam_nu) + np.array(gam_nu_sd))
    
    fig, ax = plt.subplots(figsize=(6,4))
    plt.plot(q, gam_nu, label=r"$(\gamma/\nu)(q)$", color='tab:olive', marker="s")
    ax.plot(q, low_bound, color='tab:olive', alpha=0.1)
    ax.plot(q, upp_bound, color='tab:olive', alpha=0.1)
    ax.fill_between(q, low_bound, upp_bound, color='tab:olive', alpha=0.2)
    ax.set_ylabel(r"$(\gamma/\nu)(q)$ for TCP Process")
    ax.set_xlabel('q for TCP Process')
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    plt.ylim(1.45, 2.05)
    plt.xlim(0, 0.12)
    ax.set_title(r"Finite Scaling Critical Exponent: $(\gamma/\nu)(q)$ For TCP Process")
    
    plt.plot(0.06, 1.81, color= "red", marker=(6,2,0)) #Red asterisk
    plt.plot(0.06, 1.9912, color='tab:olive', marker="s", alpha=0.55, label=r"Alt estimate for $(\gamma/\nu)$ at $q=0.06$")
    plt.legend()
    
    plt.savefig("Alt Finite Scaling Critical Exponent Gamma vs Q TCP.png", dpi=400, transparent=True)
    plt.show()
    plt.close()
    
    
    os.chdir("../../../../analysis/Mass Action/TCP")
    
    
main()
    
