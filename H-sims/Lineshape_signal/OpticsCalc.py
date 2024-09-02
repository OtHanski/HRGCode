from math import *
import numpy as np
import matplotlib.pyplot as plt

def calcW0(W = 1E-3,lam = 243E-9, z = 3.5):
    if False:
        print(repr(W)+":\n"+repr(W**4-4*(lam/pi)**2*z**2))
    W0 = sqrt((W**2-sqrt(W**4-4*(lam/pi)**2*z**2))/2)
    return W0
    
def calcz0(W0 = 300E-6, lam = 243E-9):
    z0 = pi*W0**2/lam
    return z0
    
def calcWz(z, W0 = 300E-6, lam = 243E-9):
    # W as function of Z, given beam parameters
    z0 = (pi/lam)*W0**2
    W = W0*sqrt(1+(z/z0)**2)
    return W

def simW0(Wmin = 1E-3, Wmax = 3E-3, lam = 243E-9, z = 3.5, samples = 10000, drw = True):
    Warr = np.linspace(Wmin,Wmax,samples)
    W0 = np.zeros(len(Warr))
    for i in range(len(Warr)):
        W0[i] = calcW0(W = Warr[i])
    
    fig,ax = plt.subplots()
    
    # Set axes
    ax.set_xlabel("Beam width RT, $W$ [$\mathrm{\mu m}$]")
    ax.set_ylabel("Theoretical minimum beam waist $W_0$ [$\mathrm{\mu m}$]")
    #secax = ax.secondary_yaxis("right", functions = (lambda x: (Etrans*10**9)*x, lambda x: x/(Etrans*10**9)))
    #secax.set_ylabel("$L_α$ fluorescence power [nW]")
    
    ax.plot(Warr*1E6,W0*1E6)
    logsc = False
    
    if drw:
        if logsc:
            plt.yscale("log")
            plt.grid(which = "both", linestyle = "--")
        else:
            plt.grid(linestyle = "--")
        plt.legend()
        plt.show()