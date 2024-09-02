### FOR calculating linewidth contributions ###
###     and finding such parameters that    ###
###      allow for minimizing the width     ###
from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


naive = True
samples = 10000
### Nat constants ###
# Fine structure constant, -
alph = 0.0072973525628
# Planck constant, eV/Hz
h = 4.135667696E-15
# Bohr magneton, eV/G
muB = 5.7883818060E-9
# Boltzmann constant, eV/K
kB = 8.617333262E-5
kBm = 1.380649E-23
# Gravitational accel, m/s^2
g = 9.8067
# Lightspeed, m/s^2
c = 299792458

# H constants
# PI crosssection, m^2 (Killian)
PIC = 7.9E-22
# Mass, kg
m_H = MH = 1.6735575E-27
#####################

### Instrument variables ###
# Laser power, W
P = 10E-3
# Beam waist, m
w0 = 250E-6
# Intensity at waist
I = P/(pi*w0**2/4)

# Beam wavelength, m
lam = 240E-9
f = c/lam
# Propagation constant, 1/m
k = 2*pi/lam

# Confocal parameter, m
b = k*w0**2
# Beam misalignment, rad
dTheta = 10E-4

# Trap parameters
# Starting density, 1/m^3
n = 1E18
# Temperature, K
T = 1E-4
v0 = sqrt(2*kB*T/m_H)
# Sample radius, m
r = 1E-3
# Sample length, m
l = kB*T/(m_H*g)
# Stray fields, V/m
Es = 1
# Wall magnetic field, G
Bw = 15
# Bias field, G
B0 = 3
# Field/temp coefficient, G/K
BT = 1.5E4
# Trap radius, m
rT = 8E-3
############################

def nat_brd():
    return 8.23/(2*pi)

def Zeeman_brd(Bw = Bw, w0 = w0, poles = 12, worst_case = True):
    # Zeeman shift/B, Hz/G (multiply w/B)
    sft = (alph**2*muB)/(4*h)
    
    # Temperature as defined by wall field
    T = 0.2*Bw/BT
    # B0 should be of the order of sample temp
    B0 = T*BT/2
    # Bw = B0 + c_12*rT^5, 12-pol
    if poles == 12:
        c = (Bw-B0)/rT**5
        #r_s covering region where field between 0.1 and 0.2 Bw
        rmax = ((0.2*Bw-B0)/c)**(1/5)
        #print(rmax)
    # Bw = B0 + c_8*rT^3, 8-pol
    if poles == 8:
        c = (Bw-B0)/rT**3
        #r_s covering region where field between 0.1 and 0.2 Bw
        rmax = ((0.2*Bw-B0)/c)**(1/3)
    # Bw = B0 + c_4*rT, 4-pol
    if poles == 4:
        c = (Bw-B0)/rT
        #r_s covering region where field between 0.1 and 0.2 Bw
        rmax = ((0.2*Bw-B0)/c)
        #print(rmax)
    
    def B(r):
        if poles == 12:
            return B0+c*r**5
        if poles == 8:
            return B0+abs(c*r**3)
        if poles == 4:
            return B0+abs(c*r)
    
    #r_s covering region where field between 0.1 and 0.2 Bw
    if worst_case:
        dB = B(rmax)-B0
    else:
        dB = B(w0)-B0
    
    brd = sft*dB
    if isinstance(brd,float):
        if brd < 0:
            brd = np.nan
    elif isinstance(brd, np.ndarray):
        for i in range(len(brd)):
            if brd[i] < 0: brd[i] = np.nan
    
    return brd

def PI_brd(I = I):
    if naive:
        # From thesis of Killian
        return 2*9.7E-4*I/(2*pi)
    else:
        pass
        # TODO: Replace simple addition/multiplication with integration over beam area w/ gaussian beam profile
        return 2*9.7E-4*I/(2*pi)
        
def ACStark_brd(I = I):
    if naive:
        # Killian again
        return 3.34E-4*I

def Stark_brd(Es = Es):
    # Sandberg, E in V/m
    return 4.4E6*Es**2

def CollQuench_brd(n = n, T = T):
    # From thesis of Sandberg, extrapol high T case. At low T should be less
    # n in 1/m^3
    return 5E-16*n*T**(1/3)/(2*pi)

def CCFS_brd(n = n, T = T):
    # Cold collision frequency shift, Killian, n in 1/cm3
    return 7E-9*n*(T)**(1/2)

def Doppler1_brd(dTheta = dTheta, T = T):
    # Sandberg
    return 8.8E8*dTheta*(T)**(1/2)

def Doppler2_brd(T = T):
    # Sandberg
    return 341*T

def TOF_brd(T = T, w0 = w0):
    FWHM = log(2)/pi*(2*kBm*B*T/m_H)**(1/2)/w0
    return FWHM
    

def SampleLossRate():
    pass

def I_limits():
    # Test Photoionization, ACStark
    w0 = np.linspace(220E-6, 1000E-6, samples)
    # Intensity at waist
    I = P/(pi*w0**2/4)
    
    PI = PI_brd(I)
    ACS = ACStark_brd(I)
    PIACS = PI+ACS
    
    fig, ax = plt.subplots()
    ax.plot(w0*1E6,PI, label = "Photoionization")
    ax.plot(w0*1E6,ACS, label = "AC Stark effect")
    ax.plot(w0*1E6,PIACS, label = "Combined $I_{243}$ widening")
    
    ax.set_xlabel("$w_0$ [$\mathrm{\mu m}$]")
    ax.set_ylabel("$\Delta\omega_{1S-2S}$ contribution [Hz]")
    secax = ax.secondary_xaxis("top", functions = (lambda x: P/(pi*(x/2E4)**2), lambda x: (P/(pi*x))**(1/2)*2E4))
    secax.set_xlabel("243 nm beam intensity [$\mathrm{W/cm^2}$]")
    
    ax.ticklabel_format(style = "sci", scilimits = (-2,5))
    ax.set_yticks([nat_brd(), 10, 25, 50, 75, 100, 125, 150, 175])
    secax.set_xticks([25,15,10,5,3,2])
    plt.grid(linestyle = "--")
    
    plt.legend()
    plt.show()

def E_limits():
    # Test regular Stark effect, E = V/m
    E = np.linspace(0, 0.005, samples)
    # Intensity at waist
    
    brd =  Stark_brd(E)
    
    fig, ax = plt.subplots()
    ax.plot(E,brd, label = "Stark effect")
    
    ax.set_xlabel("$E$ [$\mathrm{V/m}$]")
    ax.set_ylabel("$\Delta\omega_{1S-2S}$ contribution [Hz]")
    
    ax.ticklabel_format(style = "sci",axis = "x", scilimits = (0,0))
    ax.set_yticks([nat_brd(), 10, 20, 40, 60, 80, 100, 120])
    plt.grid(linestyle = "--")
    
    plt.legend()
    plt.show()

def Dopp_limits():
    # Test Doppler effects
    dThet = np.geomspace(1E-8, 1E-5, samples)
    T = np.linspace(0, 1.05E-3,samples)
    
    # Calculate 2D space of D1 broadening
    Dopp1_brd = np.multiply.outer(T**(1/2), 8.8E8*dThet)
    # Same for D2
    Dopp2_brd = np.multiply.outer(Doppler2_brd(T), dThet**0)
    # Sum
    Dopptot_brd = Dopp1_brd + Dopp2_brd
    
    #Plot levels
    lvs = [nat_brd(), 10, 100, 1000]
    cmap = mpl.cm.plasma
    
    fig1, ax1 = plt.subplots()
    contot = ax1.contourf(dThet,T,Dopptot_brd, levels = lvs, cmap = mpl.cm.plasma, \
                          norm = mpl.colors.BoundaryNorm(lvs, cmap.N, extend = "both"),\
                          extend = "both")
    
    cbartot = fig1.colorbar(contot)
    cbartot.ax.set_ylabel("Linewidth broadening [Hz]")
    ax1.set_title("Sum Doppler broadening")
    ax1.set_xlabel("$\Delta \Theta$ [rad]")
    ax1.set_ylabel("$T$ [K]")
    
    ax1.ticklabel_format(style = "sci",axis = "both", scilimits = (0,0))
    plt.grid(linestyle = "--")
    
    fig2, ax2 = plt.subplots()
    contot = ax2.contourf(dThet,T,Dopp1_brd, levels = lvs, cmap = mpl.cm.plasma, \
                          norm = mpl.colors.BoundaryNorm(lvs, cmap.N, extend = "both"),\
                          extend = "both")
    
    cbartot = fig2.colorbar(contot)
    cbartot.ax.set_ylabel("Linewidth broadening [Hz]")
    ax2.set_title("First order Doppler broadening")
    ax2.set_xlabel("$\Delta \Theta$ [rad]")
    ax2.set_ylabel("$T$ [K]")
    
    ax2.ticklabel_format(style = "sci",axis = "both", scilimits = (0,0))
    plt.grid(linestyle = "--")
    plt.xscale("log")
    
    fig3, ax3 = plt.subplots()
    contot = ax3.contourf(dThet,T,Dopp2_brd, levels = lvs, cmap = mpl.cm.plasma, \
                          norm = mpl.colors.BoundaryNorm(lvs, cmap.N, extend = "both"),\
                          extend = "both")
    
    cbartot = fig3.colorbar(contot)
    cbartot.ax.set_ylabel("Linewidth broadening [Hz]")
    ax3.set_title("Second order Doppler broadening")
    ax3.set_xlabel("$\Delta \Theta$ [rad]")
    ax3.set_ylabel("$T$ [K]")
    
    ax3.ticklabel_format(style = "sci",axis = "both", scilimits = (0,0))
    plt.grid(linestyle = "--")
    
    plt.show()

def Coll_limits():
    # Test collision broadening
    n = np.geomspace(1E10, 1E13, samples)
    T = np.linspace(0, 1.05E-3, samples)
    X, Y = np.meshgrid(n, T)
    
    # Calculate 2D space of Collision quenching
    CQ_brd = np.multiply.outer(T**(1/3),5E-10/(2*pi)*n)
    # Same for cold collision frequency shift 7E-9*n*(T)**(1/2)
    CCFS_brd = np.multiply.outer(T**(1/2),7E-9*n)
    # Sum
    Colltot_brd = CQ_brd + CCFS_brd
    
    #Plot levels
    lvs = [nat_brd(), 10, 100, 1000]
    cmap = mpl.cm.plasma
    
    fig1, ax1 = plt.subplots()
    contot = ax1.contourf(n,T,Colltot_brd, levels = lvs, cmap = mpl.cm.plasma,\
                          norm = mpl.colors.BoundaryNorm(lvs, cmap.N, extend = "both"),\
                          extend = "both")
    
    cbartot = fig1.colorbar(contot)
    cbartot.ax.set_ylabel("Linewidth broadening [Hz]")
    ax1.set_title("Sum collision broadening")
    ax1.set_xlabel("$n$ [$\mathrm{1/cm^3}$]")
    ax1.set_ylabel("$T$ [K]")
    
    ax1.ticklabel_format(style = "sci",axis = "both", scilimits = (0,0))
    plt.xscale("log")
    plt.grid(linestyle = "--")
    
    fig2, ax2 = plt.subplots()
    conCQ = ax2.contourf(n,T,CQ_brd, levels = lvs, cmap = mpl.cm.plasma,\
                          norm = mpl.colors.BoundaryNorm(lvs, cmap.N, extend = "both"),\
                          extend = "both")
    
    cbarCQ = fig2.colorbar(conCQ)
    cbarCQ.ax.set_ylabel("Linewidth broadening [Hz]")
    ax2.set_title("Collision quenching broadening")
    ax2.set_xlabel("$n$ [$\mathrm{1/cm^3}$]")
    ax2.set_ylabel("$T$ [K]")
    
    ax2.ticklabel_format(style = "sci",axis = "both", scilimits = (0,0))
    plt.xscale("log")
    plt.grid(linestyle = "--")
    
    fig3, ax3 = plt.subplots()
    conCCFS = ax3.contourf(X,Y,CCFS_brd, levels = lvs, cmap = mpl.cm.plasma,\
                          norm = mpl.colors.BoundaryNorm(lvs, cmap.N, extend = "both"),\
                          extend = "both")
    
    cbarCCFS = fig3.colorbar(conCCFS)
    cbarCCFS.ax.set_ylabel("Linewidth broadening [Hz]")
    ax3.set_title("Cold collision frequency shift broadening")
    ax3.set_xlabel("$n$ [$\mathrm{1/cm^3}$]")
    ax3.set_ylabel("$T$ [K]")
    
    ax3.ticklabel_format(style = "sci",axis = "both", scilimits = (0,0))
    plt.xscale("log")
    plt.grid(linestyle = "--")
    
    plt.show()

def Zeeman_limits():
    # Test Zeeman broadening
    
    twoD = False
    
    if twoD:
        samples = 2000
        Bw = np.linspace(0, 30, samples+1)
        T = np.geomspace(0, 1.05E-3, samples+1)
        # Calculate 2D space of Collision quenching
        Zm_brd = np.zeros((len(B),len(Ts)))
        for i in range(len(B)):
            for j in range(len(Ts)):
                Zm_brd[j][i] = Zeeman_brd(Bw = B[i], T = Ts[j])
        
        
        #Plot levels
        lvs = [nat_brd(), 10, 100, 1000]
        cmap = mpl.cm.plasma
        
        fig1, ax1 = plt.subplots()
        contot = ax1.contourf(B,Ts,Zm_brd, levels = lvs, cmap = mpl.cm.plasma, \
                            norm = mpl.colors.BoundaryNorm(lvs, cmap.N, extend = "both"),\
                            extend = "both")
        
        cbar = fig1.colorbar(contot)
        cbar.ax.set_ylabel("Linewidth broadening [Hz]")
        ax1.set_title("Zeeman broadening")
        ax1.set_xlabel("$B_w$ [G]")
        ax1.set_ylabel("$T$ [K]")
        ax1.ticklabel_format(style = "sci",axis = "both", scilimits = (0,0))
    
    else:
        samples = 10000
        T = np.linspace(9E-6, 1.05E-3, samples+1)
        Bw = T*BT
        Zm_brd12 = np.zeros(len(T))
        Zm_brd4 = np.zeros(len(T))
        for i in range(len(T)):
            Zm_brd12[i] = Zeeman_brd(Bw = Bw[i], poles = 12, w0 = 1000E-6, worst_case = False)
            Zm_brd4[i] = Zeeman_brd(Bw = Bw[i], poles = 4, w0 = 1000E-6, worst_case = False)
        
        
        fig1, ax1 = plt.subplots()
        ax1.plot(T,Zm_brd12, label = "12 pole")
        ax1.plot(T,Zm_brd4, label = "4 pole")
        
        ax1.set_xlabel("$T$ [K]")
        ax1.set_ylabel("$\Delta\omega_{1S-2S}$ contribution [Hz]")
        secax = ax1.secondary_xaxis("top", functions = (lambda x: x*BT*5, lambda x: x/BT/5))
        secax.set_xlabel("Wall field $B_w$ [G]")
        
        #ax1.ticklabel_format(style = "sci", scilimits = (-2,5))
        #ax.set_yticks([nat_brd(), 10, 25, 50, 75, 100, 125, 150, 175])
        #secax.set_xticks([25,15,10,5,3,2])
        plt.legend()
        plt.grid(which = "both",linestyle = "--")
        plt.xscale("log")
    
    plt.grid(linestyle = "--")
    
    plt.show()

def TOF_limits():
    # Test tof broadening
    w0 = np.linspace(200E-6, 2E-3, samples+1)
    T = np.geomspace(1E-6, 1.05E-3, samples+1)
    
    # Calculate 2D space of tof
    tof_brd = np.multiply.outer(log(2)/pi*(2*kBm*T/m_H)**(1/2),1/w0)
    
    #Plot levels
    lvs = [nat_brd(), 10, 100, 1000]
    cmap = mpl.cm.plasma
    
    fig1, ax1 = plt.subplots()
    contot = ax1.contourf(w0*1E6,T,tof_brd, levels = lvs, cmap = mpl.cm.plasma, \
                          norm = mpl.colors.BoundaryNorm(lvs, cmap.N, extend = "both"),\
                          extend = "both")
    
    cbar = fig1.colorbar(contot)
    cbar.ax.set_ylabel("Linewidth broadening [Hz]")
    ax1.set_title("TOF broadening")
    ax1.set_xlabel("$w_0$ [$\mathrm{\mu m}$]")
    ax1.set_ylabel("$T$ [K]")
    
    ax1.ticklabel_format(style = "sci",axis = "both", scilimits = (0,4))
    plt.yscale("log")
    plt.grid(linestyle = "--")
    
    plt.show()

def find_limits():
    
    intensity = 0
    stark = 0
    dopp = 0
    coll = 0
    zeeman = 1
    tof = 0
    
    if intensity:
        I_limits()
    
    if stark:
        E_limits()
    
    if dopp:
        Dopp_limits()
        
    if coll:
        Coll_limits()
    
    if zeeman:
        Zeeman_limits()
    
    if tof:
        TOF_limits()
    
    