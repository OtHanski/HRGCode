from math import *
import numpy as np
import matplotlib.pyplot as plt

# Simulation done with equation:
# 4/pi*n_0*k*D_ab^2*P_las*w_div/V_0*f(β,ω)

##################
### Parameters ###
##################
testmode = 0
multiplots = 1
lifetimes = 0

# Adjustables
w0 = [300E-6,750E-6] # [m]
n = [1E12] # 1/cm^3
l_s = [200E-3,200E-3,150E-3] # [m]
r_s = [25E-3,13E-3,2E-3] # [m]
T = [1E-4] # [K]
P = 5E-3 # [W]

# Detector solid angle percentage:
Detang = 0.05


# Const
# Laser wavelength and *** constant
lam = 243E-9 # [m]
k = 2*pi/lam # [1/m]
# Physical constants
Dab = 4.63E-4 # [m^2/J]
MH = 1.6735575E-27 # [kg]
kB = 1.380649E-23 # [kgm^2/Ks^2]
wtrans = 2.466061413187035E15 # [Hz]
Etrans = 6.62607015E-34*wtrans # [J/Hz]*[Hz]

# Plotting settings
samples = 20000
detune = 0.25E5 # max detuning, [Hz]
simtime = 1000 # tracking sample lifetime, [s]
logsc = False


#################

# Function definitions
def vr0(T):
    #rms transverse velocity, T in K
    return sqrt(2*kB*T/MH)

def simple_FWHM(v,w):
    return 1/pi*np.log(2)*v/w

def DoppFree(param):

    n,ls,w,T,detune = param[0],param[1],param[2],param[3],param[4]
    # Calculate line value for set laser detuning omeg 
    v = vr0(T)
    return (n*1E6)*ls*Dab**2*P**2/(w*v)*exp(-w*abs(detune)/v)

def NaiveLineShape():
    # Assuming lineshape widening is naively widened as sum of dOmega_part
    # R = 4/pi*n*k*D**2*P**2*arctan(l/(k*w**2))/dOmega
    pass

def calc_lineshape(param):
    
    n,ls,w,T,detune_max,samples = param[0],param[1],param[2],param[3],param[4],param[5]
    
    # define range
    x = np.linspace(-detune_max,detune_max,samples)
    # Calculate line values
    y = np.zeros(samples)
    for i in range(samples):
        y[i] = DoppFree([n,ls,w,T,x[i]])
        
    return [x,y,simple_FWHM(vr0(T),w),param]

def single_plot_lineshape(data,drw = True):
    
    x,y,FWHM,param = data[0],data[1],data[2],data[3]
    lab = f"n = {param[0]:.2g} $\\frac{{1}}{{\\mathrm{{cm^3}}}}$, T = {param[3]*1000} mK,\n$l_s$ = {param[1]*100} cm, $w_0$ = {param[2]*1E6} μm"
    
    fig,ax = plt.subplots()
    ax.plot(x/1000,y, label= lab)
    
    # Set axes
    ax.set_xlabel("$2ω-ω_{1S-2S}$ [kHz]")
    ax.set_ylabel("Excitation rate [1/s]")
    secax = ax.secondary_yaxis("right", functions = (lambda x: (Etrans*10**9)*x, lambda x: x/(Etrans*10**9)))
    secax.set_ylabel("$L_α$ fluorescence power [nW]")
    
    # Set notes
    ax.annotate(f"FWHM = {FWHM:.4} Hz",xy = (0,0), xytext = (0.03,0.95),textcoords = "axes fraction")
    
    if drw:
        if logsc:
            plt.yscale("log")
        plt.grid(linestyle = "--")
        plt.legend()
        plt.show()
    else:
        print("Returning")
        return ax
    

def multiplot(datas, drw = True):
    
    fig,ax = plt.subplots()
    
    # Set axes
    ax.set_xlabel("$2ω-ω_{1S-2S}$ [kHz]")
    ax.set_ylabel("Excitation rate [1/s]")
    secax = ax.secondary_yaxis("right", functions = (lambda x: (Etrans*10**9)*x, lambda x: x/(Etrans*10**9)))
    secax.set_ylabel("$L_α$ fluorescence power [nW]")
    FWHMtext = ""
    
    for data in datas:
        x,y,FWHM,param = data[0],data[1],data[2],data[3]
        lab = f"n = {param[0]:.2g} $\\frac{{1}}{{\\mathrm{{cm^3}}}}$, T = {param[3]*1000} mK,\n$l_s$ = {param[1]*100} cm, $w_0$ = {param[2]*1E6} μm"
        ax.plot(x/1000,y, label= lab)
        FWHMtext += f"FWHM = {FWHM:.4} Hz\n"
    
    # Set notes
    yt = 1-0.055*len(datas)
    ax.annotate(FWHMtext,xy = (0,0), xytext = (0.03,yt),textcoords = "axes fraction")
    
    if drw:
        if logsc:
            plt.yscale("log")
            plt.grid(which = "both", linestyle = "--")
        else:
            plt.grid(linestyle = "--")
        plt.legend()
        plt.show()
    return 0

def lifetime_sim(datas,drw = True):
    # datas is of form [[x,y,simple_FWHM(vr0(T),w),param],...]
    
    fig,ax = plt.subplots()
    # Set axes
    ax.set_xlabel("$t$ [s]")
    ax.set_ylabel("sample density")
    FWHMtext = ""
    
    lives = []
    for data in datas:
        param = data[3]
        n,ls,w,T,detune_max,samples = param[0],param[1],param[2],param[3],param[4],param[5]
        trange = np.linspace(0,simtime,samples)
        N = np.zeros(len(trange))
        V = r_s[2]**2*pi*ls*1E6    # cm^3
        N[0] = n*V
        
        for t in range(1,len(trange)):
            N[t] = N[t-1] - DoppFree([N[t-1]/V,ls,w,T,0,samples])*(trange[t]-trange[t-1])
        n = N/V
        
        lab = f"n = {param[0]:.2g} $\\frac{{1}}{{\\mathrm{{cm^3}}}}$, T = {param[3]*1000} mK,\n$l_s$ = {param[1]*100} cm, $w_0$ = {param[2]*1E6} μm"
        ax.plot(trange,n, label= lab)
        print(N)

    if drw:
        if logsc:
            plt.yscale("log")
            plt.grid(which = "both", linestyle = "--")
        else:
            plt.grid(linestyle = "--")
        plt.legend()
        plt.show()
    
    
def main():
    print("Running main")
    # param: [n,ls,w,T,omeg,samples]
    #param1 = [n[2],l_s[2],w0[0],T[2],detune,samples]
    if multiplots:
        datas = []
        for i in range(len(w0)):
            for j in range(len(T)):
                datas.append(calc_lineshape([n[0],l_s[2],w0[i],T[j],detune,samples]))
        
        multiplot(datas)
    
    if lifetimes:
        datas = []
        for i in range(len(w0)):
            for j in range(len(T)):
                datas.append(calc_lineshape([n[2],l_s[2],w0[i],T[j],detune,samples]))
        
        lifetime_sim(datas)

def test():
    # param: [n,ls,w,T,omeg,samples]
    param1 = [1E12,150E-3,40E-6,1E-4,4E5,20000]
    param2 = [1E12,150E-3,50E-6,1E-4,4E5,20000]
    param3 = [1E12,150E-3,60E-6,1E-4,4E5,20000]
    param4 = [1E12,150E-3,70E-6,1E-4,4E5,20000]
    param5 = [1E12,150E-3,80E-6,1E-4,4E5,20000]
    datas = [calc_lineshape(param1),calc_lineshape(param2),calc_lineshape(param3),calc_lineshape(param4),calc_lineshape(param5)]
    multiplot(datas)
    

if __name__ =="__main__":
    if testmode:
        print("Running test")
        test()
    else:
        main()



















