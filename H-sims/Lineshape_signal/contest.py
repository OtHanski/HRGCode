from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def plotter():
    samples = 1000
    
    a = np.geomspace(1E0, 1E3, samples+1)
    b = np.linspace(0, 1.05, 106)
    dat =  np.multiply.outer(b,a)
    print(dat)
    
    X, Y = np.meshgrid(a, b)
    
    lvs = [1, 10, 100]
    
    fig1, ax1 = plt.subplots()
    contot = ax1.contour(X,Y,dat, levels = lvs)
    plt.grid(linestyle = "--")
    
    cbartot = fig1.colorbar(contot)
    plt.show()

plotter()