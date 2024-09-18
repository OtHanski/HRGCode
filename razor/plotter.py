import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit
import os
from tkinter import filedialog

def choose_files():
    """
    Opens a file dialog to select .dat files.
    Returns a list of file paths.
    """
    filex = filedialog.askopenfilename(initialdir=".", filetypes=[('DAT Files', '*.dat')])
    filey = filedialog.askopenfilename(initialdir=".", filetypes=[('DAT Files', '*.dat')])
    return {"x": filex, "y": filey}

def read_dat_file(file_path):
    """
    Reads a .dat file and returns the data as a numpy array.
    Assumes the .dat file has two columns of numerical data.
    """
    data = np.loadtxt(file_path)
    return data

def gaussian(x, amp, mean, stddev):
    """
    Gaussian function.
    """
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

def plot_data(datax, datay, title):
    """
    Plots the data and its derivative using matplotlib.
    Assumes data is a 2D numpy array with two columns.
    """
    xx = datax[:, 0]
    yx = datax[:, 1]
    xy = datay[:, 0]
    yy = datay[:, 1]
    dy_dx_x = (-1)*np.gradient(yx, xx)
    dy_dx_y = (-1)*np.gradient(yy, xy)

    # Fit Gaussian to the x data
    popt_x, _ = curve_fit(gaussian, xx, dy_dx_x, p0=[1, np.mean(xx), np.std(xx)])
    mean_x, stddev_x = popt_x[1], popt_x[2]

    # Fit Gaussian to the y data
    popt_y, _ = curve_fit(gaussian, xy, dy_dx_y, p0=[1, np.mean(xy), np.std(xy)])
    mean_y, stddev_y = popt_y[1], popt_y[2]

    # Calculate FWHM for Gaussian fits
    fwhm_x = 2 * np.sqrt(2 * np.log(2)) * popt_x[2]
    fwhm_y = 2 * np.sqrt(2 * np.log(2)) * popt_y[2]

    fig = plt.figure(figsize=(18, 6))
    # Create a grid of subplots
    gs = fig.add_gridspec(2, 8)
    # x-axis
    ax1 = fig.add_subplot(gs[0,0:2])
    ax3 = fig.add_subplot(gs[0,2:4])
    # y-axis
    ax2 = fig.add_subplot(gs[1,0:2])
    ax4 = fig.add_subplot(gs[1,2:4])
    # 2D projection
    ax5 = fig.add_subplot(gs[:, 4:6])
    ax6 = fig.add_subplot(gs[:, 6:8])
    #axes = [ax1, ax2, ax3, ax4, ax5]

    ax1.plot(xx, yx, label='Data X')
    ax1.set_ylabel('P [mW]')
    ax1.set_xlabel('X [mum]')
    ax1.set_title(title+" X-axis")

    ax2.plot(xx, dy_dx_x, label='Inverse of Derivative', color='r')
    ax2.plot(xx, gaussian(xx, *popt_x), label='Gaussian Fit', color='g', linestyle='--')
    ax2.set_xlabel('X-axis')
    ax2.set_ylabel('dY/dX')
    ax2.annotate(f'Waist X = {-2*stddev_x:.2f}', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=12, verticalalignment='bottom', horizontalalignment='left')
    ax2.axvline(x=mean_x, color='r', linestyle='--', label=f'')
    ax2.axvline(x=mean_x-(stddev_x*2), color='r', linestyle='--', label=f'Waist')


    ax3.plot(xy, yy, label='Data Y')
    ax3.set_ylabel('P [mW]')
    ax3.set_xlabel('Y [mum]')
    ax3.set_title(title+" Y-axis")

    ax4.plot(xy, dy_dx_y, label='Inverse of Derivative', color='r')
    ax4.plot(xy, gaussian(xy, *popt_y), label='Gaussian Fit', color='g', linestyle='--')
    ax4.set_xlabel('Y-axis')
    ax4.set_ylabel('dY/dX')
    ax4.annotate(f'Waist Y = {2*stddev_y:.2f}', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=12, verticalalignment='bottom', horizontalalignment='left')
    ax4.axvline(x=mean_y, color='r', linestyle='--', label=f'')
    ax4.axvline(x=mean_y+(stddev_y*2), color='r', linestyle='--', label=f'Waist')

    # Create 2D Gaussian projection from raw data
    X, Y = np.meshgrid(xx, xy)
    Zraw = gaussian(X, *popt_x) * gaussian(Y, *popt_y)
    # Create 2D Gaussian projection from fits
    samples = 1000
    xx_samp = np.linspace(xx.min(), xx.max(), samples)
    xy_samp = np.linspace(xy.min(), xy.max(), samples)
    X, Y = np.meshgrid(xx_samp, xy_samp)
    Zfit = gaussian(X, *popt_x) * gaussian(Y, *popt_y)

    normal = Normalize(vmin=Zraw.min(), vmax=Zraw.max(), clip = False)

    ax5.imshow(Zraw, extent=(xx.min(), xx.max(), xy.min(), xy.max()), norm = normal, origin='lower', cmap='viridis', aspect='equal')
    ax5.set_title('2D Gaussian Projection from raw data')
    ax5.set_xlabel('X-axis')
    ax5.set_ylabel('Y-axis')

    ax6.imshow(Zfit, extent=(xx.min(), xx.max(), xy.min(), xy.max()), norm = normal, origin='lower', cmap='viridis', aspect='equal')
    ax6.set_title('2D Gaussian Projection from fit')
    ax6.set_xlabel('X-axis')
    ax6.set_ylabel('Y-axis')

    plt.tight_layout()
    plt.show()

def main():
    # Directory containing .dat files
    data_dir = '.'
    
    files = choose_files()
    data_x = read_dat_file(files["x"])
    data_y = read_dat_file(files["y"])
    plot_data(data_x, data_y, title = "Beam profile")


if __name__ == "__main__":
    main()