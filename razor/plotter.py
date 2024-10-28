import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit
import os
from tkinter import filedialog
# Used for the ellipses for now
from matplotlib.patches import Ellipse

def choose_files():
    """
    Opens a file dialog to select .dat files.
    Returns a list of file paths.
    """
    filex = filedialog.askopenfilename(initialdir=".", filetypes=[('DAT Files', '*.dat'), ('Text Files', '*.txt')])
    filey = filedialog.askopenfilename(initialdir=".", filetypes=[('DAT Files', '*.dat'), ('Text Files', '*.txt')])
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

    # Sampling
    samples = 1000
    polydepth = 8
    offsetx1, offsetx2, offsety1, offsety2 = 0, 0, -2, 0
    xx_samp = np.linspace(xx.min(), xx.max(), samples)
    xy_samp = np.linspace(xy.min(), xy.max(), samples)
    X_samp, Y_samp = np.meshgrid(xx_samp, xy_samp)


    # Fit Gaussian to the x data
    popt_x, _ = curve_fit(gaussian, xx, dy_dx_x, p0=[1, np.mean(xx), np.std(xx)])
    mean_x, stddev_x = popt_x[1], popt_x[2]

    # Fit Gaussian to the y data
    popt_y, _ = curve_fit(gaussian, xy, dy_dx_y, p0=[1, np.mean(xy), np.std(xy)])
    mean_y, stddev_y = popt_y[1], popt_y[2]

    # Deduct the gaussian fit from the data
    dy_dx_x_poly = dy_dx_x - gaussian(xx, *popt_x)
    dy_dx_y_poly = dy_dx_y - gaussian(xy, *popt_y)
    
    # Divide the data into two subsets around mean for polynomial fits
    xx1 = xx[xx < mean_x-abs(2*stddev_x)]
    xx2 = xx[xx > mean_x+abs(2*stddev_x)]
    xy1 = xy[xy < mean_y-abs(2*stddev_y)]
    xy2 = xy[xy > mean_y+abs(2*stddev_y)]
    dy_dx_x_poly1 = dy_dx_x_poly[xx < mean_x-abs(2*stddev_x)]
    dy_dx_x_poly2 = dy_dx_x_poly[xx > mean_x+abs(2*stddev_x)]
    dy_dx_y_poly1 = dy_dx_y_poly[xy < mean_y-abs(2*stddev_y)]
    dy_dx_y_poly2 = dy_dx_y_poly[xy > mean_y+abs(2*stddev_y)]

    # Fit polynomials to the x data
    if len(xx1) > 1:
        polyfit_x1 = np.polyfit(xx1, dy_dx_x_poly1, polydepth+offsetx1)
    if len(xx2) > 1:
        polyfit_x2 = np.polyfit(xx2, dy_dx_x_poly2, polydepth+offsetx2)
    # Fit polynomials to the y data
    if len(xy1) > 1:
        polyfit_y1 = np.polyfit(xy1, dy_dx_y_poly1, polydepth+offsety1)
    if len(xy2) > 1:
        polyfit_y2 = np.polyfit(xy2, dy_dx_y_poly2, polydepth+offsety2)

    polyfit_x = np.zeros(samples)
    polyfit_y = np.zeros(samples)
    # Combine the two polynomial fits in their respective regions
    if len(xx1) > 1:
        polyfit_x[xx_samp < mean_x-abs(2*stddev_x)] = np.polyval(polyfit_x1, xx_samp[xx_samp < mean_x-abs(2*stddev_x)])
    else:
        polyfit_x[xx_samp < mean_x-abs(2*stddev_x)] = 0
    if len(xx2) > 1:
        polyfit_x[xx_samp > mean_x+abs(2*stddev_x)] = np.polyval(polyfit_x2, xx_samp[xx_samp > mean_x+abs(2*stddev_x)])
    else:
        polyfit_x[xx_samp > mean_x+abs(2*stddev_x)] = 0
    if len(xy1) > 1:
        polyfit_y[xy_samp < mean_y-abs(2*stddev_y)] = np.polyval(polyfit_y1, xy_samp[xy_samp < mean_y-abs(2*stddev_y)])
    else:
        polyfit_y[xy_samp < mean_y-abs(2*stddev_y)] = 0
    if len(xy2) > 1:
        polyfit_y[xy_samp > mean_y+abs(2*stddev_y)] = np.polyval(polyfit_y2, xy_samp[xy_samp > mean_y+abs(2*stddev_y)])
    else:
        polyfit_y[xy_samp > mean_y+abs(2*stddev_y)] = 0

    polygaussian_x = gaussian(xx_samp, *popt_x) + polyfit_x
    polygaussian_y = gaussian(xy_samp, *popt_y) + polyfit_y


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

    ax2.scatter(xx, dy_dx_x, s= 5, label='Inverse of Derivative', color='purple')
    ax2.plot(xx_samp, gaussian(xx_samp, *popt_x), label='Gaussian Fit', color='g', linestyle='--')
    ax2.plot(xx_samp, polygaussian_x, label='Polynomial Fit', color='b', linestyle='--')
    ax2.set_xlabel('X-axis')
    ax2.set_ylabel('dY/dX')
    ax2.annotate(f'Waist X = {abs(2*stddev_x):.2f}\nCenter X = {mean_x:.2f}', xy=(0.05, 0.78), xycoords='axes fraction', fontsize=12, verticalalignment='bottom', horizontalalignment='left')
    ax2.axvline(x=mean_x, color='r', linestyle='--', label=f'')
    ax2.axvline(x=mean_x+abs(stddev_x*2), color='r', linestyle='--', label=f'Waist')


    ax3.plot(xy, yy, label='Data Y')
    ax3.set_ylabel('P [mW]')
    ax3.set_xlabel('Y [mum]')
    ax3.set_title(title+" Y-axis")

    ax4.scatter(xy, dy_dx_y, s= 5, label='Inverse of Derivative', color='purple')
    ax4.plot(xy_samp, gaussian(xy_samp, *popt_y), label='Gaussian Fit', color='g', linestyle='--')
    ax4.plot(xy_samp, polygaussian_y, label='Polynomial Fit', color='b', linestyle='--')
    ax4.set_xlabel('Y-axis')
    ax4.set_ylabel('dY/dX')
    ax4.annotate(f'Waist Y = {abs(2*stddev_y):.2f}\nCenter Y = {mean_y:.2f}', xy=(0.05, 0.78), xycoords='axes fraction', fontsize=12, verticalalignment='bottom', horizontalalignment='left')
    ax4.axvline(x=mean_y, color='r', linestyle='--', label=f'')
    ax4.axvline(x=mean_y+abs(stddev_y*2), color='r', linestyle='--', label=f'Waist')

    # Create 2D Gaussian projection from poly fits
    Zpoly = np.outer(polygaussian_x, polygaussian_y).transpose()

    # Create 2D Gaussian projection from gaussian fits
    Zfit = gaussian(X_samp, *popt_x) * gaussian(Y_samp, *popt_y)

    normal = Normalize(vmin=Zfit.min(), vmax = Zfit.max(), clip = False)

    ax5.imshow(Zpoly, extent=(xx.min(), xx.max(), xy.min(), xy.max()), norm = normal, origin='lower', cmap='viridis', aspect='equal')
    ax5.set_title('2D Gaussian Projection from poly fit')
    ax5.set_xlabel('X-axis')
    ax5.set_ylabel('Y-axis')
    # Add dotted lines to show the waist location
    ax5.axvline(x=mean_x, color='r', linestyle='--', label=f'')
    ax5.axhline(y=mean_y, color='r', linestyle='--', label=f'')
    # Draw the waist ellipse
    ellipsepoly = Ellipse((mean_x, mean_y), width=abs(2*stddev_x), height=abs(2*stddev_y), angle=0, fill=False, color='r', linestyle='--', label='Waist')
    ax5.add_patch(ellipsepoly)
    #ax5.add_patch(plt.Circle((mean_x, mean_y), abs(2*stddev_x), fill=False, color='r', linestyle='--', label='Waist'))

    ax6.imshow(Zfit, extent=(xx.min(), xx.max(), xy.min(), xy.max()), norm = normal, origin='lower', cmap='viridis', aspect='equal')
    ax6.set_title('2D Gaussian Projection from fit')
    ax6.set_xlabel('X-axis')
    ax6.set_ylabel('Y-axis')
    ax6.axvline(x=mean_x, color='r', linestyle='--', label=f'')
    ax6.axhline(y=mean_y, color='r', linestyle='--', label=f'')
    # Draw the waist ellipse
    ellipsegaussian = Ellipse((mean_x, mean_y), width=abs(2*stddev_x), height=abs(2*stddev_y), angle=0, fill=False, color='r', linestyle='--', label='Waist')
    ax6.add_patch(ellipsegaussian)
    #ax6.add_patch(plt.Circle((mean_x, mean_y), abs(2*stddev_x), fill=False, color='r', linestyle='--', label='Waist'))

    plt.tight_layout()
    plt.show()

def main():
    # Directory containing .dat files
    data_dir = '.'
    
    files = choose_files()
    data_x = read_dat_file(files["x"])
    data_y = read_dat_file(files["y"])
    print(f"Data read successfully, files:\n{files['x']}\n{files['y']}")
    plot_data(data_x, data_y, title = "Beam profile")


if __name__ == "__main__":
    main()