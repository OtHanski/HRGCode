import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

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

def plot_data(data, title):
    """
    Plots the data and its derivative using matplotlib.
    Assumes data is a 2D numpy array with two columns.
    """
    x = data[:, 0]
    y = data[:, 1]
    dy_dx = np.gradient(y, x)

    # Fit Gaussian to the x data
    popt_x, _ = curve_fit(gaussian, x, y, p0=[1, np.mean(x), np.std(x)])
    # Fit Gaussian to the derivative
    popt_dy_dx, _ = curve_fit(gaussian, x, dy_dx, p0=[1, np.mean(x), np.std(x)])

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 12))

    ax1.plot(x, y, label='Data')
    ax1.plot(x, gaussian(x, *popt_x), label='Gaussian Fit', color='g', linestyle='--')
    ax1.set_ylabel('Y-axis')
    ax1.set_title(title)
    ax1.legend()

    ax2.plot(x, dy_dx, label='Derivative', color='r')
    ax2.plot(x, gaussian(x, *popt_dy_dx), label='Gaussian Fit', color='g', linestyle='--')
    ax2.set_xlabel('X-axis')
    ax2.set_ylabel('dY/dX')
    ax2.legend()

    # Create 2D Gaussian projection
    X, Y = np.meshgrid(x, x)
    Z = gaussian(X, *popt_x) * gaussian(Y, *popt_dy_dx)

    ax3.imshow(Z, extent=(x.min(), x.max(), x.min(), x.max()), origin='lower', cmap='viridis')
    ax3.set_title('2D Gaussian Projection')
    ax3.set_xlabel('X-axis')
    ax3.set_ylabel('Y-axis')

    plt.tight_layout()
    plt.show()

def main():
    # Directory containing .dat files
    data_dir = '.'
    
    # Iterate over all .dat files in the directory
    for file_name in os.listdir(data_dir):
        if file_name.endswith('.dat'):
            file_path = os.path.join(data_dir, file_name)
            data = read_dat_file(file_path)
            plot_data(data, title=file_name)

if __name__ == "__main__":
    main()