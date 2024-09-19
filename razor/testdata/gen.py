import numpy as np
import matplotlib.pyplot as plt

# Parameters for the Gaussian
mean = 0
stddev = 1
num_points = 300

second_peak = True

# Generate Gaussian data
x = np.linspace(-8, 8, num_points)
y = np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

# Generate secondary smaller gaussian with offset
y2 = np.exp(-((x - mean - 7) ** 2) / (2 * (stddev/2) ** 2))*0.5
if second_peak:
    print("Second peak added")
    y += y2

# Add small random noise
noise = np.random.normal(0, 0.05, num_points)
y_noisy = y + noise

# Plot the data
plt.plot(x, y_noisy, label='Noisy Gaussian Data')
plt.plot(x, y, label='Original Gaussian', linestyle='--')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.title('Gaussian Data with Noise')
plt.show()

#Calculate cumulative integral of the data
cumulative = np.cumsum(y_noisy)
plt.plot(x, cumulative, label='Cumulative Integral')

# Write cumulative integral data to a file
data = np.column_stack(((-1)*x-1.25, cumulative))
np.savetxt('cumulative_integral.txt', data)
