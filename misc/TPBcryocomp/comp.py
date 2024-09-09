import matplotlib.pyplot as plt

x = [0,1,2,3,4]
bluecold = [8.8e-6, 7.5e-6, 9.0e-6,2.0e-6, 4.8e-6]
UVcold = [1.85e-5, 9.4e-7, 1.8e-5, 5.1e-6, 3.5e-6]
bluewarm = [2.4e-4, 2.15e-4, 2.4e-4, 5.5e-5, 1.3e-4]
UVwarm = [2.1e-4, 4.5e-6, 2.0e-4, 1.2e-5, 9.8e-6]
UVwarm28 = [9.5e-4, 2.4e-5, 8.7e-4, 7.0e-5, 4.3e-5]

# Plot the data
plt.plot(x, bluecold, label='Blue Cold')
plt.plot(x, UVcold, label='UV Cold')
plt.plot(x, bluewarm, label='Blue Warm')
plt.plot(x, UVwarm, label='UV Warm')
plt.plot(x, UVwarm28, label='UV Warm 28')
# Add a legend
plt.legend()
plt.show()

# Plot data normalized to the first value
plt.plot(x, [y/bluecold[0] for y in bluecold], label='Blue Cold')
plt.plot(x, [y/UVcold[0] for y in UVcold], label='UV Cold')
plt.plot(x, [y/bluewarm[0] for y in bluewarm], label='Blue Warm')
plt.plot(x, [y/UVwarm[0] for y in UVwarm], label='UV Warm')
plt.plot(x, [y/UVwarm28[0] for y in UVwarm28], label='UV Warm 28')
# Add a legend
plt.legend()
plt.show()

