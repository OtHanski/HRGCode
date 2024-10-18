T = 100E-6

k = 1.38E-23
m = 1.67E-27
g = 9.81

def column_height_calc(T, k, m, g):
    return (3/2*k*T)/(m*g)

print(f"{1000*column_height_calc(T, k, m, g):.3f} mm")