from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl



Bw = 15 #G
w0 = 300E-6
BT = 15E3 # G/K
poles = 8
samples = 10001
#trap wall radius
rT = 8E-3

# Fine structure constant, -
alph = 0.0072973525628
# Planck constant, eV/Hz
h = 4.135667696E-15
# Bohr magneton, eV/G
muB = 5.7883818060E-9

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

def B(r):
    if poles == 12:
        return B0+abs(c*r**5)
    if poles == 8:
        return B0+abs(c*r**3)
    if poles == 4:
        return B0+abs(c*r)

r = np.linspace(-rT,rT, samples)
Br = B(r)

rs = r[(r>-rmax) & (r<rmax)]
rsy = np.zeros(len(rs))+B(rs[0])

rbeam = np.array([-w0/2,w0/2])
beambot = np.array([0,0])
beamtop = np.array([Br[0]*1.1,Br[0]*1.1])

#For filling
Bfill = np.zeros(len(rs))
k = 0
for i in range(len(Br)):
    if r[i] in rs:
        Bfill[k] = Br[i]
        k+=1

dB = B(rmax)-B0

brd = sft*dB

#Plotting start

fig, ax = plt.subplots()
ax.plot(r*1E3,Br, color = "red", label = "B(r)")
ax.plot(rs*1E3,rsy, color = "cyan", label = "H sample")

ax.set_title(f"Sample at $T$ = {T*1E3:.2} mK, $B_w$ = {Bw} G, $w_0$ = {w0*1E6} $\mathrm{{\mu m}}$\nusing {poles}-pole trap")
ax.set_xlabel("Radial location $r$ [mm]")
ax.set_ylabel("B [G]")

ax.fill_between(1E3*rs, Bfill, rsy,color = "cyan", alpha=0.2)
ax.fill_between(1E3*rbeam, beamtop, beambot,color = "violet", alpha=0.8)
plt.vlines([-rT*1E3,rT*1E3], 0, Br[0]*1.1, colors = "black")
plt.legend()
plt.show()

