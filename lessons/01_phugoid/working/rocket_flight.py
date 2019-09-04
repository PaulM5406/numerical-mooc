# Modules
from math import pi
import numpy as np
import matplotlib.pyplot as plt

# Inputs
ms = 50
g = 9.81
rho = 1.091
r = 0.5
A = pi*r**2
ve = 325
Cd = 0.15
mp0 = 100

dt = 0.1
tmax = 60
v0 = 0
h0 = 0

# Functions
def mpp(t):
  if t < 5:
    return 20
  else:
    return 0
  
def dv(t, mp, v):
  return dt*(-g + (mpp(t)*ve - 0.5*rho*v*abs(v)*A*Cd)/(ms + mp))
  
# Main
t_list = [0]
v_list = [v0]
h_list = [h0]
mp_list = [mp0]
for t in np.arange(0, tmax, dt):
  t_list.append(t+dt)
  if t < 1e-7:
    mp = mp0
    v = v0
    h = h0
  h += v*dt
  h_list.append(h)
  v += dv(t, mp, v)
  v_list.append(v)
  mp -= mpp(t)*dt
  mp_list.append(mp)
  
plt.figure(1)
plt.title("mp")
plt.plot(t_list, mp_list)
plt.show()
plt.clf()
plt.figure(2)
plt.title("v")
plt.plot(t_list, v_list)
plt.show()
plt.clf()
plt.figure(3)
plt.title("h")
plt.plot(t_list, h_list)
plt.show()
plt.clf()

ind = t_list.index(3.2)
print("mp(3.2) =", mp_list[ind])
v_max = max(v_list)
print("v_max =", v_max)
ind = v_list.index(v_max)
print("h_v_max =", h_list[ind])
print("t_v_max =", t_list[ind])
h_max = max(h_list)
print("h_max =", h_max)
ind = h_list.index(h_max)
print("t_h_max =", t_list[ind])
ind = np.where(np.array(h_list) <= 0)[0][2]
print("t_impact =", t_list[ind])
print("v_impact =", v_list[ind])

  

