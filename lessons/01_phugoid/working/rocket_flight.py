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

t0 = 0
mp0 = 100
v0 = 0
h0 = 0

dt = 0.01
tmax = 60

# Functions
def mpp(t):
  if t < 5:
    return 20
  else:
    return 0
    
def f(u, dt, g, ve, rho, A, Cd, ms):
  t, mp, v = u[0]+dt, u[1], u[2]
  dh = v
  dmp = -mpp(t)
  dv = (-g + (mpp(t)*ve - 0.5*rho*v*abs(v)*A*Cd)/(ms + mp))
  return np.array([1, dmp, dv, dh])
  
def euler(dt, u, f, *args):
  u_new = u + dt*f(u, dt, *args)
  return u_new
  
def rk2(dt, u, f, *args):
  u_star = u + dt/2*f(u, dt/2, *args)
  u_new = u + dt*f(u_star, dt/2, *args)
  return u_new

# Initialisation
u0 = np.array([t0, mp0, v0, h0])
n = int(tmax/dt) + 1
u = np.zeros((n+1, 4))
u[0] = u0

# Main
for i in range(n):
  u[i+1] = euler(dt, u[i], f, g, ve, rho, A, Cd, ms)
  #u[i+1] = rk2(dt, u[i], f, g, ve, rho, A, Cd, ms)
  
plt.figure(1)
plt.title("mp")
plt.plot(u[:, 0], u[:, 1])
plt.show()
plt.clf()
plt.figure(2)
plt.title("v")
plt.plot(u[:, 0], u[:, 2])
plt.show()
plt.clf()
plt.figure(3)
plt.title("h")
plt.plot(u[:, 0], u[:, 3])
plt.show()
plt.clf()

ind = np.where(np.abs(u[:, 0] - 3.2) < 1e-5)[0][0]
print("mp(3.2) =", u[ind, 1])
v_max = max(u[:, 2])
print("v_max =", v_max)
ind = list(u[:, 2]).index(v_max)
print("h_v_max =", u[ind, 3])
print("t_v_max =", u[ind, 0])
h_max = max(u[:, 3])
print("h_max =", h_max)
ind = list(u[:, 3]).index(h_max)
print("t_h_max =", u[ind, 0])
ind = np.where(u[:, 3] <= 0)[0][2]
print("t_impact =", u[ind, 0])
print("v_impact =", u[ind, 2])
