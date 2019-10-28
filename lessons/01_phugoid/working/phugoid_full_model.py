# Modules
import math
import numpy
from matplotlib import pyplot

# Inputs

# Set parameters.
g = 9.81  # gravitational acceleration (m.s^{-2})
vt = 4.9  # trim velocity (m.s)
CD = 1.0 / 5  # drag coefficient
CL = 1.0  # lift coefficient

# Set initial conditions.
x0 = 0.0  # horizontal position
y0 = 1.5  # vertical position (altitude)

# Functions
def rhs_phugoid(u, CL, CD, g, vt):
    """
    Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : list or numpy.ndarray
        Solution at the previous time step
        as a list or 1D array of four floats.
    CL : float
        Lift coefficient.
    CD : float
        Drag coefficient.
    g : float
        Gravitational acceleration.
    vt : float
        Trim velocity.
    
    Returns
    -------
    rhs : numpy.ndarray
        The right-hand side of the system
        as a 1D array of four floats.
    """
    v, theta, x, y = u
    rhs = numpy.array([-g * math.sin(theta) - CD / CL * g / vt**2 * v**2,
                       -g * math.cos(theta) / v + g / vt**2 * v,
                       v * math.cos(theta),
                       v * math.sin(theta)])
    return rhs
    
def euler_step(u, f, dt, *args):
    """
    Returns the solution at the next time step using Euler's method.
    
    Parameters
    ----------
    u : numpy.ndarray
        Solution at the previous time step
        as a 1D array of floats.
    f : function
        Function to compute the right-hand side of the system.
    dt : float
        Time-step size.
    args : tuple, optional
        Positional arguments to pass to the function f.
    
    Returns
    -------
    u_new : numpy.ndarray
        The solution at the next time step
        as a 1D array of floats.
    """
    u_new = u + dt * f(u, *args)
    return u_new

# Main
dt = 0.001  # time-step size

theta0_list = numpy.arange(-90, 95, 5)*math.pi/180  # trajectory angle
v0_list = [3]*theta0_list.size  # start at the trim velocity

# Create a list to store x_end at each input parameters
x_end_list = []
for v0, theta0 in zip(v0_list, theta0_list):
  # Create a list to store the solution at each time step.
  u = []
  # Set the initial conditions.
  u.append(numpy.array([v0, theta0, x0, y0]))

  # Time integration with Euler's method.
  n = 0
  yn = u[0][3]
  while yn > 0:
    u.append(euler_step(u[n], rhs_phugoid, dt, CL, CD, g, vt))
    yn = u[n + 1][3]
    n += 1
  x_end = u[-1][2]
  x_end_list.append(x_end)
  print("Longueur atteinte après lancé :", round(x_end, 2), "m", "pour y =", round(yn, 2), "m. Inputs: v0 =", round(v0, 2), "m/s, theta0 =", round(theta0*180/math.pi, 2))
   
# Get the glider's position over the time.
x = [un[2] for un in u]
y = [un[3] for un in u]

T = n*dt

# Plot the path of the glider.
pyplot.figure(figsize=(9.0, 4.0))
pyplot.title('Path of the glider (flight time = {})'.format(T))
pyplot.xlabel('x')
pyplot.ylabel('y')
pyplot.grid()
pyplot.plot(x, y, color='C0', linestyle='-', linewidth=2)
pyplot.show()
pyplot.clf()