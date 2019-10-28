import numpy
from matplotlib import pyplot

# Set the font family and size to use for Matplotlib figures.
pyplot.rcParams['font.family'] = 'serif'
pyplot.rcParams['font.size'] = 16

# Set parameters.
nx = 150  # number of spatial discrete points
L = 2.0  # length of the 1D domain
dx = L / (nx - 1)  # spatial grid size
nt = 2500  # number of time steps
dt = 0.0001  # time-step size
c = 1.0  # convection speed

# Define the grid point coordinates.
x = numpy.linspace(0.0, L, num=nx)

# Set initial conditions with 1.0 everywhere (for now).
u0 = numpy.ones(nx)
# Get a list of indices where 0.5 <= x <= 1.0.
#mask = numpy.where(numpy.logical_and(x >= 0.5, x <= 1.0))
# Set initial condition u = 2.0 where 0.5 <= x <= 1.0.
#u0[mask] = 2.0

u0[(x >= 0.5) & (x <= 1.0)] = 2.0
#print(u0)

# Compute the solution using Euler's method and array slicing.
u = u0.copy()
for n in range(1, nt):
  #u[1:] = u[1:] - dt / dx * c * (u[1:] - u[:-1])
    u[1:] = u[1:] - dt / dx * u[1:] * (u[1:] - u[:-1])

# Plot the solution after nt time steps
# along with the initial conditions.
pyplot.figure(figsize=(4.0, 4.0))
pyplot.xlabel('x')
pyplot.ylabel('u')
pyplot.grid()
pyplot.plot(x, u0, label='Initial',
            color='C0', linestyle='--', linewidth=2)
pyplot.plot(x, u, label='nt = {}'.format(nt),
            color='C1', linestyle='-', linewidth=2)
pyplot.legend()
pyplot.xlim(0.0, L)
pyplot.ylim(0.0, 2.5)
pyplot.show()
pyplot.clf()

