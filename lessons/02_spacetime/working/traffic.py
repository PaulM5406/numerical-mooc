import numpy
import sympy
from matplotlib import pyplot
from sympy.utilities.lambdify import lambdify

# Set the font family and size to use for Matplotlib figures.
pyplot.rcParams['font.family'] = 'serif'
pyplot.rcParams['font.size'] = 16

sympy.init_printing()

# Set parameters.
# Vmax = 80
Vmax = 136
nx = 51  # number of spatial grid points
L = 11 # length of the domain
rho_max = 250
dt = 0.001 # time-step size in hours
dx = L / (nx - 1)  # spatial grid size
tmax = 6/60. # In hours
nt = int(tmax/dt)  # number of time steps to compute

#sigma = 0.1  # CFL limit
#dt = sigma * dx**2 / nu  # time-step size

rho = sympy.symbols('rho')
F = Vmax*rho*(1 - rho/rho_max)
dF_drho = F.diff(rho)
dF_drho = lambdify((rho), dF_drho)

def V(rho):
    return Vmax*(1 - rho/rho_max)

# Set initial conditions.
t = [0]
x = numpy.linspace(0, L, nx)
# rho0 = numpy.ones(nx)*10
rho0 = numpy.ones(nx)*20
rho0[10:20] = 50

# Integrate the equation in time.
rho = numpy.zeros((nx, nt+1))
rho[:, 0] = rho0
for n in range(nt):
    t.append(dt*(n+1))
    # Update all interior points.
    rho[1:, n+1] = rho[1:, n] - \
               dF_drho(rho[1:, n])*dt/dx*(rho[1:, n] - rho[:-1, n])
    # Update boundary points.
    # rho[0, n+1] = 10
    rho[0, n+1] = 20
       
Vmin_t0 = min(V(rho[:, 0]))
print(f"Vmin_t0 = {Vmin_t0/3.6:.2f}m/s")
Vmean_t180 = numpy.mean(V(rho[:, t.index(3/60)]))
print(f"Vmean_t180 = {Vmean_t180/3.6:.2f}m/s")
Vmin_t180 = min(V(rho[:, t.index(3/60)]))
print(f"Vmin_t180 = {Vmin_t180/3.6:.2f}m/s")
Vmin_t360 = min(V(rho[:, t.index(6/60)]))
print(f"Vmin_t360 = {Vmin_t360/3.6:.2f}m/s")

# Plot the numerical solution along with the analytical solution.
pyplot.figure(figsize=(6.0, 4.0))
pyplot.xlabel('x')
pyplot.ylabel('rho')
pyplot.grid()
pyplot.plot(x, rho[:, 0], label="Initial state",
            color='C0', linestyle='-', linewidth=2)
ntp = t.index(6/60)
pyplot.plot(x, rho[:, ntp], label=f"nt={ntp}",
            color='C1', linestyle='-', linewidth=2)
pyplot.legend()
pyplot.xlim(0.0, L)
#pyplot.ylim(0.0, 10.0);
pyplot.show()
pyplot.clf()
