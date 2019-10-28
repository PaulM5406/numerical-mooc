import sympy
import numpy
from sympy.utilities.lambdify import lambdify

x = sympy.symbols('x')
phi = sympy.cos(x)**2*sympy.sin(x)**3/4/x**5/sympy.exp(x)
phi_diff = phi.diff(x)

phi_lamb = lambdify(x, phi_diff)
print(phi_lamb(2.2))