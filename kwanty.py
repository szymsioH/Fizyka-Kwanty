import sympy as sp

sp.init_printing()

x = sp.Symbol('x', real = True)
k = sp.Symbol('k', real = True, positive = True)
a = sp.Symbol('a', real = True)
n = sp.Symbol('n', integer = True, nonzero = True)
h = sp.Symbol('hbar', real = True, positive = True)
psi = sp.Function('\psi')
p = sp.Symbol('p')

eq = sp.Eq(psi(x).diff(x, x), -k**2 * psi(x))
sol = sp.dsolve(eq, psi(x), ics = {psi(0): 0})
sol = sol.subs(k, n * sp.pi / a)
sp.Symbol('C1')
sol = sol.subs(sp.Symbol('C1'), sp.sqrt(2 / a))

p=h * sp.I*-1
C = psi(x)* p * psi(x).diff(x)

C = C.subs(psi(x), sol.rhs)

ep = sp.integrate(C,(x,0,a)).simplify()

C = psi(x) * p**2 * psi(x).diff(x,x)

C = C.subs(psi(x), sol.rhs)

ep2 = sp.integrate(C,(x,0,a)).simplify()

print("Sprawdzenie czy zasada nieoznaczoności jest spełniona:")

ex = sp.integrate(sol.rhs*x*sol.rhs,(x,0,a)).simplify()
ex2 = sp.integrate(sol.rhs*x**2*sol.rhs,(x,0,a)).simplify()
a = (ex2-ex**2)*(ep2-ep**2)

(a>= h/2).simplify()