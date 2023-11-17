using SymPy
using Plots

@vars x z x0 z0 E_0 k a


sympy.init_printing(pretty_print=False, use_unicode=False, wrap_line=False, no_global=True)

# Clear screen
print("\033[2J")

# Tomlinson potential V(x,z)
Vpotencial = E_0/2 * (1-cos(2π*x/a)) + k/2 * (z-x)^2

# Derivative of the Tomlinson potential with respect to x
V_x = diff(Vpotencial,x).subs(x,x0).subs(z,z0)

# Derivative of the Tomlinson potential with respect to z
V_z = diff(Vpotencial,z).subs(x,x0).subs(z,z0)

# Second derivative of the Tomlinson potential with respect to x
V_xx = diff(Vpotencial,x,x).subs(x,x0).subs(z,z0)

# Second derivative of the Tomlinson potential with respect to z
V_zz = diff(Vpotencial,z,z).subs(x,x0).subs(z,z0)

# Second derivative of the Tomlinson potential with respect to x and z
V_xz = diff(Vpotencial,x,z).subs(x,x0).subs(z,z0)

# Third derivative of the Tomlinson potential with respect to x
V_xxx = diff(Vpotencial,x,x,x).subs(x,x0).subs(z,z0)

# Third derivative of the Tomlinson potential with respect to z
V_zzz = diff(Vpotencial,z,z,z).subs(x,x0).subs(z,z0)

# Taylor approximation of Tomlinson potential
Vp_taylor = Vpotencial.subs(x,x0).subs(z,z0) + V_x*(x-x0) + V_z*(z-z0) + 1/2 * (V_xx*(x-x0)^2 + 2*V_xz*(x-x0)*(z-z0) + V_zz*(z-z0)^2) + V_xxx*(x-x0)^3/6 + V_zzz*(z-z0)^3/6
fn = lambdify(Vp_taylor.subs(E_0,1).subs(k,1).subs(a,1).subs(x0,0.1).subs(z0,1).subs(z,0.1))
xs = range(-1, stop=1.5, length=256)
plot(xs, fn.(xs), ylims=(-5,6), lw=4)

# Taylor approx  with z = 0.1   and   x0 = 0.1  and   z0 = 1 
print("\n\n\n\n\n=====================================================\n")
Vp_taylor.subs(z,0).subs(x,0)

simplify(Vp_taylor.subs(E_0,1).subs(k,1).subs(a,1).subs(x0,0.1).subs(z0,1).subs(z,0.1))
# Con esta expresión, podemos convertirla a 3/2  factorizando x^2 y multiplicando por x^1/2