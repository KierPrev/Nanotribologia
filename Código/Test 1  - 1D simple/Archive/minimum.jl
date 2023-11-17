using Plots

# Define the parameters
E_0 = 1.0
k = 2.0
a = 1.0
x0 = 0.5
z0 = 0.5
z = 1.0

# Define the function V(x,z)
V(x, z) = E_0/2 * (1 - cos(2π*x/a)) + k/2 * (z-x)^2

# Define the function for the minimum of V(x)
xmin(z) = (π^2 * E_0 / a^2 * x0 - k*(z-z0)) / (π^2 * E_0 / a^2 + k*(z-z0))

# Plot the function V(x,z)
xrange = range(0, stop=a, length=100)
zrange = range(0, stop=2, length=100)
plot(xrange, zrange, V, st=:surface, xlabel="x", ylabel="z", zlabel="V(x,z)", title="Function V(x,z)")

# Plot the minimum of V(x) for z=1.0
x_min = xmin(1.0)
scatter!([x_min], [1.0], [V(x_min, 1.0)], label="minimum of V(x)", color=:red)
