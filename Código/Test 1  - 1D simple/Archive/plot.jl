using Plots
#plotly()
gr()

# Define the parameters
E_0 = 1
k = 1
a = 1

# Define the potential function
function V(x, z)
    return E_0/2 * (1-cos(2π*x/a)) + k/2 * (z-x)^2
end

# Create a plot of V as a function of x
x = range(-3*a, 3*a, length=100)
V_x = V.(x, 0.5a) # Evaluate V for a fixed value of z
plot(x, V_x, xlabel="x", ylabel="V(x)", label="", legend=:bottomright,layout=(1,1), size=(700, 700), lw=3, tickfontsize=20, guidefontsize=20, legendfontsize=20)

#=
# Create a plot of V as a function of z
z = range(-3*a, 3*a, length=100)
V_z = V.(0.5a, z) # Evaluate V for a fixed value of x
plot(z, V_z, xlabel="z", ylabel="V(z)", label="", legend=:bottomright,layout=(1,1), size=(700, 700), lw=3,tickfontsize=20, guidefontsize=20, legendfontsize=20)
=#

#=
#### 3D surface of V as a function of x and z
V_xz = [V(xi, zi) for xi in x, zi in z]
surface(x, z, V_xz, xlabel="x", ylabel="z", zlabel="V(x,z)", label="", legend=:bottomright,layout=(1,1), size=(700, 700))
=#


# Animate the plot from z = 0 to z = 2a
anim = @animate for zi in range(0, 2a, length=100)
    V_x = V.(x, zi)
    plot!(x, V_x, xlabel="x", ylabel="V(x)", label="", legend=:bottomright,layout=(1,1), size=(700, 700), lw=6, tickfontsize=20, guidefontsize=20, legendfontsize=20)

end
gif(anim, "V_x.gif", fps=15)


