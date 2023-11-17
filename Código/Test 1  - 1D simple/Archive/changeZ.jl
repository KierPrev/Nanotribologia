using Plots
gr()

# Define the parameters
E_0 = 1
k = 1
a = 1

# Define the potential function
function V(x, z)
    return E_0/2 * (1-cos(2π*x/a)) + k/2 * (z-x)^2
end

# Define a function that creates a plot of V(x,z)
function plot_V_z(z)
    x = range(-3*a, 3*a, length=100)
    V_x = V.(x, z)
    p = plot(x, V_x, xlabel="x", ylabel="V(x)", label="", legend=:bottomright, 
             layout=(1,1), size=(700, 700), lw=3, tickfontsize=20, 
             guidefontsize=20, legendfontsize=20, ylims=(0,6))
    return p
end

#=
# Create an animation of V(x,z) for z varying from 0 to 3
anim = @animate for z in range(0, 3, length=100)
    plot_V_z(z)
end

gif(anim, "V_x.gif", fps=24)
=#

# Plot V(x,z), with z=0
plot_V_z(0)

V0(x) = E_0/2 * (1-cos(2π*x/a))
plot!(V0,lw=3)
