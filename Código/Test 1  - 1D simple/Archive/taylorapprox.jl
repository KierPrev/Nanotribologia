using Plots

# Define constants
E_0 = 1.0
a = 1.0
k = 1.0
x0 = 0.0
z0 = 0.0

# Define the function for V(x,z)
function V(x, z)
    # Compute the terms of the Taylor series explicitly
    term1 = E_0/2 * (1 - cos(2π*x0/a))
    term2 = π*E_0/a * sin(2π*x0/a) * (x - x0)
    term3 = k/2 * (z0 - x0)^2
    term4 = (2π/a)^2 * E_0/2 * cos(2π*x0/a) * (x - x0)^2
    term5 = -k * (x - x0) * (z - z0)
    term6 = k/2 * (z - z0)^2

    # Sum the terms to approximate V(x,z)
    V_approx = term1 + term2 + term3 + term4 + term5 + term6
    return V_approx
end

# Create a range of x values to plot over
x_range = -2a:0.01:2a

# Create a function that generates the plot for a given z value
function plot_V(z)
    # Plot the function for x ranging from -2a to 2a
    plot(x_range, x -> V(x, z), xlabel="x", ylabel="V(x,z)", label="V(x,z)",ylims=(0,600))
end

# Create an animation of the function for z ranging from 0 to 3
anim = @animate for z in 0:1:30
    plot_V(z)
end

# Save the animation as an mp4 file
mp4(anim, "V_animation.mp4")
