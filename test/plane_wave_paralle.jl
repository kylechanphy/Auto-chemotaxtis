using Pkg
Pkg.activate(".")
include("../src/AutoChemo.jl")

using JLD2
using Base.Threads


# Set the number of threads to the maximum available on your system
Threads.nthreads()



# Define the grid for the 2D sinusoidal wave
x = range(-50, 50, length=1000)
y = range(-50, 50, length=1000)

# Preallocate memory for the 2D sinusoidal wave
wave = Array{Float64}(undef, length(y), length(x))

# Function for the 2D sinusoidal wave
function sinusoidal_wave(amplitude, wavelength, kx, ky, phase,  x, y)
    return amplitude * sin(kx * x + ky * y + phase)
end

function getField!(wave, x, y)
    # Define the parameters for the 2D sinusoidal wave
    amplitude = 1.0
    wavelength = 10.0
    kx = 2 * pi / wavelength
    ky = 2 * pi / wavelength
    ky = 0
    phase = 0.0
    # Use multithreading to create the 2D sinusoidal wave
    for j in 1:length(y)
        for i in 1:length(x)
            wave[j, i] = amplitude * sin(kx * x[i] + ky * y[j] + phase)
    end
end
    return wave
end


getField!(wave, x, y)
# Plot the 2D sinusoidal wave
heatmap(x, y, wave, color=:viridis, aspect_ratio =1, 
        xlabel="X", ylabel="Y", title="2D Sinusoidal Plane Wave")
