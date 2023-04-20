using Pkg
Pkg.activate(".")
include("../src/AutoChemo.jl")

using JLD2

function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))

    return X, Y
end
L = 20
A = 1.0   # Amplitude
v = 2/L
k = 2 * pi* v    # Wave vector
w = 2 * pi / 0.2  # Angular frequency
phi0 = 0.0
t = 2

# Define the x and y coordinates for the grid
x = y = range(-10, stop=10, length=100)

# Create a meshgrid from the x and y coordinates
X, Y = meshgrid(x, y)

# Calculate the value of the wave at each point in the grid
Z = A * sin.(k * X .+ k*w*t .+ phi0)

# Plot the wave as a 3D surface
# surface(x, y, Z)

# Create a heatmap of the wave
heatmap(x, y, Z, aspect_ratio=:equal)