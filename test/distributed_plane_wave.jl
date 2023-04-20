using Pkg
Pkg.activate(".")
include("../src/AutoChemo.jl")

using Distributed
addprocs(4) # add 4 worker processes

# Define the grid for the 2D sinusoidal wave
x = range(-50, 50, length=1000)
y = range(-50, 50, length=1000)

# Preallocate memory for the 2D sinusoidal wave
wave = Array{Float64}(undef, length(y), length(x))

# Function for the 2D sinusoidal wave
@everywhere function sinusoidal_wave(amplitude, wavelength, kx, ky, phase, x, y)
    return amplitude * sin.(kx .* x .+ ky .+ phase .+ y')
end

@everywhere function distributed_sinusoidal_wave(amplitude, wavelength, kx, ky, phase, x, y, start_idx, end_idx)
    wave = 0.0
    for i in start_idx:end_idx
        wave += sum(sinusoidal_wave(amplitude, wavelength, kx, ky, phase, x[i], y[i, :]))
    end
    return wave
end

function getField(wave, x, y)
    # Define the parameters for the 2D sinusoidal wave
    amplitude = 1.0
    wavelength = 10.0
    kx = 2 * pi / wavelength
    ky = 0
    phase = 0.0
    # Use distributed computing to create the 2D sinusoidal wave
    n = length(x)
    np = nworkers() # get number of worker processes
    chunksize = ceil(Int, n / np)
    futures = Vector{Future}()
    for i in 1:np
        start_idx = (i - 1) * chunksize + 1
        end_idx = min(i * chunksize, n)
        future_wave = @spawnat i distributed_sinusoidal_wave(amplitude, wavelength, kx, ky, phase, x, y, start_idx, end_idx)
        push!(futures, future_wave)
    end
    wave[:] = vcat(fetch.(futures)...)
    return wave
end

wave = similar(wave)
getField(wave, x, y)
