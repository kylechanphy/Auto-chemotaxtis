function sinPlaneWave2D(sysPara)
    @unpack ny, nx, dy, dx = sysPara
    L = dx * (nx - 1)
    A = 5.0   # Amplitude
    v = 2 / L
    k = 2 * pi * v    # Wave vector
    w = 2 * pi / 0.2  # Angular frequency
    phi0 = 0.0
    t = 0

    # Define the x and y coordinates for the grid
    x = (1:nx) * dx
    y = (1:ny) * dy

    # Create a meshgrid from the x and y coordinates
    X, Y = meshgrid(x, y)

    # Calculate the value of the wave at each point in the grid
    Z = A * sin.(k * X .+ k * w * t .+ phi0)

    return Z
end

function sinPlaneWave2D!(exField ,sysPara, t)
    @unpack ny, nx, dy, dx = sysPara
    L = dx * (nx - 1)
    A = 5.0   # Amplitude
    v = 2 / L
    k = 2 * pi * v    # Wave vector
    w = 2 * pi / L # Angular frequency
    phi0 = 0.0
    # Define the x and y coordinates for the grid
    x = (1:nx) * dx
    y = (1:ny) * dy

    # Create a meshgrid from the x and y coordinates
    X, Y = meshgrid(x, y)

    # Calculate the value of the wave at each point in the grid
    exField[:,:] = A * sin.(k * X .+ k * w * t .+ phi0)
end
