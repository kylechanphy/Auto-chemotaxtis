using PyCall
# import GLMakie as gm

@pyimport matplotlib
matplotlib.use("TkAgg") # Required 
@pyimport matplotlib.pyplot as plt
@pyimport numpy as np


function streamline(field::Vector{Vector{SVector{2,Float64}}})
    nx = length(field)
    ny = length(field[1])

    x, y = np.arange(0, nx), np.arange(0, ny)

    U, V = zeros(nx, ny), zeros(nx, ny)
    for i in 1:nx
        U[i, :] = [v[1] for v in field[i]]
        V[i, :] = [v[2] for v in field[i]]
    end

    U, V = transpose(U), transpose(V)
    speed = sqrt.(U .^ 2 + V .^ 2)
    fig, axs = plt.subplots(1, 1)
    axs.set_aspect("equal", "box")
    h = axs.streamplot(x, y, U, V, color=speed, cmap="seismic")
    cbar = fig.colorbar(h.lines, ax=axs)
    plt.show()
    return fig, axs
end


function streamline(field::Vector{Matrix{Float64}})
    dims, nx, ny = np.shape(field)
    x, y = np.arange(0, nx), np.arange(0, ny)
    U, V = field
    U, V = transpose(U), transpose(V)
    speed = sqrt.(U .^ 2 + V .^ 2)
    fig, axs = plt.subplots(1, 1)
    axs.set_aspect("equal", "box")
    h = axs.streamplot(x, y, U, V, color=speed, cmap="seismic")
    cbar = fig.colorbar(h.lines, ax=axs)

    plt.show(fig)
    # return fig, axs
end


function update(frame)
    field = randn(10, 10)
    ∇ = np.gradient(field)
    x, y = np.arange(0, nx), np.arange(0, ny)
    ax.streamplot(x, y, field[1], field[2])
end


function pymoive()
    fig, ax = plt.subplots()
    anim = FuncAnimation(fig, func=update, frames=1:20)
    plt.show()

end


