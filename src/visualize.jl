
circle(R, x, y) = (θ = LinRange(0, 2π, 30);
(x .+ R .* cos.(θ), y .+ R .* sin.(θ)))

function viz(u, all_p, sysPara, part)
    Nx = sysPara.nx
    Ny = sysPara.ny
    dx = sysPara.dx
    dy = sysPara.dy

    hmx = (0:Nx-1) * dy
    hmy = (0:Ny-1) * dx
    hm = heatmap(hmx, hmy, transpose(u),
        xlims=(0, (Nx) * dx), ylims=(0, (Ny) * dy), aspect_ratio=1)
    x = [v[1] for v in all_p]
    y = [v[2] for v in all_p]
    plot!(circle(part.R, part.pos[1], part.pos[2]), label="", c=:green, fill=1)
    plot!(hm, x, y, label="", c=:white)
    return hm
end

function viz(u, all_p, sysPara)
    Nx = sysPara.nx
    Ny = sysPara.ny
    dx = sysPara.dx
    dy = sysPara.dy

    hmx = (0:Nx-1) * dy
    hmy = (0:Ny-1) * dx
    hm = heatmap(hmx, hmy, transpose(u),
        xlims=(0, (Nx) * dx), ylims=(0, (Ny) * dy), aspect_ratio=1)
    x = [v[1] for v in all_p]
    y = [v[2] for v in all_p]
    # plot!(circle(part.R, part.pos[1], part.pos[2]), label="", c=:green, fill=1)
    plot!(hm, x, y, label="", c=:white)
    return hm
end

function viz(pos::Vector; N=length(pos))
    x = [v[1] for v in pos]
    y = [v[2] for v in pos]

    plot(x[1:N],y[1:N], label="",
        aspect_ratio = 1)
end



function viz(pos::Vector{SVector{3, Float64}}, sysPara)
    dt = sysPara.dt
    x = [v[1] for v in pos]
    y = [v[2] for v in pos]
    z = [v[3] for v in pos]
    N = length(x)
    t = 0:dt:dt*(N-1)
  
    plot(x,y,z, label="",
        xlabel="x", ylabel="y", zlabel="z",
        aspect_ratio=1,
        line_z=t,
        color=:viridis)
    # scatter!([x[1]], [y[1]], [z[1]])
end

function vizXY(pos::Vector{SVector{3,Float64}}, sysPara)
    x = [v[1] for v in pos]
    y = [v[2] for v in pos]
    N = length(x)
    dt = sysPara.dt
    t = 0:dt:dt*(N-1)
    
    plot(x,y, label="",
        xlabel="x", ylabel="y", left_margin=2mm,
        aspect_ratio=1,
        line_z=t,
        color=:viridis)
end

function vizXZ(pos::Vector{SVector{3,Float64}}, sysPara)
    x = [v[1] for v in pos]
    z = [v[3] for v in pos]
    N = length(x)
    dt = sysPara.dt
    t = 0:dt:dt*(N-1)
    
    plot(x,z, label="",
        xlabel="x", ylabel="z", left_margin=2mm,
        aspect_ratio=1,
        line_z=t,
        color=:viridis)
end

function vizYZ(pos::Vector{SVector{3,Float64}}, sysPara)
    z = [v[3] for v in pos]
    y = [v[2] for v in pos]
    N = length(z)
    dt = sysPara.dt
    t = 0:dt:dt*(N-1)

    plot(y, z, label="",
        xlabel="y", ylabel="z", left_margin=2mm,
        aspect_ratio=1,
        line_z=t,
        color=:viridis)
end


function plot_color(data, xlab="", ylab="")
    # t = range(0.0,stop=2π,length=1000)
    # y = sin.(t)
    plot(data, xlabel=xlab, ylabel=ylab,
        linewidth=5,
        linez=range(0.0, stop=1.0, length=length(data)),
        c=:batlow,
        legend=false,
        colorbar=false,)
    # px = [v[1] for v in data]
    # py = [v[2] for v in data]
    scatter!(data,
        markersize=5,
        zcolor=range(0.0, stop=1.0, length=length(data)),
        c=:batlow)
end


function viz_color(pos::Vector ;marker=false, cbar=false)
    x = [v[1] for v in pos]
    y = [v[2] for v in pos]

    fig = plot(x, y, label="",
        aspect_ratio=1,
        linez=range(0.0, stop=1.0, length=length(y)),
        linewidth=5,
        c=:lightrainbow,
        legend=false,
        colorbar=cbar)

    if marker == true
        px = [v[1] for v in pos]
        py = [v[2] for v in pos]
        scatter!(fig, px, py,
            markersize=5,
            zcolor=range(0.0, stop=1.0, length=length(py)),
            c=:batlow)
    end
    return fig
end

function plot_Fc(logger, sysPara, part, sacle)
    fig = viz(logger.field, logger.pos, sysPara, part)
    x, y = logger.pos[end]
    u, v = logger.Fc[end]
    vec = 5 * [u, v] / norm([u, v])
    quiver!(fig, [x], [y], quiver=([vec[1]], [vec[2]]), color=:cyan)

    return fig
end
