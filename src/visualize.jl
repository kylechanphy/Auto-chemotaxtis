
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

    plot(x[1:N], y[1:N], label="",
        aspect_ratio=1)
end


function set_equla_aspect!(fig)
    x12, y12, z12 = xlims(fig), ylims(fig), zlims(fig)
    d = maximum([diff([x12...]), diff([y12...]), diff([z12...])])[1] / 2
    xm, ym, zm = mean(x12), mean(y12), mean(z12)

    Plots.plot!(fig, xlims=(xm - d, xm + d), ylims=(ym - d, ym + d), zlims=(zm - d, zm + d))
end


function viz(pos::Vector{SVector{3, Float64}}, sysPara)
    dt = sysPara.dt
    x = [v[1] for v in pos]
    y = [v[2] for v in pos]
    z = [v[3] for v in pos]
    N = length(x)
    t = 0:dt:dt*(N-1)
  
    fig = plot(x,y,z, label="",
        xlabel="x", ylabel="y", zlabel="z",
        aspect_ratio=:equal,
        line_z=t,
        color=:viridis)
    # plot!(aspect_ratio=:equal)
    set_equla_aspect!(fig)
    # scatter!([x[1]], [y[1]], [z[1]])
end



# function plot_iso3d(xs, ys, zs; lw=3, label=false)
#     # condition data for nearly isometric 3D plot 
#     x12, y12, z12 = extrema(xs), extrema(ys), extrema(zs)
#     d = maximum([diff([x12...]), diff([y12...]), diff([z12...])])[1] / 2
#     xm, ym, zm = mean(x12), mean(y12), mean(z12)

#     # plot data
#     p = Plots.plot(; xlabel="x", ylabel="y", zlabel="z", aspect_ratio=:equal, grid=:true)
#     Plots.plot!(xlims=(xm - d, xm + d), ylims=(ym - d, ym + d), zlims=(zm - d, zm + d))
#     # Plots.plot!(; camera=(45,30))    #(azimuth,elevation) ???
#     Plots.plot!(xs, ys, zs, line_z = t, c=:viridis, label=label)
#     # Plots.plot!(xs, ys, zlims(p)[1] .+ 0 * zs, lw=1, label=false)
#     # Plots.plot!(xs, ylims(p)[2] .+ 0 * ys, zs, lw=1,  label=false)
#     # Plots.plot!(xlims(p)[1] .+ 0 * xs, ys, zs, lw=1, lc=:lightgray, label=false)
# end

# function viz(pos::Vector{SVector{3, Float64}}, sysPara)
#     GLMakie.activate!()
#     dt = sysPara.dt
#     x = [v[1] for v in pos]
#     y = [v[2] for v in pos]
#     z = [v[3] for v in pos]
#     N = length(x)
#     t = 0:dt:dt*(N-1)

#     aspect = (1, 1, 1)
#     perspectiveness = 0.5
#     # the figure
#     fig = Figure()
#     ax1 = Axis3(fig[1, 1]; aspect = (1,1,1))
#     l = lines!(ax1, x, y, z; color=t, linewidth=3)
#     Colorbar( fig[1,2], l; label="time")
#     scale!(ax1.scene, (1,1,1)...)

#     fig
# end

function vizXY(pos::Vector{SVector{3,Float64}}, sysPara)
    x = [v[1] for v in pos]
    y = [v[2] for v in pos]
    N = length(x)
    dt = sysPara.dt
    t = 0:dt:dt*(N-1)
    
    plot(x,y, label="",
        xlabel="x", ylabel="y", left_margin=2mm,
         aspect_ratio=:equal,
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
        aspect_ratio=:equal,
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
        aspect_ratio=:equal,
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
