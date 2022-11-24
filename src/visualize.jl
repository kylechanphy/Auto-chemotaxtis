
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
        xlims=(0, (Nx - 1) * dx), ylims=(0, (Ny - 1) * dy), aspect_ratio=1)
    x = [v[1] for v in all_p]
    y = [v[2] for v in all_p]
    plot!(circle(part.R, part.pos[1], part.pos[2]), label="", c=:green, fill=1)
    plot!(hm, x, y, label="", c=:white)
    return hm
end