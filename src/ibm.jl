using Pkg
using DrWatson
@quickactivate
Pkg.instantiate()
Pkg.resolve()
using ActiveMatter

using Revise
using Parameters
using Plots
using BenchmarkTools
using StatsBase
using StaticArrays
using LinearAlgebra
using Interpolations
# using Parameters
# using BasicInterpolators:BicubicSplineInterpolator
include("Pyplot.jl")
# include("diffuse_test.jl")

@with_kw mutable struct SysPara
    dx::Float64 = 1
    dy::Float64 = 1
    dt::Float64 = 0.05
    nx::Int64 = 1000
    ny::Int64 = 1000
    Nstep::Int64 = 10

    npoly::Int64 = 64

    value = 0
end


@with_kw mutable struct Particle
    v0::Float16 = 10
    ω0::Float64 = 1
    Dr::Float64 = 0
    D::Float64 = 1
    R::Float64 = 10
    src::Float64 = 1.0
    α::Float64 = 1.0

    pos::SVector = SA[0.0, 0.0]
    vel::SVector = SA[0.0, 0.0]
    ϕ::Float64 = 0.0
end

# @with_kw mutable struct Particle
#     pos::SVector = SA[0., 0.]
#     ϕ::Float64 = 0.0
# end

function constSrc!(du, sysPara, part)
    @unpack nx, ny, dt = sysPara
    @unpack pos, R, src = part
    Threads.@threads for i in 1:nx
        for j in 1:ny
            if norm(pos - SA[i-1, j-1]) < R
                du[i, j] = src
            end
        end
    end
end

function constSrc2!(du, sysPara, part)
    @unpack nx, ny, dx, dy, dt = sysPara
    @unpack pos, R, src = part

    Lscale = SA[dx, dy]
    x, y = pos
    Threads.@threads for i in 1:nx
        for j in 1:ny
            if norm(pos - SA[i-1, j-1] .* Lscale) < R

                du[i, j] = src * ibm4c(abs(x - (i - 1))) * ibm4c(abs(y - (j - 1)))
            end
        end
    end
end

function constFlux2!(du, sysPara, part)
    @unpack nx, ny, dx, dy, dt = sysPara
    @unpack pos, R, src = part
    # R = R * 1.5
    Lscale = SA[dx, dy]
    x, y = pos 
    
    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)


    Threads.@threads for i in xlimlo:xlimup
        for j in ylimlo:ylimup
            # @show i,j
            du[i, j] = du[i, j] + dt * src * ibm4c(abs(x - (i - 1) * dx) / dx) * ibm4c(abs(y - (j - 1) * dy) / dy) / dx^2
        end
    end
end

function constFlux3!(du, sysPara, part)
    @unpack nx, ny, dx, dy, dt = sysPara
    @unpack pos, R, src, vel, ϕ= part
    # R = R * 1.5
    Lscale = SA[dx, dy]
    if vel == SA[0.,0.]
        vel = SA[cos(ϕ), sin(ϕ)]
    end
    x, y = pos .- R*0.8* vel ./ norm(vel)

    R = R * 1.5

    xlimlo = floor(Int, (x - 2R) / dx + 1)
    xlimup = ceil(Int, (x + 2R) / dx + 1)
    ylimlo = floor(Int, (y - 2R) / dy + 1)
    ylimup = ceil(Int, (y + 2R) / dy + 1)


    Threads.@threads for i in xlimlo:xlimup
        for j in ylimlo:ylimup
            # @show i,j
            du[i, j] = du[i, j] + dt * src * ibm4c(abs(x - (i - 1) * dx) / dx) * ibm4c(abs(y - (j - 1) * dy) / dy) / dx^2
        end
    end
end


function updataGrid!(u, du, sysPara, part)
    @unpack dx, dy, dt, nx, ny = sysPara
    @unpack pos, R, D = part
    Threads.@threads for i in 2:nx-1
        for j in 2:ny-1
            # if norm(pos - (SA[i-1, j-1])) < R
            #     du[i, j] = u[i, j]

            # else
            #  !isInside(R, pos, SA[i, j])
            du[i, j] = u[i, j] + dt * D * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / dx^2
                                           +
                                           (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / dy^2)
            # end
        end
    end
end

function updataGridAdvection!(u, du, flow, sysPara, part)
    @unpack dx, dy, dt, nx, ny = sysPara
    @unpack pos, R, D = part
    Threads.@threads for i in 2:nx-1
        for j in 2:ny-1

            du[i, j] = u[i, j] + dt * D * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / dx^2 - ((u[i+1, j] - u[i-1, j] / 2dx)) * flow[i][j][1]
                                           +
                                           (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / dy^2 - ((u[i, j+1] - u[i, j-1] / 2dy)) * flow[i][j][2])
        end
    end
end
function updataGrid4Order!(u, du, sysPara, part)
    @unpack dx, dy, dt = sysPara
    @unpack pos, R, D = part
    Threads.@threads for i in 3:sysPara.nx-3
        for j in 3:sysPara.nx-3
            # if norm(pos - (SA[i-1, j-1])) < R
            #     du[i, j] = u[i, j]

            # else
            #  !isInside(R, pos, SA[i, j])
            du[i, j] = u[i, j] + dt / 1 * D * ((-u[i+2, j] + 16u[i+1, j] - 30u[i, j] + 16u[i-1, j] - u[i-2, j]) / (12 * dx^2)
                                               +
                                               (-u[i, j+2] + 16u[i, j+1] - 30u[i, j] + 16u[i, j-1] - u[i, j-2]) / (12dy^2))
            # end
        end
    end
end

function updataGridAdvection_4order!(u, du, flow, sysPara, part)
    @unpack dx, dy, dt = sysPara
    @unpack pos, R, D = part
    for _ in 1:1
        Threads.@threads for i in 3:sysPara.nx-3
            for j in 3:sysPara.nx-3
                # if norm(pos - (SA[i-1, j-1])) < R
                #     du[i, j] = u[i, j]

                # else
                #  !isInside(R, pos, SA[i, j])
                du[i, j] = u[i, j] + dt / 1 * D * ((-u[i+2, j] + 16u[i+1, j] - 30u[i, j] + 16u[i-1, j] - u[i-2, j]) / (12 * dx^2)
                                                   -
                                                   ((-u[i+2, j] + 8u[i+1, j] - 8u[i-1, j] + u[i-2, j]) / (12dx)) * flow[i][j][1]
                                                   +
                                                   (-u[i, j+2] + 16u[i, j+1] - 30u[i, j] + 16u[i, j-1] - u[i, j-2]) / (12dy^2)
                                                   -
                                                   ((-u[i, j+2] + 8u[i, j+1] - 8u[i, j-1] + u[i, j-2]) / (12dy)) * flow[i][j][2])
                # end
            end
        end
    end
end


function DirichletBoundary!(du, sysPara)
    @unpack nx, ny, value = sysPara

    @views du[:, 1] .= value
    @views du[:, ny] .= value

    @views du[1, :] .= value
    @views du[nx, :] .= value
    # for i in 1:nx
    #     du[i, 1] = value
    #     du[i, ny] = value
    # end

    # for j in 1:ny
    #     du[1, j] = value
    #     du[ny, j] = value
    # end
end


function diffusion(u, du, flow, sysPara, part)
    for _ in 1:1
        # constSrc!(u, sysPara, part)
        # constSrc2!(u, sysPara, part)
        constFlux2!(u, sysPara, part)
        # constFlux3!(u, sysPara, part)
        # u[250, 250] = 1
        # updataGrid!(u, du, sysPara, part)
        updataGridAdvection!(u, du, flow, sysPara, part)
        # updataGrid4Order!(u, du, sysPara, part)
        # updataGridAdvection_4order!(u, du, flow, sysPara, part)
        DirichletBoundary!(du, sysPara)
        # du[250,250] = 1
        # constSrc!(du, sysPara, part)
        u, du = du, u
    end
    return u, du
end


"""
Return the nearest 4 grid 
"""
function nearGrid(pos)
    ### bottom left corner of the nearest grid
    # @show pos
    p0 = floor.(Int, pos) .+ 1 ### julia array start from 1
    p1 = p0 + SA[1, 0]
    p2 = p0 + SA[1, 1]
    p3 = p0 + SA[0, 1]
    return SA[p0, p1, p2, p3]
end

function getChemForce(field, sysPara, part)
    @unpack pos, R = part
    @unpack dx, dy, npoly = sysPara


    x, y = pos
    θ = range(0, 2π - (2π / npoly), npoly)
    # angle = collect(θ)
    # append!(angle, -collect(θ))
    dθ = θ[2] - θ[1]
    bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ]
    # bound = [SA[x+R*cos(i), y+R*sin(i)] for i in angle[1:end]]
    force = SA[0.0, 0.0]

    # grid = nearGrid(bound[1])
    # i0, j0 = grid[1]
    # refpoint = [field[i, j] for i in i0:i0+1, j in j0:j0+1]
    # interpol = LinearInterpolation((i0:i0+1, j0:j0+1), refpoint)
    # i0 = 
    # j0 = 1g
    # @show bound
    for pt in bound
        x, y = pt
        xlimlo = floor(Int, x / dx - 2) - 1
        xlimup = ceil(Int, x / dx + 2) - 1
        ylimlo = floor(Int, y / dy - 2) - 1
        ylimup = ceil(Int, y / dy + 2) - 1

        refpoint = [field[i, j] for i in xlimlo:xlimup, j in ylimlo:ylimup]
        # interpol = interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell()))))
        xlist = round((xlimlo) * dx, digits=3):dx:round((xlimup) * dx, digits=3)
        ylist = round((ylimlo) * dy, digits=3):dy:round((ylimup) * dy, digits=3)
        sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)

        F = Interpolations.gradient(sitp, x, y)
        vec = pos - pt
        force = force .+ proj(F, vec)
    end

    return force / npoly
end


function getChemForce2(field, sysPara, part)
    @unpack pos, R = part
    @unpack dx, dy, npoly = sysPara

    
    x, y = pos
    θ = range(0, 2π - (2π / npoly), npoly)
    dθ = 2π / npoly
    bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ]
    # bound = [SA[x+R*cos(i), y+R*sin(i)] for i in angle[1:end]]
    force = SA[0.0, 0.0]

    xlimlo = floor(Int, (x - 1.5R) / dx + 1)
    xlimup = ceil(Int, (x + 1.5R) / dx + 1)
    ylimlo = floor(Int, (y - 1.5R) / dy + 1)
    ylimup = ceil(Int, (y + 1.5R) / dy + 1)


    # i0, j0 = grid[1]
    refpoint = [field[i, j] for i in xlimlo:xlimup, j in ylimlo:ylimup]
    # for i in 1:xlimup-xlimlo+1
    #    for j in 1:ylimup-ylimlo+1
    #    refpoint[i,j] = field[xlimlo+i-1, ylimlo+j-1]
    #    end
    # end
    # interpol = interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell()))))
    # interpol = LinearInterpolation((xlimlo:xlimup, ylimlo:ylimup), refpoint)
    # interpol = CubicSplineInterpolation((round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3),
    #                                     round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)),
    #                                     refpoint)

    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)
    # pts = (xlist, ylist)
    # interpol = interpolate(pts, refpoint, Gridded(Linear()))
    # i0 = 
    # j0 = 1g
    # @show bound
    for pt in bound
        # F = Interpolations.gradient(interpol, pt[1], pt[2])
        F = Interpolations.gradient(sitp, pt[1], pt[2])
        vec = pos - pt
        force = force .+ proj(F, vec)
        # force += F
    end

    # return 2*force /npoly
    return force / npoly
end

function getChemForce3(field, sysPara, part)
    @unpack pos, R = part
    @unpack dx, dy, npoly = sysPara


    x, y = pos
    θ = range(0, 2π - (2π / npoly), npoly)
    dθ = 2π / npoly
    bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ]
    # bound = [SA[x+R*cos(i), y+R*sin(i)] for i in angle[1:end]]
    force = SA[0.0, 0.0]

    xlimlo = floor(Int, (x - 2R) / dx + 1)
    xlimup = ceil(Int, (x + 2R) / dx + 1)
    ylimlo = floor(Int, (y - 2R) / dy + 1)
    ylimup = ceil(Int, (y + 2R) / dy + 1)

    # xlimlo = floor(Int, x - R) - 2
    # xlimup = ceil(Int, x + R) + 2
    # ylimlo = floor(Int, y - R) - 2
    # ylimup = ceil(Int, y + R) + 2

    # i0, j0 = grid[1]
    refpoint = [field[i, j] for i in xlimlo:xlimup, j in ylimlo:ylimup]
    interpol = interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell()))))
    # interpol = LinearInterpolation((xlimlo:xlimup, ylimlo:ylimup), refpoint)
    # interpol = CubicSplineInterpolation((round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3),
    #                                     round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)),
    #                                     refpoint)

    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)
    sitp = scale(interpol, xlist, ylist)


    # force = Interpolations.gradient(interpol, x, y)
    force = Interpolations.gradient(sitp, x, y)
    # force = force .+ getNormal(F, pt, pos)
    # force += F
    # end

    return force
end


function getChemForce4(field, sysPara, part)
    @unpack pos, R = part
    @unpack dx, dy, npoly = sysPara


    x, y = pos
    θ = range(0, 2π - (2π / npoly), npoly)

    dθ = θ[2] - θ[1]
    bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ]
    force = SA[0.0, 0.0]


    for pt in bound
        i, j = nearesrPoint(pt, sysPara)
        F = SA[(field[i+1, j]-field[i-1, j])/2dx, (field[i, j+1]-field[i, j-1])/dy]
        vec = pos - pt
        force = force .+ proj(F, vec)
    end

    return force / npoly
end

function getChemForce5(field, sysPara, part)
    @unpack pos, R, = part
    @unpack dx, dy, npoly = sysPara

    
    x, y = pos
    θ = range(0, 2π - (2π / npoly), npoly)
    dθ = 2π / npoly
    bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ]
    # bound = [SA[x+R*cos(i), y+R*sin(i)] for i in angle[1:end]]
    force = SA[0.0, 0.0]

    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)


    # i0, j0 = grid[1]
    refpoint = [field[i, j] for i in xlimlo:xlimup, j in ylimlo:ylimup]

    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)

    len = Int(length(bound)/2)
    for i in 1:len
        pt1 = bound[i] 
        pt2 = bound[len+i]
        F = ( sitp(pt1[1], pt1[2]) - sitp(pt2[1], pt2[2]) ) / 2R
        F *= SA[cos(θ[i]), sin(θ[i])]
        # F = Interpolations.gradient(sitp, pt[1], pt[2])
        # vec = pos - pt
        # force = force .+ proj(F, vec)
        force += F
    end

    # return 2*force /npoly
    return force / npoly
    # return force / npolys
end


function proj(v1, v2)
    v2hat = v2 ./ norm(v2)
    return (v1 ⋅ v2hat) * v2hat
end

function ibm4c(r)
    ϕ = 0
    if r < 1
        ϕ = 1 - 0.5 * r - r^2 + 0.5 * r^3
    elseif 1 <= r < 2
        ϕ = 1 - 11r / 6 + r^2 - (r^3) / 6
    elseif 2 <= r
        ϕ = 0
    end

    return ϕ
end


function getNormal(F, v1, v2)
    vec = v2 - v1
    vhat = vec ./ (norm(vec))
    return abs(F ⋅ vhat) .* vhat
end


function nearesrPoint(pos, sysPara)
    # @unpack dx, dy = sysPara
    x, y = pos

    idx_i = round(Int, x / sysPara.dx) + 1
    idx_j = round(Int, y / sysPara.dy) + 1

    return SA[idx_i, idx_j]
end

function nearestImage(pos0, pos, sysPara)
    @unpack dx, dy, nx, ny = sysPara

    Δx, Δy = pos0 - pos
    Lx, Ly = dx * nx, dy * ny

    if Δx > Lx * 0.5
        Δx = Δx - Lx
        # @show pos0, pos
        # @show Δx, dy

    elseif Δx <= -Lx * 0.5
        Δx = Δx + Ly
        # @show pos0, pos
        # @show Δx, dy
    end

    if Δy > Ly * 0.5
        Δy = Δy - Ly
        # @show pos0, pos
        # @show dx, Δy

    elseif Δy <= -Ly * 0.5
        Δy = Δy + Ly
        # @show pos0, pos
        # @show dx, Δy
    end

    return Δx, Δy
end


function farfield2!(field, sysPara, part)
    @unpack pos, vel, ω0, R = part
    @unpack nx, ny, dx, dy = sysPara
    # field = [[SV(0.0, 0.0) for _ in 1:nx] for _ in 1:ny]
    Threads.@threads for i in 1:nx
        # for i in 1:nx
        for j in 1:ny
            r = norm(pos - (SA[i-1, j-1]) .* SA[dx, dy])
            if r == 0
                field[i][j] = SA[0.0, 0.0]
            elseif r < 1
                pos0 = (SA[i, j] .- 1.0) .* SA[dx, dy]
                # field[i][j] = SA[0.0, 0.0]
                # flow = dipole2D(vel, pos0, pos, part) 
                flow = dipole2D(vel, pos0, pos, part) + rotlet(ω0, pos0, pos, part)
                field[i][j] = flow * r

            else
                pos0 = (SA[i, j] .- 1.0) .* SA[dx, dy]
                # field[i][j] += dipole2D(v, pos, p, para) + rotlet(ω, pos, p, para) + f_dipole(v, pos, p, para)
                field[i][j] = dipole2D(vel, pos0, pos, part) + rotlet(ω0, pos0, pos, part)
                # field[i][j] = dipole2D(vel, pos0, pos, part)
                # field[i][j] = myflow(vel, pos0, pos, part)

                # field[i][j] = SA[0.,0.]
            end
        end
    end
    # end
    # @show field[ii][jj]
    # @show ii, jj
    # field[ii][jj] = SV(0.0, 0.0)
    # field[ii][jj] = v
    # @show field
    # return field
end

function myflow(vel, pos0, pos, part)
    x, y = pos0 - pos
    e0 = atand(vel[2], vel[1])
    r = norm(SA[x, y])

    return -1vel / (2π * r^2)

end

function dipole2D(vel, pos0, pos, part)
    # nx, ny = para.nx, para.ny
    # dx, dy = para.dx, para.dy

    stokesflow = SA[0.0, 0.0]
    D = norm(vel)
    # D = 1
    # @show D
    # D = 0
    e0 = atand(vel[2], vel[1])
    # x, y = nearestImage(pos, pos0, sysPara) ./ part.R
    x, y = pos0 - pos
    # @show x,y
    r = norm(SA[x, y])
    # @show r
    e = atand(y, x)

    vr = D * cosd(e - e0) / (2π * r^2)
    ve = D * sind(e - e0) / (2π * r^2)

    vx = vr * cosd(e) - ve * sind(e)
    vy = vr * sind(e) + ve * cosd(e)
    stokesflow = SA[vx, vy]
    # @show stokesflow
    return stokesflow
end

function rotlet(ω, pos0, pos, part)

    # x, y = nearestImage(pos, pos0, sysPara) ./ part.R
    x, y = pos0 - pos
    r = norm(SA[x, y])

    rvec = SA[x, y, 0.0]
    flow = cross(SA[0.0, 0.0, 1.0], rvec) ./ r^2
    flow = ω .* flow ./ (4π)
    flow = SA[flow[1], flow[2]]
    return flow
end


function setField!(u, sysPara, part)
    f(x, y, θ) = 1x * tand(90 - θ) + 1y

    @unpack pos = part
    @unpack dx, dy, nx, ny = sysPara

    for i in 1:nx
        for j in 1:ny
            u[i, j] = f((i - 1) * dx, (j - 1) * dy, 45)
        end
    end
end

function chemo(u, du, sysPara, part)
    @unpack pos, ϕ, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    dpos = copy(pos)
    dϕ = copy(ϕ)
    all_pos = [dpos for i in 1:Nstep]
    all_F = [SA[0.0, 0.0] for i in 1:Nstep]

    flow_field = [[SV(0.0, 0.0) for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    # setField!(u, sysPara, part)
    for j in 2:Nstep
        # farfield2!(flow_field, sysPara, part) 
        u, du = diffusion(u, du, flow_field, sysPara, part)
        F = getChemForce5(u, sysPara, part)
        # F = SA[0,0]
        all_F[j] = copy(F.* α)
        vel = SA[cos(ϕ), sin(ϕ)] * v0

        part.vel = vel + α * F
        dpos = pos + part.vel * dt
        # dpos = pos + vel * dt 
        # @show all_pos[j]
        all_pos[j] = copy(dpos)
        dϕ = ϕ + ω0 * dt + sqrt(2 * Dr * dt) * randn()

        pos, dpos = dpos, pos
        ϕ, dϕ, = dϕ, ϕ

        part.pos = pos
        part.ϕ = ϕ

    end
    # constSrc!(u, src, pos, R)
    return u, all_pos, all_F, flow_field
end

part = Particle()
sysPara = SysPara(nx=400, ny=400)

sysPara.dt = 0.001
sysPara.dx = 0.2
sysPara.dy = 0.2
# part.pos = (SA[sysPara.nx/3, sysPara.ny/2] .+ randn(2)) .* SA[sysPara.dx, sysPara.dy]
part.pos = (SA[sysPara.nx/2, sysPara.ny/2]) .* SA[sysPara.dx, sysPara.dy]
# part.ϕ = rand()*2π
# part.R = sysPara.dx*6
part.R = 1
# part.v0 = copy(part.R)
part.v0 = 1
part.ω0 = 0.5
part.Dr = 0.0
part.D = 0.1
# part.α =5 * 10^2
part.α = -80
# part.α = 0

sysPara.npoly = 180
sysPara.Nstep = 200_000

u = zeros(sysPara.nx, sysPara.ny)
du = copy(u)

@time u, all_pos, all_F, field = chemo(u, du, sysPara, part)

"""
test diffusion
"""


# function diffuseTest()
#     part.pos = (SA[sysPara.nx/2, sysPara.ny/2]) .* SA[sysPara.dx, sysPara.dy]
#     p0 = copy(part.pos)
#     part.α = 0
#     part.v0 = 1.5
#     part.ω0 = 
#     part.Dr = 0
#     part.D = 1
#     u = zeros(sysPara.nx, sysPara.ny)
#     du = copy(u)
#     @time u, all_pos, all_F, field = chemo(u, du, sysPara, part)

#     testPoints = [SA[xi, 2] for xi in -10:sysPara.dx:20]
#     len = length(testPoints)
#     result = zeros(len)
#     sim = copy(result)
#     for i in 1:len
#         if i % 10 == 0
#             x, y = round.(Int, (p0 .+ testPoints[i])./ sysPara.dx) .+ 1
#             sim[i] = u[x, y]
#         end
#         result[i] = integal(f, testPoints[i], sysPara.dt * sysPara.Nstep, sysPara.dt, r_prime, part.v0)
#         # x, y = Int.((p0 .+ testPoints[i])./ sysPara.dx) .+ 1
#         # comp[i] = u[x, y]
#     end
#     # for i in 1:10*sysPara.dx:len
#     #     x, y = Int.((p0 .+ testPoints[Int(i)])./ sysPara.dx) .+ 1
#     #     comp[i] = u[x, y]
#     # end
#     return result, sim, testPoints, u, all_pos
# end


"""
gradient disk size test
"""
function gradientDistSizeTest(sysPara, part)
    Rlist = 1:6
    dxlist = (1, 0.5, 0.2)
    Flist = []
    for dx in dxlist
        for R in Rlist
            part.R = dx * R
            sysPara.dx = dx

            part.pos = (SA[sysPara.nx/5, sysPara.ny/2] .+ randn(2)) .* SA[sysPara.dx, sysPara.dy]
            u = zeros(sysPara.nx, sysPara.ny)
            du = copy(u)
            println("dx=$(dx) R=$R")
            @time u, all_pos, all_F, field = chemo(u, du, sysPara, part)
            F = norm.(all_F)
            append!(Flist, [F])
        end
    end
    Flist = reshape(Flist, length(Rlist), length(dxlist))

    return Flist, Rlist, dxlist
end


function gradientDistSizeTestPlot(F, Rlist, dxlist)
    pics = Array{Any}(nothing, length(dxlist))
    for i in 1:length(dxlist)
        pic = plot(F[4, i], label="R=$(Rlist[4])dx")
        # pic = plot(F[1,i], label="R=$(Rlist[1])dx")
        for j in 6:2:length(Rlist)
            # dx = Rlist[j]
            # @show dx
            plot!(pic, F[j, i], label="R=$(Rlist[j])dx")
        end
        plot!(pic, xlabel="step", ylabel="force")
        pics[i] = pic
    end
    return pics
end
# F, Rlist, dxlist = gradientDistSizeTest(sysPara, part)







# du, u = diffusion(u, du, sysPara, part)
# constSrc2!(u, sysPara, part)
function viz(u, all_p, sysPara)
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
    plot!(circle(part.R, part.pos[1], part.pos[2]), label="", c=:green)
    plot!(hm, x, y, label="", c=:white)
    return hm
end
circle(R, x, y) = (θ = LinRange(0, 2π, 30);
(x .+ R .* cos.(θ), y .+ R .* sin.(θ)))

function plot_traj(all_pos)
    plot(all_pos)
    plot!(circle(part.R, all_pos[end][1], all_pos[end][2]), label="", c=:red)
end