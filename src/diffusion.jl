#=
* Solve the chemical field by diffusion
=#
function diffusion!(u, du, flow, sysPara, part)
    constFlux!(u, sysPara, part) #* point sources with constant rate
    updataGridAdvection!(u, du, flow, sysPara, part)
    DirichletBoundary!(du, sysPara)
    # return u, du
end

function diffusion!(u, du, flow, sysPara, part::Particle3D)
    constFlux!(u, sysPara, part) #* point sources with constant rate
    updataGrid!(u, du, sysPara, part)
    # DirichletBoundary!(du, sysPara)
    NeumannBoundary!(u, du, sysPara)
    # return u, du
end

#=
* implementation of point sources
=#
function constFlux!(du, sysPara, part)
    @unpack nx, ny, dx, dy, dt = sysPara
    @unpack pos, R, src = part

    x, y = pos #! physical positin

    #! grid coordinate, julia array start from 1
    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)


    Threads.@threads for i in xlimlo:xlimup
        for j in ylimlo:ylimup
            du[i, j] = du[i, j] + dt * src * ibm4c(abs(x - (i - 1) * dx) / dx) * ibm4c(abs(y - (j - 1) * dy) / dy) / dx^2
        end
    end
end

function constFlux!(du, sysPara, part::Particle3D)
    @unpack nx, ny, dx, dy, dz, dt = sysPara
    @unpack pos, R, src = part

    x, y, z = pos #! physical positin

    #! grid coordinate, julia array start from 1
    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)
    zlimlo = floor(Int, (z - 1.2R) / dz + 1)
    zlimup = ceil(Int, (z + 1.2R) / dz + 1)


    Threads.@threads for i in xlimlo:xlimup
        for j in ylimlo:ylimup
            for k in zlimlo:zlimup
            du[i, j, k] = du[i, j, k] + dt * src * ibm4c(abs(x - (i - 1) * dx) / dx) * ibm4c(abs(y - (j - 1) * dy) / dy)* ibm4c(abs(y - (k - 1) * dz) / dz) / dx^3
            end
        end
    end
end


#=
* ibm 4-point cubic kernal
! r in distant between grid and particle
=#
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


#= 
* update grid in diffusion equation with 2-order method
=#
function updataGrid!(u, du, sysPara, part::Particle3D)
    @unpack dx, dy, dt, nx, ny = sysPara
    @unpack pos, R, D = part
    Threads.@threads for i in 2:nx-1
        for j in 2:ny-1

            du[i, j] = u[i, j] + dt * D * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / dx^2
                                           +
                                           (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / dy^2)
            # end
        end
    end
end


function updataGrid!(u, du, sysPara, part::Particle3D)
    @unpack dx, dy, dz, dt, nx, ny, nz = sysPara
    @unpack pos, R, D = part
    for i in 2:nx-1
       Threads.@threads for j in 2:ny-1
            for k in 2:nz-1
                du[i, j, k] = u[i, j, k] + dt * D * ((u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) / dx^2
                                            +
                                            (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) / dy^2
                                            +
                                            (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) / dz^2)
            end
        end
    end
end
#= 
* update grid in advecton-diffusion equation with 2-order method
=#
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

#! unchange
function updataGridAdvection!(u, du, flow, sysPara, part::Particle3D)
    @unpack dx, dy, dz, dt, nx, ny, nz = sysPara
    @unpack pos, R, D = part
    for i in 2:nx-1
        Threads.@threads for j in 2:ny-1
            for k in 2:nz-1
                du[i, j, k ] = u[i, j, k] + dt * D * ((u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) / dx^2 - ((u[i+1, j, k] - u[i-1, j, k] / 2dx)) * flow[i][j][k][1]
                                            +
                                            (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) / dy^2 - ((u[i, j+1, k] - u[i, j-1, k] / 2dy)) * flow[i][j][k][2]
                                            +
                                            (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) / dz^2 - ((u[i, j, k+1] - u[i, j, k-1] / 2dz)) * flow[i][j][k][3])
            end
        end
    end
end


#= 
* update grid in diffusion equation with 4-order method
=#
function updataGrid4Order!(u, du, sysPara, part)
    @unpack dx, dy, dt = sysPara
    @unpack pos, R, D = part
    Threads.@threads for i in 3:sysPara.nx-3
        for j in 3:sysPara.nx-3
            du[i, j] = u[i, j] + dt / 1 * D * ((-u[i+2, j] + 16u[i+1, j] - 30u[i, j] + 16u[i-1, j] - u[i-2, j]) / (12 * dx^2)
                                               +
                                               (-u[i, j+2] + 16u[i, j+1] - 30u[i, j] + 16u[i, j-1] - u[i, j-2]) / (12dy^2))
            # end
        end
    end
end


#= 
* update grid in advecton-diffusion equation with 4-order method
=#
function updataGridAdvection_4order!(u, du, flow, sysPara, part)
    @unpack dx, dy, dt = sysPara
    @unpack pos, R, D = part
    for _ in 1:1
        Threads.@threads for i in 3:sysPara.nx-3
            for j in 3:sysPara.nx-3

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
end

function DirichletBoundary!(du::Array{Float64,3}, sysPara)
    @unpack nx, ny, nz, value = sysPara

    @views du[:, 1, :] .= value
    @views du[:, ny, :] .= value

    @views du[1, :, :] .= value
    @views du[nx, :, :] .= value

    @views du[:, :, 1] .= value
    @views du[:, :, nz] .= value
end


function NeumannBoundary!(u::Array{Float64,3}, du, sysPara)
    @unpack nx, ny, nz = sysPara

    #* X-Y plane
    Threads.@threads for i in 1:nx
        for j in 1:ny
            du[i, j, 1] = u[i, j, 3]
            du[i, j, nz] = u[i, j, nz-2]
        end
    end
    
    #* Y-Z plane
    Threads.@threads for j in 1:ny
        for k in 1:nz
            du[1, j, k] = u[3, j, k]
            du[nx, j, k] = u[nx-2, j, k]
        end
    end

    #* X-Z plane
    Threads.@threads for i in 1:nx
        for k in 1:nz
            du[i, 1, k] = u[i, 3, k]
            du[i, ny, k] = u[i, ny-2, k]
        end
    end
end