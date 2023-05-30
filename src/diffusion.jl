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

function diffusion!(u, du, sysPara, partSet::Vector{Particle})
    for part in partSet
        constFlux_periodic!(u, sysPara, part)
    end #* point sources with constant rate
    updataGrid!(u, du, sysPara, partSet[1])
    # DirichletBoundary!(du, sysPara)
    NeumannBoundary!(u, du, sysPara)
    # PeriodicBoundary!(du, u, sysPara, partSet[1])
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

function constFlux_periodic!(du, sysPara, part)
    @unpack nx, ny, dx, dy, dt = sysPara
    @unpack pos, R, src = part
    Lx, Ly = (nx - 1) * dx, (ny - 1) * dy
    x, y = pos #! physical positin


    #! grid coordinate, julia array start from 1
    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)

    Threads.@threads for i in xlimlo:xlimup
        for j in ylimlo:ylimup
            # @show i, j, idxPeriodic(i,nx), idxPeriodic(j,ny)
            # i, j = idxperiodic(i,nx), idxperiodic(j, ny)
            δx = distPeriodic(x - idxPeriodic(i - 1, nx) * dx, Lx)
            δy = distPeriodic(y - idxPeriodic(j - 1, ny) * dy, Ly)
            # du[i, j] = du[i, j] + dt * src * ibm4c(abs(x - (i - 1) * dx) / dx) * ibm4c(abs(y - (j - 1) * dy) / dy) / dx^2
            du[idxPeriodic(i, nx), idxPeriodic(j, ny)] = du[idxPeriodic(i, nx), idxPeriodic(j, ny)] + dt * src * ibm4c(δx / dx) * ibm4c(δy / dy) / dx^2
        end
    end
end


function constFlux!(du, sysPara, part::Particle3D)
    @unpack nx, ny, nz, dx, dy, dz, dt = sysPara
    @unpack pos, R, src = part

    _dx, _dy, _dz = 1/dx, 1/dy, 1/dz
    _dx3 = 1/dx^3
    x, y, z = pos #! physical positin


    #! grid coordinate, julia array start from 1
    xlimlo = floor(Int, (x - 1.2R) * _dx + 1)
    xlimup = ceil(Int, (x + 1.2R) * _dx + 1)
    ylimlo = floor(Int, (y - 1.2R) * _dy + 1)
    ylimup = ceil(Int, (y + 1.2R) * _dy + 1)
    zlimlo = floor(Int, (z - 1.2R) * _dz + 1)
    zlimup = ceil(Int, (z + 1.2R) * _dz + 1)

    if xlimup > nx || ylimup > ny || zlimup > nz
        print("constFlux out of bound")
    elseif xlimlo < 0 || ylimlo < 0 ny || zlimlo < 0
        print("constFlux out of bound")
    end
        

    Threads.@threads for k in zlimlo:zlimup
        for j in ylimlo:ylimup
            for i in xlimlo:xlimup
            @inbounds du[i, j, k] = du[i, j, k] + dt * src * ibm4c(abs(x - (i - 1) * dx) * _dx) * ibm4c(abs(y - (j - 1) * dy) * _dy) * ibm4c(abs(z - (k - 1) * dz) * _dz) * _dx3
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
        #* 1/6 = 0.1666666667
        ϕ = 1 - 11r * 0.1666666667 + r^2 - (r^3) * 0.1666666667
    elseif 2 <= r
        ϕ = 0
    end

    return ϕ
end


#= 
* update grid in diffusion equation with 2-order method
=#
function updataGrid!(u, du, sysPara, part)
    @unpack dx, dy, dt, nx, ny = sysPara
    @unpack  D = part
    Threads.@threads for i in 2:nx-1
        for j in 2:ny-1

            @inbounds du[i, j] = u[i, j] + dt * D * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / dx^2
                                           +
                                           (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / dy^2)
        end
    end
end



function updataGrid!(u, du, sysPara, part::Particle3D)
    @unpack dx, dy, dz, dt, nx, ny, nz = sysPara
    @unpack pos, R, D = part

    _dx2, _dy2, _dz2 = 1/dx^2, 1/dy^2, 1/dz^2
    
    Threads.@threads for k in 2:nz-1
         for j in 2:ny-1
            for i in 2:nx-1
                @inbounds du[i, j, k] = u[i, j, k] + dt * D * ((u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) * _dx2
                                            +
                                            (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) * _dy2
                                            +
                                            (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) * _dz2)
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



function PeriodicBoundary!(u, du, sysPara, part)
    @unpack nx, ny, dx, dy, dt = sysPara
    @unpack D = part

    #* updata the top and bottom along x
    for i in 2:nx-1
        du[i, 1] = u[i, 1] + dt * D * ((u[i+1, 1] - 2 * u[i, 1] + u[i-1, 1]) / dx^2
                                    +
                                    (u[i, 2] - 2 * u[i, 1] + u[i, ny]) / dy^2)

        du[i, ny] = u[i, ny] + dt * D * ((u[i+1, ny] - 2 * u[i, ny] + u[i-1, ny]) / dx^2
                                        +
                                        (u[i, 1] - 2 * u[i, ny] + u[i, ny-1]) / dy^2)
    end

    #* updata the left and right along y
    for j in 2:ny-1
            du[1, j] = u[1, j] + dt * D * ((u[2, j] - 2 * u[1, j] + u[nx, j]) / dx^2 
                                           +
                                           (u[1, j+1] - 2 * u[1, j] + u[1, j-1]) / dy^2)

            du[nx, j] = u[nx, j] + dt * D * ((u[1, j] - 2 * u[nx, j] + u[nx-1, j]) / dx^2 
                                            +
                                            (u[nx, j+1] - 2 * u[nx, j] + u[nx, j-1]) / dy^2)                              
    end

    #* updata 4 corner
    du[1,1] = u[1,1] + dt * D *((u[2,1] - 2*u[1,1] + u[nx,1]) / dx^2 + (u[1,2] - 2*u[1,1] + u[1,ny]) / dy^2)
    du[nx,ny] = u[nx,ny] + dt * D *((u[1,ny] - 2*u[nx,ny] + u[nx-1,ny]) / dx^2 + (u[nx,1] - 2*u[nx,ny] + u[nx,ny-1]) / dy^2)
    du[1,ny] = u[1,ny] + dt * D *((u[2,ny] - 2*u[1,ny] + u[nx,ny]) / dx^2 + (u[1,1] - 2*u[1,ny] + u[1,ny-1]) / dy^2)
    du[nx,1] = u[nx,1] + dt * D *((u[1,1] - 2*u[nx,1] + u[nx-1,1]) / dx^2 + (u[nx,2] - 2*u[nx,1] + u[nx,ny]) / dy^2)
end

# function NeumannBoundary!(u, du, sysPara)
#     @unpack nx, ny, = sysPara



#     Threads.@threads for j in 1:ny
#             du[1, j] = (u[3, j] - 4 * u[2, j]) / 3
#             du[nx, j] = (u[nx-2, j] - 4 * u[nx-1, j]) / 3
#     end

#     Threads.@threads for i in 1:nx
#             du[i, 1] = (u[i, 3] - 4 * u[i, 2]) / 3
#             du[i, ny] = (u[i, ny-2] - 4 * u[i, ny-1]) / 3
#     end
# end

function NeumannBoundary!(u, du, sysPara)
    @unpack nx, ny, = sysPara

    Threads.@threads for j in 1:ny
        du[1, j] = u[2, j]
        du[nx, j] = u[nx-1, j]
    end

    Threads.@threads for i in 1:nx
        du[i, 1] = u[i,2]
        du[i, ny] = u[i,ny-1]
    end
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

# function NeumannBoundary!(u::Array{Float64,3}, du, sysPara)
#     @unpack nx, ny, nz = sysPara

#     #* X-Y plane
#     Threads.@threads for i in 1:nx
#         for j in 1:ny
#             du[i, j, 1] = (u[i, j, 3] - 4*u[i,j,2]) / 3
#             du[i, j, nz] = (-u[i, j, nz-2] + 4*u[i, j, nz-1])/3
#         end
#     end
    
#     #* Y-Z plane
#     Threads.@threads for j in 1:ny
#         for k in 1:nz
#             du[1, j, k] = (u[3, j, k] - 4*u[2,j,k])/3
#             du[nx, j, k] = (-u[nx-2, j, k] + 4*u[nx-1, j, k])/3
#         end
#     end

#     #* X-Z plane
#     Threads.@threads for i in 1:nx
#         for k in 1:nz
#             du[i, 1, k] = (u[i, 3, k] - 4*u[i,2,k])/3
#             du[i, ny, k] = (-u[i, ny-2, k] + 4*u[i, ny-1, k])/3
#         end
#     end
# end

function NeumannBoundary!(u::Array{Float64,3}, du, sysPara)
    @unpack nx, ny, nz = sysPara

    #* X-Y plane
    Threads.@threads for j in 1:ny
        for i in 1:nx
           @views du[i, j, 1] = (u[i, j, 3] - 4*u[i,j,2]) / 3
           @views du[i, j, nz] = (-u[i, j, nz-2] + 4 * u[i, j, nz-1]) / 3
        end
    end
    
    #* Y-Z plane
    Threads.@threads for k in 1:nz
        for j in 1:ny
           @views du[1, j, k] = (u[3, j, k] - 4*u[2,j,k])/3
           @views du[nx, j, k] = (-u[nx-2, j, k] + 4*u[nx-1, j, k])/3
        end
    end

    #* X-Z plane
    Threads.@threads for k in 1:nz
        for i in 1:nx
           @views du[i, 1, k] = (u[i, 3, k] - 4*u[i,2,k])/3
           @views du[i, ny, k] = (-u[i, ny-2, k] + 4*u[i, ny-1, k])/3
        end
    end
end