#=
* Solve the chemical field by diffusion
=#
function diffusion!(u, du, flow, sysPara, part)
    constFlux!(u, sysPara, part) #* point sources with constant rate
    updataGridAdvection!(u, du, flow, sysPara, part)
    DirichletBoundary!(du, sysPara)
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
function updataGrid!(u, du, sysPara, part)
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

