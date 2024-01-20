
"""
Calculate the chemical force by green function 
"""
#* green function of diffusion equation
function kernal3D(r, r_prime, t, t_prime, D, dim)
    d = dim
    dr = norm(r - r_prime)
    return exp(-(dr)^2 / (4 * D * (t - t_prime))) / (4π * D * (t - t_prime))^(d/2)
end


"""
Kernal of gradient green function

return a static array SA[Fx, Fy] for 2D and SA[Fx, Fy, Fz] for 3D

# Arguments
- 'r': target position
- 'r_prime: postion of droplet at in t_prime 
- 't': current time
- 't_prime': past of time
"""
function gradient_kernal2D(r, r_prime, t, t_prime, D)
    dr = norm(r - r_prime)
    dt = t - t_prime
    return -exp(-(dr)^2 / (4*D*dt)) * (2/(π*(4D*dt)^2)) * (r - r_prime)
end

function gradient_kernal3D(r, r_prime, t, t_prime, D)
    dr = norm(r - r_prime)
    dt = t - t_prime
    return -exp(-(dr)^2 / (4 * D * dt)) * (2 / (π^(3/2) * (4D * dt)^(5/2))) * (r - r_prime)
end


function getChemForceGreen(kernal, sysPara, part::Particle, logger, surface_vec, Step, force_cache)
    @unpack pos, R, D = part
    @unpack dx, dy, nx, ny, npoly, dt = sysPara
    path = logger.pos
    t = dt*Step
    
    for i in eachindex(force_cache)
        force_cache[i] = SA[0.0, 0.0]
    end 
    unit_vec, ϕ, dϕ = surface_vec
    #* find the gradient across a finte size particle
    Threads.@threads for j in eachindex(force_cache)
        # id = Threads.threadid()
        n = unit_vec[j]
        r = n .* R + pos 
        ∇C = force_cache[j]
        for i in 1:Step-1
            ∇C += kernal(r, path[i], t, i * dt, D)*dt
        end
        ∇C = ∇C .- dot(∇C, n) .* n
        force_cache[j] = ∇C * R * dϕ
        # @show force_cache[j]
    end


    return sum(force_cache) / (2π * R)
end

function getChemForceGreen2(kernal, sysPara, part::Particle, logger, surface_vec, Step, force_cache)
    @unpack pos, R, D = part
    @unpack dx, dy, nx, ny, npoly, dt = sysPara
    path = logger.pos
    t = dt * Step

    t_cut = 5000
    ∇C_h = SA[0., 0.] 
    for i in eachindex(force_cache)
        force_cache[i] = SA[0.0, 0.0]
    end
    unit_vec, ϕ, dϕ = surface_vec
    if Step > t_cut
        r = pos
        for i in 1:Step-t_cut
            ∇C_h += kernal(r, path[i], t, i * dt, D) * dt
        end
        for j in eachindex(force_cache)
            n = unit_vec[j]
            r = n .* R + pos
            ∇C = force_cache[j]
            for i in t_cut:Step-1
                ∇C += kernal(r, path[i], t, i * dt, D) * dt
            end
            ∇C += ∇C_h
            ∇C = ∇C .- dot(∇C, n) .* n
            force_cache[j] = ∇C * R * dϕ
        end
    else
        #* find the gradient across a finte size particle
        Threads.@threads for j in eachindex(force_cache)
            # id = Threads.threadid()
            n = unit_vec[j]
            r = n .* R + pos
            ∇C = force_cache[j]
            for i in 1:Step-1
                ∇C += kernal(r, path[i], t, i * dt, D) * dt
            end
            ∇C = ∇C .- dot(∇C, n) .* n
            force_cache[j] = ∇C * R * dϕ
            # @show force_cache[j]
        end
    end

    return sum(force_cache) / (2π * R)
end

function test_green(pos, path, Step, dt, D, kernal, unit_vec, force_cache)
    r = pos
    t = dt*Step
    ∇C_h = SA[0., 0.]
    # for i in 1:Step-1
    @tturbo warn_check_args = false for i in 1:Step-1
        ∇C_h += kernal(r, path[i], t, i * dt, D) * dt
    end

    return ∇C_h
end



function getChemForceGreen(kernal, sysPara, part::Particle3D, logger, surface_vec, Step, force_cache)
    @unpack pos, R, D = part
    @unpack dx, dy, dz, dt,  nx, ny, nz = sysPara
    path = logger.pos
    t = dt * Step

    unit_vec, θ, ϕ, dθ, dϕ = surface_vec

    @Threads.threads for j in eachindex(ϕ)
        for i in eachindex(θ)
            ∇C = force_cache[i,j]
            ∇C = SA[0., 0., 0.]
            n = unit_vec[i, j]
            r = n .* R + pos

            for i in 1:Step-1
                ∇C += kernal(r, path[i], t, i * dt, D) * dt
            end

            ∇C = ∇C .- dot(∇C, n) .* n
            force_cache[i, j] = ∇C * R^2 * sin(θ[i])dθ * dϕ
        end
    end


    return sum(force_cache) / (4π * R^2)
end


function getChemForce_periodic(field, sysPara, part)
    @unpack pos, R, = part
    @unpack dx, dy, nx, ny, npoly = sysPara

    force = SA[0.0, 0.0]

    x, y = pos
    θ = range(0, 2π - (2π / npoly), npoly)
    bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ] #! in physical postion

    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)

    #* values for interpolation 
    # refpoint = [field[i, j] for i in xlimlo:xlimup, j in ylimlo:ylimup]
    refpoint = [field[idxPeriodic(i,nx), idxPeriodic(j,ny)] for i in xlimlo:xlimup, j in ylimlo:ylimup]

    #* change physical position to grid point
    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)

    #* interpolation
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)

    #* find the gradient across a finte size particle
    len = Int(length(bound) / 2)
    for i in 1:len
        pt1 = bound[i]
        pt2 = bound[len+i] #! always opposite to pt1

        F = (sitp(pt1[1], pt1[2]) - sitp(pt2[1], pt2[2])) / 2R #! It is a scaler
        F *= SA[cos(θ[i]), sin(θ[i])] #! times unit vector
        force += F
    end

    return force / len
end

function getChemForce_periodic_static(field, static_field, sysPara, part)
    @unpack pos, R, = part
    @unpack dx, dy, nx, ny, npoly = sysPara
    total_field = field .+ static_field

    force = SA[0.0, 0.0]

    x, y = pos
    θ = range(0, 2π - (2π / npoly), npoly)
    bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ] #! in physical postion

    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)

    #* values for interpolation 
    # refpoint = [field[i, j] for i in xlimlo:xlimup, j in ylimlo:ylimup]
    refpoint = [total_field[idxPeriodic(i, nx), idxPeriodic(j, ny)] for i in xlimlo:xlimup, j in ylimlo:ylimup]

    #* change physical position to grid point
    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)

    #* interpolation
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)

    #* find the gradient across a finte size particle
    len = Int(length(bound) / 2)
    for i in 1:len
        pt1 = bound[i]
        pt2 = bound[len+i] #! always opposite to pt1

        F = (sitp(pt1[1], pt1[2]) - sitp(pt2[1], pt2[2])) / 2R #! It is a scaler
        F *= SA[cos(θ[i]), sin(θ[i])] #! times unit vector
        force += F
    end

    return force / len
end


"""
    Calculate chemotactic force based on the surface gradient of concentration in 2D

    return a static array SA[Fx, Fy] 
"""
function getChemForce2(field, sysPara, part, surface_vec)
    @unpack pos, R, = part
    @unpack dx, dy, nx, ny, npoly = sysPara

    force = SA[0.0, 0.0]

    x, y = pos
    # θ = range(0, 2π - (2π / npoly), npoly)
    # bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ] #! in physical postion

    unit_vec, ϕ, dϕ = surface_vec
    #* indexs
    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)


    #* values for interpolation 
    refpoint = [field[i, j] for i in xlimlo:xlimup, j in ylimlo:ylimup]

    #* change the grid to physical position
    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)

    #* interpolation
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)

    #* find the gradient across a finte size particle
    for j in 1:length(ϕ)
        # unit_vec[i, j] = SA[sin(θ[i])cos(ϕ[j]), sin(θ[i])sin(ϕ[j]), cos(θ[i])]
        n = unit_vec[j]
        pt1 = n .* R
        ∇C = Interpolations.gradient(sitp, x + pt1[1], y + pt1[2])
        ∇C = ∇C .- dot(∇C, n) .* n
        force += ∇C * R * dϕ
    end


    return force / (2π * R)
end


"""
    Calculate chemotactic force based on the surface gradient of concentration in 3D

    return a static array SA[Fx, Fy, Fz] 
"""
# function getChemForce(field, sysPara, part::Particle3D, bound_vec)
#     @unpack pos, R, = part
#     @unpack dx, dy, dz, npoly = sysPara

#     _dx, _dy, _dz = 1/dx, 1/dy, 1/dz
    
#     force = SA[0.0, 0.0, 0.0]

#     x, y, z = pos

#     xlimlo = floor(Int, (x - 1.2R) * _dx + 1)
#     xlimup = ceil(Int, (x + 1.2R) * _dx + 1)
#     ylimlo = floor(Int, (y - 1.2R) * _dy + 1)
#     ylimup = ceil(Int, (y + 1.2R) * _dy + 1)
#     zlimlo = floor(Int, (z - 1.2R) * _dz + 1)
#     zlimup = ceil(Int, (z + 1.2R) * _dz + 1)


#     #* values for interpolation 
#     refpoint = [field[i, j, k] for i in xlimlo:xlimup, j in ylimlo:ylimup, k in zlimlo:zlimup]

#     #* match physical position to grid point
#     xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
#     ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)
#     zlist = round((zlimlo - 1) * dz, digits=3):dz:round((zlimup - 1) * dz, digits=3)
#     #* interpolation
#     sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist, zlist)

#     for i in 1:length(θ)
#         for j in 1:length(ϕ)    
#             # unit_vec[i, j] = SA[sin(θ[i])cos(ϕ[j]), sin(θ[i])sin(ϕ[j]), cos(θ[i])]
#             @inbounds n = unit_vec[i,j]
#             pt1 = n .* R
#             @inbounds ∇C = Interpolations.gradient(sitp, x + pt1[1], y + pt1[2])
#             ∇C = ∇C .- dot(∇C, n) .* n
#             force += ∇C * R * dϕ     
#         end
#     end
#     return force / (length(bound_vec[1]) * 7)
# end

#* get force by surface gradient
function getChemForce2(field, sysPara, part::Particle3D, surface_vec)
    @unpack pos, R, = part
    @unpack dx, dy, dz, nx, ny, nz = sysPara

    unit_vec, θ, ϕ, dθ, dϕ = surface_vec

    _dx, _dy, _dz = 1 / dx, 1 / dy, 1 / dz

    force = SA[0.0, 0.0, 0.0]

    x, y, z = pos

    xlimlo = floor(Int, (x - 1.2R) * _dx + 1)
    xlimup = ceil(Int, (x + 1.2R) * _dx + 1)
    ylimlo = floor(Int, (y - 1.2R) * _dy + 1)
    ylimup = ceil(Int, (y + 1.2R) * _dy + 1)
    zlimlo = floor(Int, (z - 1.2R) * _dz + 1)
    zlimup = ceil(Int, (z + 1.2R) * _dz + 1)


    """ bound check """
    dims = 0
    if xlimup > nx
        println("out of bound in upper limit of x axis")
        dims = 1
    elseif ylimup > ny
        println("out of bound in upper limit of y axis")
        dims = 2
    elseif zlimup > nz
        println("out of bound in upper limit of z axis")
        dims = 3
    elseif xlimlo < 1
        println("out of bound in lower limit of z axis")
        dims = -1
    elseif ylimlo < 1
        println("out of bound in lower limit of y axis")
        dims = -2
    elseif zlimlo < 1
        println("out of bound in lower limit of z axis")
        dims = -3
    end
    """ """

    #* values for interpolation 
    refpoint = [field[i, j, k] for i in xlimlo:xlimup, j in ylimlo:ylimup, k in zlimlo:zlimup]

    #* match physical position to grid point
    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)
    zlist = round((zlimlo - 1) * dz, digits=3):dz:round((zlimup - 1) * dz, digits=3)


    #* interpolation
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist, zlist)

    for j in eachindex(ϕ)
        for i in eachindex(θ)
            # unit_vec[i, j] = SA[sin(θ[i])cos(ϕ[j]), sin(θ[i])sin(ϕ[j]), cos(θ[i])]
            n = unit_vec[i, j]
            pt1 = n .* R
            ∇C = Interpolations.gradient(sitp, x + pt1[1], y + pt1[2], z + pt1[3])
            ∇C = ∇C .- dot(∇C, n) .* n
            force += ∇C * R^2 * sin(θ[i])dθ * dϕ
        end
    end


    return force / (4π*R^2)
end


function WCA_force(all_pos, sysPara, part, i)
    L = (sysPara.nx*sysPara.dx, sysPara.ny*sysPara.dy)
    R = part.R
    cutoff = R*2^(1/6)
    force = SA[0.0, 0.0] 
    for j in eachindex(all_pos)
        if j != i 
            r = distanceVec(all_pos[j], all_pos[i], L) #! form j to i 
            r_norm = norm(r)
            r_hat = r ./ r_norm
            force += WCA(R, r_norm, cutoff) .*r_hat
        end
    end
    return force
end



function soft_force(all_pos, sysPara, part, i)
    L = (sysPara.nx * sysPara.dx, sysPara.ny * sysPara.dy)
    R = part.R
    cutoff =2*R
    force = SA[0.0, 0.0]
    k = 5
    for j in eachindex(all_pos)
        if j != i
            r = distanceVec(all_pos[j], all_pos[i], L) #! form j to i 
            r_norm = norm(r)
            r_hat = r ./ r_norm
            force += softCollision(R, r_norm, k, cutoff) .* r_hat
        end
    end
    return force
end


function getBoundForce(pos, sysPara, part, i)
    force = SA[0., 0.]
    R = part.R
    cutoff = 2*R
    L = (sysPara.nx*sysPara.dx, sysPara.ny*sysPara.dy)
    p = fold(pos, L)
    x , y = p 
    lx, ly = L
    r1 = p - SA[0, y]
    r2 = p - SA[x, 0]
    r3 = p - SA[lx, y]
    r4 = p - SA[x, ly]

    for r in (r1, r2, r3, r4)
        r_norm = norm(r)
        r_hat = r ./ r_norm
        force += boundForce(R, r_norm, cutoff) .* r_hat
    end

    return force
end


function boundForce(r0, r, cutoff)
    F = 0
    k = 20
    if r > cutoff
        return F = 0
    else 
        return k*(2r0 - r)
    end
end
    



function WCA(r0, r, cutoff)
    # cutoff = r0*2^(1/6)
    F = 0
    if r > cutoff
        F =  0
    else
        -48(r0/r)^14 + 24(r0/r)^8
    end

    return F
end

function softCollision(r0, r, k, cutoff)
    if r < cutoff
        return k*(2*r0 - r)
    else 
        return 0
    end
end   




#* normal vector on the surface
function genBoundVec2(part::Particle)
    ϕ = LinRange(0, 2π - (2π/30), 30)
    dϕ = ϕ[2] - ϕ[1]
    unit_vec = Array{SVector{2,Float64},1}(undef, length(ϕ))

    for j in eachindex(ϕ)
        unit_vec[j] = SA[cos(ϕ[j]), sin(ϕ[j])]
    end


    return unit_vec, ϕ, dϕ
end

function genBoundVec2(part::Particle3D)
    θ = LinRange(0, π, 20)
    ϕ = LinRange(0, 2π - (2π/38), 38)
    dϕ = ϕ[2] - ϕ[1]
    dθ = θ[2] - θ[1]

    # vecs = [[SA[0.0, 0.0, 0.0,] for i in 1:npoly] for j in 1:7]
    unit_vec = Array{SVector{3,Float64},2}(undef, length(θ), length(ϕ))
    
    for j in eachindex(ϕ)
        for i in eachindex(θ)
            unit_vec[i, j] = SA[sin(θ[i])cos(ϕ[j]), sin(θ[i])sin(ϕ[j]), cos(θ[i])]
        end
    end

    return unit_vec, θ, ϕ, dθ, dϕ
end