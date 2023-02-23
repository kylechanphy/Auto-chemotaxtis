function getChemForce(field, sysPara, part)
    @unpack pos, R, = part
    @unpack dx, dy, npoly = sysPara

    force = SA[0.0, 0.0]

    x, y = pos
    θ = range(0, 2π - (2π / npoly), npoly)
    bound = [SA[x+R*cos(i), y+R*sin(i)] for i in θ] #! in physical postion

    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)

    #* values for interpolation 
    refpoint = [field[i, j] for i in xlimlo:xlimup, j in ylimlo:ylimup]

    #* change physical position to grid point
    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)

    #* interpolation
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)

    #* find the gradient across a finte size particle
    len = Int(length(bound)/2)
    for i in 1:len
        pt1 = bound[i] 
        pt2 = bound[len+i] #! always opposite to pt1

        F = (sitp(pt1[1], pt1[2]) - sitp(pt2[1], pt2[2]) ) / 2R #! It is a scaler
        F *= SA[cos(θ[i]), sin(θ[i])] #! times unit vector
        force += F
    end

    return force / len
end

function getChemForce(field, sysPara, part::Particle3D, bound_vec)
    @unpack pos, R, = part
    @unpack dx, dy, dz, npoly = sysPara

    force = SA[0.0, 0.0, 0.0]

    x, y, z = pos

    xlimlo = floor(Int, (x - 1.2R) / dx + 1)
    xlimup = ceil(Int, (x + 1.2R) / dx + 1)
    ylimlo = floor(Int, (y - 1.2R) / dy + 1)
    ylimup = ceil(Int, (y + 1.2R) / dy + 1)
    zlimlo = floor(Int, (z - 1.2R) / dz + 1)
    zlimup = ceil(Int, (z + 1.2R) / dz + 1)

    #* values for interpolation 
    refpoint = [field[i, j, k] for i in xlimlo:xlimup, j in ylimlo:ylimup, k in zlimlo:zlimup]

    #* match physical position to grid point
    xlist = round((xlimlo - 1) * dx, digits=3):dx:round((xlimup - 1) * dx, digits=3)
    ylist = round((ylimlo - 1) * dy, digits=3):dy:round((ylimup - 1) * dy, digits=3)
    zlist = round((zlimlo - 1) * dz, digits=3):dz:round((zlimup - 1) * dz, digits=3)
    #* interpolation
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist, zlist)

    #* find the gradient across a finte size particle
    
    len = Int(length(bound_vec[1]) / 2)
    for j in 1:7    
        for i in 1:len
            pt1 = bound_vec[j][i].*R
            pt2 = bound_vec[j][len+i].*R #! always opposite to pt1

            F = (sitp(x+pt1[1], y+pt1[2], z+pt1[3]) - sitp(x+pt2[1], y+pt2[2], z+pt2[3])) / 2R #! It is a scaler
            F *= bound_vec[j][i] #! times unit vector
            force += F
        end
    end

    return force / (length(bound_vec[1])*7)
end



function genBoundVec(npoly)
    θ = range(0, 2π - (2π / npoly), npoly)
    ϕ = range(π / 4, π - (π / 4), 3)
    vecs = [[SA[0.0, 0.0, 0.0,] for i in 1:npoly] for j in 1:7]
    unit_vec = [SA[cos(i), sin(i), 0] for i in θ] #! in physical postion
    vecs[1] = unit_vec

    RotYM = [RotY(i) for i in ϕ]
    RotXM = [RotX(i) for i in ϕ]
    lenYM = length(RotYM)
    lenXM = length(RotXM)
    
    for i in 1:lenYM
        vecs[i+1] = (RotYM[i],) .* unit_vec
    end

    for i in 1:lenXM
        vecs[i+lenYM+1] = (RotXM[i],) .* unit_vec
    end

    return vecs
end