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

