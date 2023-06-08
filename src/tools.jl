
function savedir(part, sysPara)
    if part.Dr != 0
        dir = @sprintf("Dr%.5f/Pe%s/a%s_dx%.3f_nx%s_N%s",
            part.Dr, part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
    else
        if sysPara.flow
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("flow/Pe%s/a%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
        else
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("Pe%s/a%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
        end
    end
    # @show dir
    dir = logset.prefix * updataFileVersion(dir)
    @show dir
    return dir
end



function savedir(part::Particle3D, sysPara, logset)
    if part.Dr != 0
        dir = @sprintf("Dr%.5f/Pe%s/a%s_dx%.3f_nx%s_N%s",
            part.Dr, part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
    else

        if sysPara.flow
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("flow/Pe%s/a%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
        else
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("Pe%s/a%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
        end

    end

    dir = logset.prefix * updataFileVersion(dir)
    @show dir

    return  dir
end

function savedir(partSet::Vector{Particle}, sysPara)
    part = partSet[1]
    num = length(partSet)
    if part.Dr != 0
        dir = @sprintf("Multiple/raw6/Dr%.5f/Pe%s/a%s_Np%s_dx%.3f_nx%s_N%s",
            part.Dr, part.Pe, part.α, num,sysPara.dx, sysPara.nx, sysPara.Nstep)
    else

        if sysPara.flow
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("Multiple/raw6/flow/Pe%s/a%s_Np%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, num,sysPara.dx, sysPara.nx, sysPara.Nstep)
        else
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("Multiple/raw6/Pe%s/a%s_Np%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, num, sysPara.dx, sysPara.nx, sysPara.Nstep)
        end

    end

    dir = logset.prefix * updataFileVersion(dir)
    @show dir

    return dir
end


function updataFileVersion(dir)
    if ispath(dir)
        # println("isdir")
        str = split(dir, "_")[end]

        if split(str, "")[1] != "v"
            dir = dir * "_v1"
            dir = updataFileVersion(dir)
            # mkpath(dir)
        else
            n = parse(Int64, (filter(isdigit, str)))
            dir = replace(dir, str => "v$(n+1)")
            dir = updataFileVersion(dir)
            # mkpath(dir)
        end
    else
        # mkpath(dir)
    end


    return dir
    # @show dir
end

function savedata!(field, all_pos, all_F, flow, part, sysPara)
    # if sysPara.flow
    #     savedir = "raw/v$(part.v0)_w0$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    # else
    #     savedir = "raw/v$(part.v0)_w0$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    # end
    data = Dict("field" => field, "pos" => all_pos, "force" => all_F, "flow" => flow, "part" => part, "sysPara" => sysPara)
    save((sysPara.dir * "/data.jld2"), data)

    # return savedir
end

function savedata!(logger, partSet::Vector{Particle}, sysPara)
    # if sysPara.flow
    #     savedir = "raw/v$(part.v0)_w0$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    # else
    #     savedir = "raw/v$(part.v0)_w0$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    # end
    data = Dict("field" => logger.field, "pos" => logger.all_pos, "parts"=>partSet, "sysPara" => sysPara)
    save((sysPara.dir * "/data.jld2"), data)

    # return savedir
end


function dumpTxt(obj, dir)
    open(dir * "/input.txt", "w") do io
        # write(io, "vel1=$(mean(norm.(vel)))\n")
        T = typeof(obj)
        for name in fieldnames(T)
            write(io, "$name = $(getfield(obj, name))\n")
        end
    end
end

function dumpTxt(objs::Vector, dir)
    open(dir * "/input.txt", "w") do io
        for obj in objs
            T = typeof(obj)
            for name in fieldnames(T)
                write(io, "$name = $(getfield(obj, name))\n")
            end
        end
    end
end



#=
logging data
=#

#* Define logger
@with_kw mutable struct Logger
    pos::Vector{SVector{2,Float64}} = [SA[0.0, 0.0]]
    Fc::Vector{SVector{2,Float64}} = [SA[0.0, 0.0]]
    field::Matrix{Float64} = zeros(2, 2)
    flow::Vector{Vector{SVector{2,Float64}}} = [[SA[0.0, 0.0] for _ in 1:2]]
end


@with_kw mutable struct Logger3D
    pos::Vector{SVector{3,Float64}} = [SA[0.0, 0.0, 0.0]]
    v::Vector{SVector{3,Float64}} = [SA[0.0, 0.0, 0.0]]
    Fc::Vector{SVector{3,Float64}} = [SA[0.0, 0.0, 0.0]]
    field::Array{Float64,3} = zeros(2, 2, 2)
    flow::Array{SVector{3,Float64},3} = Array{SVector{3,Float64},3}(undef, 2, 2, 2)
end


@with_kw mutable struct LoggerInteracting
    all_pos::Vector{Vector{SVector{2,Float64}}} = Array{Vector{SVector{2, Float64}}}(undef, 2)
    field::Matrix{Float64} = zeros(2, 2)
end
#* Initalse logger 
function initLogger(part::Particle, sysPara)
    logger = Logger()
    @unpack pos, ϕ, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara


    chem_field = zeros(sysPara.nx, sysPara.ny)

    all_pos = [pos for _ in 1:Nstep]
    all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]

    logger.pos = all_pos
    logger.Fc = all_F
    logger.field = chem_field


    flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    logger.flow = flow_field


    return logger
end


#* Initalse logger 3D
function initLogger(part::Particle3D, sysPara)
    logger = Logger3D()
    @unpack pos, ϕ, θ, v0, ω0, α, Dr = part
    @unpack dt, Nstep, nx, ny, nz = sysPara

    v_head = SA[cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)]
    chem_field = zeros(sysPara.nx, sysPara.ny, sysPara.nz)

    all_pos = [pos for _ in 1:Nstep]
    all_v = [v_head for _ in 1:Nstep]
    all_F = [SA[0.0, 0.0, 0.0] for _ in 1:Nstep] #* Chemical force

    if sysPara.flow == false
        flow_field = Array{SVector{3,Float64},3}(undef, 2, 2, 2)
    else
        flow_field = Array{SVector{3,Float64},3}(undef, nx, ny, nz)
    end


    logger.pos = all_pos
    logger.v = all_v
    logger.Fc = all_F
    logger.field = chem_field

    #! ummodify flow 
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    logger.flow = flow_field


    return logger
end

#* Initalse Logger interacting 2D
function initLogger(partSet::Vector{Particle}, sysPara)
    logger = LoggerInteracting()
    @unpack v0, ω0, α, Dr = partSet[1]
    @unpack dt, Nstep = sysPara
    num = length(partSet)

    chem_field = zeros(sysPara.nx, sysPara.ny)

    all_pos = [[SA[0.0, 0.0] for _ in 1:Nstep] for i in 1:num]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    for i in eachindex(partSet)
        all_pos[i][1] = partSet[i].pos
    end

    logger.all_pos = all_pos
    # logger.Fc = all_F
    logger.field = chem_field


    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    # logger.flow = flow_field


    return logger
end

function dumping(logger, s, part, sysPara, logset)
    if logset.dump_flow == true
        save((sysPara.dir * "/flow/flow_$s.jld2"), Dict("flow" => logger.flow))
    end

    if logset.dump_field == true
        save((sysPara.dir * "/field/field_$(s).jld2"), Dict("field" => logger.field))
    end
end

function dumping(logger::Logger3D, s, part::Particle3D, sysPara, logset)
    if logset.dump_flow == true
        save((sysPara.dir * "/flow/flow_$s.jld2"), Dict("flow" => logger.flow))
    end

    if logset.dump_field == true
        save((sysPara.dir * "/field/field_$(s).jld2"), Dict("field" => logger.field))
    end
end

function dumping(logger::LoggerInteracting, s, part, sysPara, logset)
    if logset.dump_flow == true
        save((sysPara.dir * "/flow/flow_$s.jld2"), Dict("flow" => logger.flow))
    end

    if logset.dump_field == true
        save((sysPara.dir * "/field/field_$(s).jld2"), Dict("field" => logger.field))
    end
end

function dumping(static_field, logger::LoggerInteracting, s, part, sysPara, logset)
    if logset.dump_flow == true
        save((sysPara.dir * "/flow/flow_$s.jld2"), Dict("flow" => logger.flow))
    end

    if logset.dump_field == true
        save((sysPara.dir * "/field/field_$(s).jld2"), Dict("field" => logger.field .+ static_field))
    end
end



function idxPeriodic(idx, lim)
    if idx <= 0
        idx = lim + idx
        idx = idxPeriodic(idx, lim)
    elseif idx > lim
        idx = idx - lim
        idx = idxPeriodic(idx, lim)
    else
        return idx
    end
end

function distPeriodic(dist, L)
    dl = abs(dist)
    if dl > L/2
        return distPeriodic(L - dl, L)
        # return L - dl
    else
        return dl
    end
end


#* distance vector form v1 poiting to v2
function distanceVec(v1, v2, L)
    vec1 = fold(v1, L)
    vec2 = fold(v2, L)
    lx, ly = L

    dx, dy = vec2 - vec1

    return SA[dx-round(dx / lx)*lx, dy - round(dy/ly)*ly]
end

function fold(x, L)
    x = x - floor(x / L)* L 
    return x
end

function fold(v::SVector, L::Tuple)
    return  SA[fold(v[1], L[1]),fold(v[2], L[2])]
end

function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))


    return transpose(X), transpose(Y)
end
# x(i) = x(i) - floor(x(i) / x_size) * x_size  ! For a box with the origin at the lower left vertex
# ! Works for x's lying in any image.
# dx = x(j) - x(i)
# dx = dx - nint(dx / x_size) * x_size



#* print all elements in struct
function printStruct(obj)
    T = typeof(obj)
    for name in fieldnames(T)
        println("$name = $(getfield(obj, name))")
    end
end



function CartesianToSpherical(v)
    x,y,z = v
    r = norm(v)
    θ = acos(z/r)
    ϕ = sign(y)*acos(x / sqrt(x^2 + y^2))
    
    return SA[r, θ, ϕ]
end

function SphericalToCartesian(v)
    r, θ, ϕ = v

    return r * SA[cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)]
end

"""
Expanding simulatiom box
"""
function expandBox(u, du, dims, sysPara, part::Particle3D, logger)
    @unpack nx, ny, nz, dx, dy, dz, dt = sysPara
    println("Simulation box expanded in dims=$dims")

    x, y, z = part.pos
    expand = 100
    if dims == 1
        new_u = zeros(nx + expand, ny, nz)
        new_du = similar(new_u)
        @views new_u[1:nx, :, :] = u
        sysPara.nx = nx + expand

        
    elseif dims == -1
        new_u = zeros(nx + expand, ny, nz)
        new_du = similar(new_u)
        @views new_u[expand+1:nx+expand, :, :] = u
        sysPara.nx = nx + expand
        
        logger.pos .+= (SA[dx*expand, 0.0, 0.0],)
        part.pos += SA[dx*expand, 0.0, 0.0]
    elseif dims == 2
        new_u = zeros(nx, ny + expand, nz)
        new_du = similar(new_u)
        @views new_u[:, 1:ny, :] = u
        sysPara.ny = ny + expand
    
    elseif dims == -2
        new_u = zeros(nx, ny + expand, nz)
        new_du = similar(new_u)
        @views new_u[:, expand+1:ny+expand, :] = u
        sysPara.ny = ny + expand
    
        logger.pos .+= (SA[0.0, dy*expand, 0.0],)
        part.pos += SA[0.0, dy*expand, 0.0]

    elseif dims == 3
        new_u = zeros(nx, ny, nz + expand)
        new_du = similar(new_u)
        @views new_u[:, :, 1:nz] = u
        sysPara.nz = nz + expand
      
    elseif dims == -3
        new_u = zeros(nx, ny, nz + expand)
        new_du = similar(new_u)
        @views new_u[:, :, expand+1:nz+expand] = u
        sysPara.nz = nz + expand

        logger.pos .+= (SA[0.0, 0.0, dz*expand], )
        part.pos += SA[0.0, 0.0, dz*expand]
    end
    @show sysPara.nx, sysPara.ny, sysPara.nz

    return new_u, new_du
end


function checkbound(u, du, sysPara, part::Particle3D, logger)
    @unpack nx, ny, nz, dx, dy, dz, dt = sysPara
    @unpack pos, R, src = part

    buffer = 100
    # _dx, _dy, _dz = 1 / dx, 1 / dy, 1 / dz
    # _dx3 = 1 / dx^3
    x, y, z = pos #! physical positin


    #! grid coordinate, julia array start from 1
    xlimlo = floor(Int, (x - R) / dx + 1)
    xlimup = ceil(Int, (x + R) / dx + 1)
    ylimlo = floor(Int, (y - R) / dy + 1)
    ylimup = ceil(Int, (y + R) / dy + 1)
    zlimlo = floor(Int, (z - R) / dz + 1)
    zlimup = ceil(Int, (z + R) / dz + 1)

    dims = 0
    if xlimup > nx - buffer
        println("out of bound")
        dims = 1
    elseif ylimup > ny - buffer
        println("out of bound")
        dims = 2
    elseif zlimup > nz - buffer
        println("out of bound")
        dims = 3
    elseif xlimlo < buffer
        println("out of bound")
        dims = -1
    elseif ylimlo < buffer
        println("out of bound")
        dims = -2
    elseif zlimlo < buffer
        println("out of bound")
        dims = -3
    end

    if dims != 0
        u, du = expandBox(u, du, dims, sysPara, part, logger)
    end

    return u, du
end




"""
curvature of chemodroplet 
"""

## Auto-correlation function 
function ACF(data, lags)
    lag = [v for v in 0:lags]
    return lag, autocor(data, lag)
end


function ACF(vel::Vector{SVector{2,Float64}}, lags)
    lag = [v for v in 0:lags]
    out = zeros(lags + 1)
    ns = length(vel) - lags
    # for k in 1:ns
    zz = dot(view(vel, 1:ns), view(vel, 1:ns)) / ns
    for k in 0:lags
        # for j in 1:ns
        out[k+1] = dot(view(vel, 1:ns), view(vel, 1+k:ns+k))
        # out += [@views (dot(vel[k], vel[k+i]) ./ (norm(vel[k])^2)) for i in 0:lags]
        # end
    end
    out = out / zz / ns
    # end
    return lag, out
end
function ACF(vel::Vector{SVector{3,Float64}}, lags)
    lag = [v for v in 0:lags]
    out = zeros(lags + 1)
    ns = length(vel) - lags
    # for k in 1:ns
    zz = dot(view(vel, 1:ns), view(vel, 1:ns)) / ns
    for k in 0:lags
        # for j in 1:ns
        out[k+1] = dot(view(vel, 1:ns), view(vel, 1+k:ns+k))
        # out += [@views (dot(vel[k], vel[k+i]) ./ (norm(vel[k])^2)) for i in 0:lags]
        # end
    end
    out = out / zz / ns
    # end
    return lag, out
end

# function ACF2(vel::Vector{SVector{2,Float64}}, lags)
#     lag = [v for v in 0:lags]
#     out = zeros(lags + 1)
#     ns = length(vel) - lags
#     # for k in 1:ns
#     # zz = dot(view(vel, 1:ns), view(vel, 1:ns)) / ns
#     for k in 1:ns
#         # for j in 1:ns
#         # out[k+1] = dot(view(vel, 1:ns), view(vel, 1+k:ns+k))
#         out += [@views (dot(vel[k], vel[k+i]) ./ (dot(vel[k], vel[k]))) for i in 0:lags]
#         # end
#     end
#     # out = out / zz / ns
#     # end
#     return lag, out / ns
# end


function FFT(signal, para::Dict)
    N = length(signal)
    dt = para["dt"]
    t = dt:dt:dt*N
    freq, F = FFT(signal, t)

    return freq, F
end

function FFT(signal, para::SysPara)
    N = length(signal)
    dt = para.dt
    t = dt:dt:dt*N
    freq, F = FFT(signal, t)

    return freq, F
end

function FFT(signal, dt::Float64)
    N = length(signal)
    # dt = para.dt
    t = dt:dt:dt*N
    freq, F = FFT(signal, t)

    return freq, F
end
function FFT(signal, t)
    dt = t[2] - t[1]
    N = length(t)
    freq = rfftfreq(N, 1 / dt)
    F = abs.(rfft(signal))

    return freq, F
end




moving_average(vs, n) = [sum(@view vs[i:(i+n-1)]) / n for i in 1:(length(vs)-(n-1))]
# moving_average(vs, n) = [sum(@view vs[i:(i+n-1)]) / n for i in 1:(length(vs)-(n-1))]