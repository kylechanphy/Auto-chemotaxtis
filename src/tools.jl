
function savedir(part, sysPara)
    if part.Dr != 0
        dir = @sprintf("raw3/Dr%.3f/Pe%s/a%s_dx%.3f_nx%s_N%s",
            part.Dr, part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
    else
        if sysPara.flow
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("raw3/flow/Pe%s/a%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
        else
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("raw3/Pe%s/a%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
        end
    end
    @show dir
    return dir
end


function savedir(part::Particle3D, sysPara)
    if part.Dr != 0
        dir = @sprintf("3D/raw3/Dr%.3f/Pe%s/a%s_dx%.3f_nx%s_N%s",
            part.Dr, part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
    else
        if sysPara.flow
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("3D/raw3/flow/Pe%s/a%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
        else
            # dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
            # dir = @printf("raw2/Pe$(1/part.D)_w$(part.ω0)/D%.4f_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)", part.D)
            dir = @sprintf("3D/raw3/Pe%s/a%s_dx%.3f_nx%s_N%s",
                part.Pe, part.α, sysPara.dx, sysPara.nx, sysPara.Nstep)
        end
    end
    @show dir
    return dir
end
function savedata!(field, all_pos, all_F, flow, part, sysPara)
    # if sysPara.flow
    #     savedir = "raw/v$(part.v0)_w0$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    # else
    #     savedir = "raw/v$(part.v0)_w0$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    # end
    data = Dict("field" => field, "pos" => all_pos, "force" => all_F, "flow" => flow, "part" => part, "sysPara" => sysPara)
    save(string(savedir(part, sysPara) * "/data.jld2"), data)

    # return savedir
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
    flow::Vector{Vector{SVector{2,Float64}}} = [[SA[0.0, 0.0] for _ in 1:2]]
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
    @unpack pos, ϕ0, θ0, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    v_head = SA[cos(ϕ0)sin(θ0), sin(ϕ0)sin(θ0), cos(θ0)]
    chem_field = zeros(sysPara.nx, sysPara.ny, sysPara.nz)

    all_pos = [pos for _ in 1:Nstep]
    all_v = [v_head for _ in 1:Nstep]
    all_F = [SA[0.0, 0.0, 0.0] for _ in 1:Nstep] #* Chemical force

    if sysPara.flow == false
        flow_field = [[SA[0.0, 0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:2]
    else
        flow_field = [[SA[0.0, 0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    end


    logger.pos = all_pos
    logger.v = all_v
    logger.Fc = all_F
    logger.field = chem_field

    #! ummodify flow 
    flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    logger.flow = flow_field


    return logger
end

function dumping(logger, s, part, sysPara, logset)
    if logset.dump_flow == true
        save(string(savedir(part, sysPara) * "/flow/flow_$s.jld2"), Dict("flow" => logger.flow))
    end

    if logset.dump_field == true
        save(string(savedir(part, sysPara) * "/field/field_$(s).jld2"), Dict("field" => logger.field))
    end
end

function dumping(logger::Logger3D, s, part::Particle3D, sysPara, logset)
    if logset.dump_flow == true
        save(string(savedir(part, sysPara) * "3D/flow/flow_$s.jld2"), Dict("flow" => logger.flow))
    end

    if logset.dump_field == true
        save(string(savedir(part, sysPara) * "3D/field/field_$(s).jld2"), Dict("field" => logger.field))
    end
end





#* print all elements in struct
function printStruct(obj)
    T = typeof(obj)
    for name in fieldnames(T)
        println("$name = $(getfield(obj, name))")
    end
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

function FFT(signal, t)
    dt = t[2] - t[1]
    N = length(t)
    freq = rfftfreq(N, 1 / dt)
    F = abs.(rfft(signal))

    return freq, F
end


moving_average(vs, n) = [sum(@view vs[i:(i+n-1)]) / n for i in 1:(length(vs)-(n-1))]
# moving_average(vs, n) = [sum(@view vs[i:(i+n-1)]) / n for i in 1:(length(vs)-(n-1))]