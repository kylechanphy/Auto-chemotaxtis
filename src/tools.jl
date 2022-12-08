
function savedir(part, sysPara)
    if sysPara.flow
        dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/flow_D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    else
        dir = "raw2/Pe$(1/part.D)_w$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    end
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

#* Initalse logger 
function initLogger(part, sysPara)
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

function dumping(logger, s, part, sysPara, logset)
    if logset.dump_flow == true
        save(string(savedir(part, sysPara) * "/flow/flow_$s.jld2"),Dict("flow"=>logger.flow))
    end
   
    if logset.dump_field == true
        save(string(savedir(part, sysPara) * "/field/field_$(s).jld2"), Dict("field"=>logger.field))
    end
end


"""
curvature of chemodroplet 
"""

### Auto-correlation function 
function ACF(data, lags)
    lag = [v for v in 0:lags]
    return lag, autocor(data, lag)
end

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
    freq = rfftfreq(N, 1/dt)
    F = abs.(rfft(signal))
    
    return freq, F
end


moving_average(vs, n) = [sum(@view vs[i:(i+n-1)]) / n for i in 1:(length(vs)-(n-1))]