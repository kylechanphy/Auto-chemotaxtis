using Pkg
Pkg.activate(".")
include("../src/AutoChemo.jl")

using JLD2
using DelimitedFiles 


function Simulation_ABP3D_test(sysPara, part::Particle3D, logset)
    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    # if part.pos != SA[0.0, 0.0, 0.0]
    #     @show pos
    # end

    logger = initLogger(part::Particle3D, sysPara)
    


    if logset.savedata == true
        dir = savedir(part, sysPara)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        inputs = [sysPara, part]
        dumpTxt(inputs, sysPara.dir)
    end

    ω_head = SA[cos(ϕ_ω)sin(θ_ω), sin(ϕ_ω)sin(θ_ω), cos(θ_ω)]
    v_head = SA[cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)]
    # chem_field = logger.field
    # dchem_field = copy(chem_field)
    # chem_field = zeros(sysPara.nx, sysPara.ny)
    # dchem_field = copy(chem_field)


    dpos = copy(pos)
    # dω_head = copy(dω_head)
    dv_head = copy(v_head)


    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    # flow_field = logger.flow
    # bound_vec = genBoundVec(sysPara.npoly)
    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep
   
        vel = v_head * v0

        part.vel = vel
        dpos = pos + part.vel * dt


        ωx, ωy, ωz = ω0 * dt * ω_head .+ SA[ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr)]
        torque = RotationVec(ωx, ωy, ωz)

        rot = torque

        dv_head = rot * v_head
        dω_head = rot * ω_head

        dv_head = @fastmath dv_head ./ norm(dv_head)
        dω_head = @fastmath dω_head ./ norm(dω_head)
        
        v_head, dv_head = dv_head, v_head
        ω_head, dω_head = dω_head, ω_head

        pos, dpos = dpos, pos

        part.pos = pos
        part.v = v_head
        logger.pos[j] = copy(pos)
        logger.v[j] = copy(v_head)
        # logger.Fc[j] = copy(F .* α)

        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end

    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end

function Simulation_ABP2D_test(sysPara, part::Particle, logset)
    @unpack pos, ϕ, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    logger = initLogger(part::Particle, sysPara)


    if logset.savedata == true
        dir = savedir(part, sysPara)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        inputs = [sysPara, part]
        dumpTxt(inputs, sysPara.dir)
    end
    # chem_field = logger.field
    # dchem_field = copy(chem_field)
    # chem_field = zeros(sysPara.nx, sysPara.ny)
    # dchem_field = copy(chem_field)


    dpos = copy(pos)
    # dω_head = copy(dω_head)


    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    # flow_field = logger.flow
    # bound_vec = genBoundVec(sysPara.npoly)
    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep
        # flowField!(flow_field, sysPara, part)
        # diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        # F = getChemForce(dchem_field, sysPara, part, bound_vec)
        # chem_field, dchem_field = dchem_field, chem_field

        F = SA[0.0, 0.0]
        vel = SA[cos(ϕ), sin(ϕ)] * v0

        part.vel = vel + α * F
        dpos = pos + part.vel * dt
        # all_pos[j] = copy(dpos)
        dϕ = ϕ + ω0 * dt + sqrt(2 * Dr * dt) * randn()
        pos, dpos = dpos, pos
        ϕ, dϕ, = dϕ, ϕ

        part.pos = pos
        part.ϕ = ϕ

        logger.pos[j] = copy(pos)
        logger.Fc[j] = copy(F .* α)
        # logger.Fc[j] = copy(F .* α)

        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end

    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end

function Simulation_CAP3D_test(sysPara, part::Particle3D, logset)
    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    # if part.pos != SA[0.0, 0.0, 0.0]
    #     @show pos
    # end

    logger = initLogger(part::Particle3D, sysPara)



    if logset.savedata == true
        dir = savedir(part, sysPara)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        inputs = [sysPara, part]
        dumpTxt(inputs, sysPara.dir)
    end

    ω_head = SA[cos(ϕ_ω)sin(θ_ω), sin(ϕ_ω)sin(θ_ω), cos(θ_ω)]
    v_head = SA[cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)]
    # chem_field = logger.field
    # dchem_field = copy(chem_field)
    # chem_field = zeros(sysPara.nx, sysPara.ny)
    # dchem_field = copy(chem_field)


    dpos = copy(pos)
    # dω_head = copy(dω_head)
    dv_head = copy(v_head)

    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep

        vel = v_head * v0

        part.vel = vel
        dpos = pos + part.vel * dt


        ωx, ωy, ωz = ω0 * dt * ω_head .+ SA[ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr)]
        torque = RotationVec(ωx, ωy, ωz)

        rot = torque

        dv_head = rot * v_head
        dω_head = rot * ω_head

        dv_head = @fastmath dv_head ./ norm(dv_head)
        dω_head = @fastmath dω_head ./ norm(dω_head)

        v_head, dv_head = dv_head, v_head
        ω_head, dω_head = dω_head, ω_head

        pos, dpos = dpos, pos

        part.pos = pos
        part.v = v_head
        logger.pos[j] = copy(pos)
        logger.v[j] = copy(v_head)
        # logger.Fc[j] = copy(F .* α)

        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end

    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end

function squard_dis(pos::Vector{SVector{3,Float64}})
    dr = pos .- (pos[1],)
    dr2 = sum.(map(x -> x .^ 2, dr))

    return dr2
end


function squard_disXY(pos::Vector{SVector{3,Float64}})
    dr = pos .- (pos[1],)
    dr2 = map(x -> x .^ 2, dr)
    dr2 = [v[1] + v[2] for v in dr2]

    return dr2
end

function squard_disZ(pos::Vector{SVector{3,Float64}})
    dr = pos .- (pos[1],)
    dr2 = map(x -> x .^ 2, dr)
    dr2 = [v[3] for v in dr2]

    return dr2
end



function MSD2(pos, sysPara)
    dr = [pos[i+1] - pos[i] for i in 1:length(pos)-1]
    dr2 = sum.(map(x -> x .^ 2, dr))
    msd = mean(dr2 / sysPara.dt)

    return msd
end



ξ(dt, Dr) = sqrt(2 * Dr * dt) * randn()
part = Particle3D()
part2d = Particle()
sysPara = SysPara()
logset = LoggerSetting()

# sysPara.nx = 150
# sysPara.ny = 150
# sysPara.nz = 250

# sysPara.dx = 0.2
# sysPara.dy = 0.2
# sysPara.dz = 0.2

""" particle 2d """

part2d.pos = SA[0.0, 0.0]
# part.ϕ0 = 0.0 #* polar angle
# part.θ0 = π/2
# part.ϕ_ω = 0.0
# part.θ_ω = 0.0 

part2d.ϕ = 0 #* polar angle



part2d.R = 1
# part.v0 = copy(part.R)
part2d.v0 = 2
part.ω0 = 0
# lp = 500
part2d.Dr = 1

# part.Pe = 50
# part.D = 1 / part.Pe
# part.α =5 * 10^2
# part.α = -3
part2d.α = 0
""" ----------------------"""

sysPara.dt = 0.0025

# part.pos = (SA[sysPara.nx/2, sysPara.ny/2, sysPara.nz/2]) .* SA[sysPara.dx, sysPara.dy, sysPara.dz]
part.pos = SA[0.0, 0.0, 0.0]
# part.ϕ0 = 0.0 #* polar angle
# part.θ0 = π/2
# part.ϕ_ω = 0.0
# part.θ_ω = 0.0 

part.ϕ = 0 #* polar angle
part.θ = π / 2
part.ϕ_ω = 0
part.θ_ω = 0

part.R = 1
# part.v0 = copy(part.R)
part.v0 = 2
part.ω0 = 0
# lp = 500
part.Dr = 1

# part.Pe = 50
# part.D = 1 / part.Pe
# part.α =5 * 10^2
# part.α = -3
part.α = 0

# sysPara.npoly = 40 #! Must be a multiple of 4
sysPara.Nstep = 20_000



dt_list = [0.01, 0.008, 0.005, 0.003, 0.001]
function compare(dt_list)
    tspan = 100
    D_list = zeros(length(dt_list))
    err_list = copy(D_list)

    for j in 1:length(dt_list)
        dt = dt_list[j]
        sysPara.Nstep = floor(Int, tspan / dt)
        sysPara.dt = dt
        t = 0:sysPara.dt:sysPara.dt*(sysPara.Nstep-1)
        N = 1000
        MSD = zeros(sysPara.Nstep)
        all_slope = zeros(N)
        @time for i in 1:N
            part.pos = SA[0.0, 0.0, 0.0]
            logger = Simulation_ABP3D_test3(sysPara, part, logset)
            MSD = squard_dis(logger.pos)
            slope = t \ MSD
            all_slope[i] = slope
        end

        # MSD = MSD ./ N
        # slope = t \ MSD
        err = std(all_slope) / sqrt(N)
        err = err / 6
        D = mean(all_slope) / 6

        D_list[j] = D
        err_list[j] = err
    end

    return D_list, err_list
end



function Deff_sim2(v0_list)
    N = 100
    tspan = 50
    N_test = length(v0_list)
    all_slope = zeros(N_test)
    dt = sysPara.dt
    sysPara.Nstep = floor(Int, tspan / dt)
    t = 0:sysPara.dt:sysPara.dt*(sysPara.Nstep-1)
    MSD = zeros(sysPara.Nstep)
    logger = initLogger(part, sysPara)

    for j in 1:N_test
        MSD = zeros(sysPara.Nstep)
        # sysPara.Nstep = floor(Int, tspan / dt)
        # t = 0:sysPara.dt:sysPara.dt*(sysPara.Nstep-1)
        @show j
        # part.v0 = v0_list[j]
        @time @sync for i in 1:N
            part.pos = SA[0.0, 0.0, 0.0]
            part.v0 = v0_list[j]
            part.ω0 = 0
            # loggers[tid] 
            logger = Simulation_ABP3D_test(sysPara, part, logset)
            # logger = Simulation_CABP3D_test(sysPara, part, logset)
            MSD += squard_dis(logger.pos)
        end
    
        # MSD = sum(MSD_th)
        MSD = MSD ./ N
        slope = t \ MSD
        all_slope[j] = slope
    end

    return all_slope, logger, t
end





function Deff_sim(Dr_list)
    N = 500
    tspan = 1200
    N_test = length(Dr_list)
    all_slope = zeros(N_test)
    dt = sysPara.dt
    sysPara.Nstep = floor(Int, tspan / dt)

    MSD = zeros(sysPara.Nstep)

    t = 0:sysPara.dt:sysPara.dt*(sysPara.Nstep-1)
    nth = Threads.nthreads()
    tid = Threads.threadid
    # nth = 1
    # tid = 1
    # MSD_th = [zeros(sysPara.Nstep) for i in 1:nth]

    ps = Array{Particle3D,1}(undef, nth)
    loggers = Array{Logger3D,1}(undef, nth)
    for n in 1: nth
        ps[n] = Particle3D()
        # ps[n].ω0 = 0
        ps[n].ϕ = 0 #* polar angle [0, π]
        ps[n].θ = 0 #* [0, 2π]
        ps[n].ϕ_ω = π/2
        ps[n].θ_ω = 0
        ps[n].R = 1

        # loggers[n] = initLogger(ps[n], sysPara)
    end
    for j in 1:N_test
        @show nth, tid()
        MSD_th = [zeros(sysPara.Nstep) for i in 1:nth]
        ps[tid()].Dr = Dr_list[j]
        # sysPara.Nstep = floor(Int, tspan / dt)
        # t = 0:sysPara.dt:sysPara.dt*(sysPara.Nstep-1)
        @show j
        # part.v0 = v0_list[j]
        @time @sync for i in 1:N
            Threads.@spawn begin
                # @show tid
                ps[tid()].Dr = Dr_list[j]
                ps[tid()].pos = SA[0.0, 0.0, 0.0]
                ps[tid()].ω0 = 1
                ps[tid()].v0 = 1
                # loggers[tid] 
                loggers[tid()] = Simulation_ABP3D_test(sysPara, ps[tid()], logset)
                # loggers[tid()] = Simulation_CAP3D_test(sysPara, ps[tid()], logset)
                
                MSD_th[tid()] += squard_dis(loggers[tid()].pos)
                # MSD_th[tid()] += squard_disXY(loggers[tid()].pos)
                # MSD_th[tid()] += squard_disZ(loggers[tid()].pos)
            end
        end
        MSD = sum(MSD_th)
        MSD = MSD ./ N 
        slope = t \ MSD
        all_slope[j] = slope
    end

    return all_slope, loggers, t
end





function Rot(dt, Dr)
    dx, dy, dz = sqrt(2Dr * dt) * @SVector randn(3)
    idM = SMatrix{3,3}(1I)
    θ = sqrt(dx^2 + dy^2 + dz^2)
    θx = SMatrix{3,3}([0 -dz dy
        dz 0 -dx
        -dy dx 0])


    return idM + θx * sin(θ) / θ + θx^(2) * (1 - cos(θ)) / θ^2
end

Deff(v0, Dr, d=3) = (v0^2) / (d * (d - 1) * Dr)

# Chiral_Deff(v0, Dr, w0=1) = (Dr / (Dr^2 + w0^2)) * (v0^2) / 2
# Chiral_Deff(v0, Dr, w0=1) = (Dr / (4*Dr^2 + w0^2)) * (v0^2) * (2/3)
# Chiral_Deff3D(v0, Dr, w0=1) = (1/6)*(v0^2/Dr)*((Dr^2 + w0^2 /12 )/(Dr^2 + w0^2 /4))
# Chiral_Deff3D_lim(v0, Dr) = (1 / 6) * (v0^2 / Dr) ./ 3
Dr_list = [0.1, 0.05, 0.01, 0.005, 0.001]
# @show v0_list
# v0_list = [1,2]
# D_sim2 , t2 = Deff_sim2(v0_list)
D_sim, ps, t = Deff_sim(Dr_list)




# v_range = LinRange(v0_list[1], v0_list[end], 50)
Dr_range = LinRange(Dr_list[1], 0.0005, 50)
# Dr_range = LinRange(0.5, 0.0005, 50)
# D_ex = Chiral_Deff3D.(1, Dr_range)
# D_lim = Chiral_Deff3D_lim.(1, Dr_range)
# D_ex = Chiral_Deff.(v_range, 1, 3)
# plot(Dr_range, D_ex, label="fixed z-axis")
# plot!(Dr_range, D_lim, label="ω/Dr >> 1")
scatter!(Dr_list , D_sim ./ 6, label="simulatiom")
plot!(xlabel="D_r", ylabel="Deff", legend=:topright)



# v_range = LinRange(0, 4, 50)
# D_ex = Chiral_Deff.(v_range, 1)
# plot(v_range, D_ex, label="theory")
# scatter!(D_sim./4, label="simulatiom")
# plot!(xlabel="v_0", ylabel="Deff", legend=:topleft)

