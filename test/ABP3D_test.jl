using Pkg
Pkg.activate(".")
include("../src/AutoChemo.jl")

using JLD2

function Simulation_ABP3D_test(sysPara, part::Particle3D, logset)
    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

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
        # flowField!(flow_field, sysPara, part)
        # diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        # F = getChemForce(dchem_field, sysPara, part, bound_vec)
        # chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        vel = v_head * v0

        part.vel = vel 
        dpos = pos + part.vel * dt


        dv_head =  v_head + ω0 * cross(ω_head, v_head) * dt + sqrt(2 * Dr * dt) * cross(randn(3), v_head)
        dv_head = @fastmath dv_head ./ norm(dv_head)
        pos, dpos = dpos, pos
        # ϕ, dϕ, = dϕ, ϕ
        # θ, dθ = dθ, θ

        v_head, dv_head = dv_head, v_head



        part.pos = pos
        part.v = v_head
        # part.ϕ = ϕ
        # part.θ = θ

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

function Simulation_ABP3D_test2(sysPara, part::Particle3D, logset)
    @unpack pos, ϕ0, θ0, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

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
    v_head = SA[cos(ϕ0)sin(θ0), sin(ϕ0)sin(θ0), cos(θ0)]


    dpos = copy(pos)
    # dω_head = copy(dω_head)
    dv_head = copy(v_head)

    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep

        vel = v_head * v0

        part.vel = vel
        dpos = pos + part.vel * dt

        dϕ = ϕ0 + sqrt(2 * Dr * dt)*randn()
        dθ = θ0 + sqrt(2*Dr*dt)*randn()

        v_head = SA[cos(dϕ)sin(dθ), sin(dϕ)sin(dθ), cos(dθ)]

        # dv_head = @fastmath v_head + ω0 * cross(ω_head, v_head) * dt + sqrt(2 * Dr * dt) * cross(randn(3), v_head)
        # dv_head = @fastmath dv_head ./ norm(dv_head)
        pos, dpos = dpos, pos
        ϕ0, dϕ, = dϕ, ϕ0
        θ0, dθ = dθ, θ0

        # v_head, dv_head = dv_head, v_head



        part.pos = pos
        part.v = v_head
        # part.ϕ = ϕ
        # part.θ = θ

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

function Simulation_ABP3D_test3(sysPara, part::Particle3D, logset)
    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

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
    chem_field = logger.field
    dchem_field = copy(chem_field)
    # chem_field = zeros(sysPara.nx, sysPara.ny)
    # dchem_field = copy(chem_field)


    dpos = copy(pos)
    dω_head = copy(ω_head)
    dv_head = copy(v_head)


    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    flow_field = logger.flow
    bound_vec = genBoundVec2(part)
    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep
        # flowField!(flow_field, sysPara, part)
        # diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        # F = getChemForce2(dchem_field, sysPara, part, bound_vec)
        # chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        F = SA[0.0, 0.0, 0.0]
        vel = v_head * v0

        part.vel = vel + α * F
        dpos = pos + part.vel * dt

        
        noise = RotationVec(ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr) )

        ωx, ωy, ωz = ω0*ω_head
        torque = RotationVec(ωx, ωy, ωz)
        rot = noise*torque
        dv_head = rot * v_head
        dω_head = rot * ω_head
        # v_head_update = v_head + dv_head

        # r, dθ, dϕ = CartesianToSpherical(dv_head)
        # dϕ_ω = ϕ_ω + dϕ
        # dθ_ω = θ_ω + dθ
        # ω_head = SA[cos(dϕ_ω)sin(dθ_ω), sin(dϕ_ω)sin(dθ_ω), cos(dθ_ω)]

        # ω_head_updata = ω_head + dω_head


        # dω_head = @fastmath ω_head +  cross(ξ, ω_head)
        dv_head = @fastmath dv_head ./ norm(dv_head)
        dω_head = @fastmath dω_head ./ norm(dω_head)
        pos, dpos = dpos, pos
        # dϕ_ω, ϕ_ω = ϕ_ω, dϕ_ω
        # dθ_ω, θ_ω = θ, dθ_ω

        v_head, dv_head = dv_head, v_head
        ω_head, dω_head = dω_head, ω_head


        part.pos = pos
        # # part.v = v_head
        part.v = ω_head
        # part.ϕ = ϕ
        # part.θ = θ

        logger.pos[j] = copy(pos)
        logger.v[j] = copy(v_head)
        logger.Fc[j] = copy(F .* α)

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

function squard_dis(pos::Vector{SVector{3, Float64}})
    dr = pos .- (pos[1],)
    dr2 = sum.(map(x->x.^2, dr))

    return dr2
end


ξ(dt, Dr) = sqrt(2*Dr*dt)*randn()
part = Particle3D()
sysPara = SysPara()
logset = LoggerSetting()

# sysPara.nx = 150
# sysPara.ny = 150
# sysPara.nz = 250

# sysPara.dx = 0.2
# sysPara.dy = 0.2
# sysPara.dz = 0.2

sysPara.dt = 0.005

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



dt_list = [0.01, 0.008,0.005, 0.003, 0.001]
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
        err = err /6 
        D = mean(all_slope) /6 

        D_list[j] = D
        err_list[j] = err
    end
    
    return D_list, err_list 
end 


function Rot(dt, Dr)
dx, dy, dz =  sqrt(2Dr*dt) * @SVector randn(3)
idM  = SMatrix{3,3}(1I)
θ = sqrt(dx^2 + dy^2 + dz^2)
θx = SMatrix{3,3}([0 -dz dy;
                    dz 0 -dx;
                    -dy dx 0 ])
    

    return idM + θx*sin(θ)/θ + θx^(2) * (1-cos(θ))/θ^2
end