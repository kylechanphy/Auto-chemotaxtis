using Parameters
using Plots
using BenchmarkTools
using StaticArrays
using LinearAlgebra
using Interpolations
using JLD2
using ProgressMeter
using StatsBase
using FFTW
using Printf
using Rotations
using Measures
# using GLMakie
# using Revise



include("parameters.jl")
include("tools.jl")
include("intercation.jl")
include("hydrodynamic.jl")
include("diffusion.jl")
include("visualize.jl")



function Simulation(sysPara, part::Particle, logset)
    @unpack pos, ϕ, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    logger = initLogger(part, sysPara)

    if logset.savedata == true
        dir = savedir(part, sysPara)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        inputs = [sysPara, part]
        dumpTxt(inputs, dir)
    end

    chem_field = logger.field
    dchem_field = copy(chem_field)
    # chem_field = zeros(sysPara.nx, sysPara.ny)
    # dchem_field = copy(chem_field)


    dpos = copy(pos)
    dϕ = copy(ϕ)

    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    flow_field = logger.flow

    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep
        flowField!(flow_field, sysPara, part)
        diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        F = getChemForce(dchem_field, sysPara, part)
        chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
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



function Simulation(sysPara, part::Particle3D, logset)
    @unpack pos, ϕ0, θ0, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    logger = initLogger(part::Particle3D, sysPara)

    if logset.savedata == true
        #* get save path
        dir = savedir(part, sysPara)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        inputs = [sysPara, part]
        dumpTxt(inputs, dir)
    end

    ω_head = SA[cos(ϕ_ω)sin(θ_ω), sin(ϕ_ω)sin(θ_ω), cos(θ_ω)]
    v_head = SA[cos(ϕ0)sin(θ0), sin(ϕ0)sin(θ0), cos(θ0)]
    chem_field = logger.field
    dchem_field = copy(chem_field)
    # chem_field = zeros(sysPara.nx, sysPara.ny)
    # dchem_field = copy(chem_field)


    dpos = copy(pos)
    # dω_head = copy(dω_head)
    dv_head = copy(v_head)


    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    flow_field = logger.flow
    bound_vec = genBoundVec(sysPara.npoly)
    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep
        flowField!(flow_field, sysPara, part)
        diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        F = getChemForce(dchem_field, sysPara, part, bound_vec)
        chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        vel = v_head * v0

        part.vel = vel + α * F
        dpos = pos + part.vel * dt
        # all_pos[j] = copy(dpos)
        # dϕ = ϕ + ω0 * dt + sqrt(2 * Dr * dt) * randn()
        # dθ = θ + sqrt(2 * Dr * dt) * randn()
        # dϕ = ϕ + ω0*sin(π/2)*sin(π/2 - θ)*dt + Dr/tan(ϕ)*dt + sqrt(2 * Dr * dt) * randn()
        # dθ = θ + ω0*(cos(π/2) - sin(π/2)cos(π/2 - θ)/tan(ϕ))*dt + sqrt(2 * Dr * dt)/sin(ϕ) * randn()
        
        dv_head = @fastmath v_head +  ω0*cross(ω_head, v_head)*dt + sqrt(2 * Dr * dt)*cross(randVec(), v_head)
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


function randVec()
    v = randn(3)
    return @fastmath v ./ norm(v)
end