using Pkg
Pkg.activate(".")

using LoopVectorization
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
using LaTeXStrings
using DelimitedFiles



include("parameters.jl")
include("tools.jl")
include("ExternalField.jl")
include("intercation.jl")
include("hydrodynamic.jl")
include("diffusion.jl")
include("visualize.jl")



ξ(dt, Dr) = sqrt(2 * Dr * dt) * randn()

"""
Run simulation for a single particle in 2D

# Arguments
- 'sysPara': a struct of systerm parameters
- 'Part:Particle': a struct of single 2D particle
- 'logset': a struct of logger setting 

"""
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
        sysPara.dir = dir
        inputs = [sysPara, part]
        dumpTxt(inputs, sysPara.dir)
    end

    chem_field = logger.field
    dchem_field = copy(chem_field)
    surface_vec = genBoundVec2(part)
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
        # F = getChemForce(dchem_field, sysPara, part)
        F = getChemForce2(dchem_field, sysPara, part, surface_vec)
        chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        heat = SA[cos(ϕ), sin(ϕ)]
        vel = heat * v0

        part.vel = vel + α * F
        dpos = pos + part.vel * dt
        # all_pos[j] = copy(dpos)
        dϕ = ϕ + ω0 * dt + sqrt(2 * Dr * dt) * randn()
        pos, dpos = dpos, pos
        ϕ, dϕ, = dϕ, ϕ

        part.pos = pos
        part.ϕ = ϕ

        logger.pos[j] = copy(pos)
        logger.phi[j] = copy(ϕ)
        logger.Fc[j] = copy(F .* α)

        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end

    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.phi, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end


"""
Run simulation for a single particle in 2D by green function

# Arguments
- 'sysPara': a struct of systerm parameters
- 'Part:Particle': a struct of single 2D particle
- 'logset': a struct of logger setting 

"""
function Simulation_green(sysPara, part::Particle, logset)
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
        sysPara.dir = dir
        inputs = [sysPara, part]
        dumpTxt(inputs, sysPara.dir)
    end

    chem_field = logger.field
    dchem_field = copy(chem_field)
    surface_vec = genBoundVec2(part)
    # chem_field = zeros(sysPara.nx, sysPara.ny)
    # dchem_field = copy(chem_field)


    dpos = copy(pos)
    dϕ = copy(ϕ)

    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    flow_field = logger.flow
    force_cache = [SA[0.0, 0.0] for i in eachindex(surface_vec[1])]
    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep
        # flowField!(flow_field, sysPara, part)
        # diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        # F = getChemForce(dchem_field, sysPara, part)
        # F = getChemForce2(dchem_field, sysPara, part, surface_vec)
        F = getChemForceGreen(gradient_kernal2D, sysPara, part, logger, surface_vec, j, force_cache)
        # chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        heat = SA[cos(ϕ), sin(ϕ)]
        vel = heat * v0

        part.vel = vel + α * F
        dpos = pos + part.vel * dt
        # all_pos[j] = copy(dpos)
        dϕ = ϕ + ω0 * dt + sqrt(2 * Dr * dt) * randn()
        pos, dpos = dpos, pos
        ϕ, dϕ, = dϕ, ϕ

        part.pos = pos
        part.ϕ = ϕ

        logger.pos[j] = copy(pos)
        logger.phi[j] = copy(ϕ)
        logger.Fc[j] = copy(F .* α)

        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end

    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.phi, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end


"""
Run simulation for a single particle in 3D

# Arguments
- 'sysPara': a struct of systerm parameters
- 'Part:Particle': a struct of single 3D particle
- 'logset': a struct of logger setting 

"""
function Simulation(sysPara, part::Particle3D, logset)
    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    logger = initLogger(part::Particle3D, sysPara)


    if logset.savedata == true
        dir = savedir(part, sysPara, logset)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        @show sysPara.dir
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

    logger.pos[1] = copy(pos)
    logger.v[1] = copy(v_head)
    # logger.v[1] = copy(ω_head)
    logger.Fc[1] = copy(SA[0., 0., 0.])

    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    flow_field = logger.flow
    bound_vec = genBoundVec2(part)
    if logset.dump_field || logset.dump_flow
        dumping(logger, 1, part, sysPara, logset)
    end

    if ω0 != 0
        T = 3 * 2π / ω0
        NT = minimum([floor(Int64, T / dt), Nstep])
    else
        NT = 1000
    end

    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2 : Nstep

        #* calucalte chemotactic force
        flowField!(flow_field, sysPara, part)
        # chem_field, dchem_field = diffusion!(chem_field, dchem_field, flow_field, sysPara, part, logger)
        diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        F = getChemForce2(dchem_field, sysPara, part, bound_vec)
        chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        # F = SA[0.0, 0.0, 0.0]
        vel = v_head * v0

        part.vel = vel + α * F
        dpos = part.pos + part.vel * dt


        # noise = RotationVec(ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr))

        ωx, ωy, ωz = ω0 * dt * ω_head .+ SA[ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr)]
        torque = RotationVec(ωx, ωy, ωz)

        rot = torque

        dv_head = rot * v_head
        dω_head = rot * ω_head

        dv_head = @fastmath dv_head ./ norm(dv_head)
        dω_head = @fastmath dω_head ./ norm(dω_head)

        v_head, dv_head = dv_head, v_head
        ω_head, dω_head = dω_head, ω_head

        part.pos, dpos = dpos, part.pos

        part.v = v_head


        # chem_field, dchem_field = checkbound(chem_field, dchem_field, sysPara, part, logger)


        logger.pos[j] = copy(part.pos)
        logger.vhead[j] = copy(v_head)
        logger.dwhead[j] = copy(ω_head)
        # logger.v[j] = copy(ω_head)
        logger.Fc[j] = copy(F .* α)


        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end


    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.vhead, logger.dwhead, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end


#! Benchmark only!!!!
function Simulation(sysPara, part::Particle3D, logset, Noise)
    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    logger = initLogger(part::Particle3D, sysPara)


    if logset.savedata == true
        dir = savedir(part, sysPara, logset)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        @show sysPara.dir
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

    logger.pos[1] = copy(pos)
    logger.v[1] = copy(v_head)
    # logger.v[1] = copy(ω_head)
    logger.Fc[1] = copy(SA[0.0, 0.0, 0.0])

    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    flow_field = logger.flow
    bound_vec = genBoundVec2(part)
    if logset.dump_field || logset.dump_flow
        dumping(logger, 1, part, sysPara, logset)
    end

    if ω0 != 0
        T = 3 * 2π / ω0
        NT = minimum([floor(Int64, T / dt), Nstep])
    else
        NT = 1000
    end

    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep

        #* calucalte chemotactic force
        flowField!(flow_field, sysPara, part)
        # chem_field, dchem_field = diffusion!(chem_field, dchem_field, flow_field, sysPara, part, logger)
        diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        F = getChemForce2(dchem_field, sysPara, part, bound_vec)
        chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        # F = SA[0.0, 0.0, 0.0]
        vel = v_head * v0

        part.vel = vel + α * F
        dpos = part.pos + part.vel * dt


        # noise = RotationVec(ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr))

        # ωx, ωy, ωz = ω0 * dt * ω_head .+ SA[ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr)]
        ωx, ωy, ωz = ω0 * dt * ω_head .+ Noise[j]
        torque = RotationVec(ωx, ωy, ωz)

        rot = torque

        dv_head = rot * v_head
        dω_head = rot * ω_head

        dv_head = @fastmath dv_head ./ norm(dv_head)
        dω_head = @fastmath dω_head ./ norm(dω_head)

        v_head, dv_head = dv_head, v_head
        ω_head, dω_head = dω_head, ω_head

        part.pos, dpos = dpos, part.pos

        part.v = v_head


        # chem_field, dchem_field = checkbound(chem_field, dchem_field, sysPara, part, logger)


        logger.pos[j] = copy(part.pos)
        logger.vhead[j] = copy(v_head)
        logger.dwhead[j] = copy(ω_head)
        # logger.v[j] = copy(ω_head)
        logger.Fc[j] = copy(F .* α)


        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end


    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.vhead, logger.dwhead, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end


"""
Run simulation for a single particle in 3D by green function

# Arguments
- 'sysPara': a struct of systerm parameters
- 'Part:Particle': a struct of single 2D particle
- 'logset': a struct of logger setting 

"""
function Simulation_green(sysPara, part::Particle3D, logset)
    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    logger = initLogger(part::Particle3D, sysPara)


    if logset.savedata == true
        dir = savedir(part, sysPara, logset)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        @show sysPara.dir
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

    logger.pos[1] = copy(pos)
    logger.v[1] = copy(v_head)
    # logger.v[1] = copy(ω_head)
    logger.Fc[1] = copy(SA[0.0, 0.0, 0.0])

    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    flow_field = logger.flow
    bound_vec = genBoundVec2(part)

    force_cache = copy(bound_vec[1])

    if logset.dump_field || logset.dump_flow
        dumping(logger, 1, part, sysPara, logset)
    end

    if ω0 != 0
        T = 3 * 2π / ω0
        NT = minimum([floor(Int64, T / dt), Nstep])
    else
        NT = 1000
    end
 
    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep

        #* calucalte chemotactic force
        # flowField!(flow_field, sysPara, part)
        # chem_field, dchem_field = diffusion!(chem_field, dchem_field, flow_field, sysPara, part, logger)
        # diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        # F = getChemForce2(dchem_field, sysPara, part, bound_vec)
        F = getChemForceGreen(gradient_kernal3D, sysPara, part, logger, bound_vec, j, force_cache)
        # chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        # F = SA[0.0, 0.0, 0.0]
        vel = v_head * v0

        part.vel = vel + α * F
        dpos = part.pos + part.vel * dt


        # noise = RotationVec(ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr))

        ωx, ωy, ωz = ω0 * dt * ω_head .+ SA[ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr)]
        torque = RotationVec(ωx, ωy, ωz)

        rot = torque

        dv_head = rot * v_head
        dω_head = rot * ω_head

        dv_head = @fastmath dv_head ./ norm(dv_head)
        dω_head = @fastmath dω_head ./ norm(dω_head)

        v_head, dv_head = dv_head, v_head
        ω_head, dω_head = dω_head, ω_head

        part.pos, dpos = dpos, part.pos

        part.v = v_head


        # chem_field, dchem_field = checkbound(chem_field, dchem_field, sysPara, part, logger)


        logger.pos[j] = copy(part.pos)
        logger.vhead[j] = copy(v_head)
        logger.dwhead[j] = copy(ω_head)
        # logger.v[j] = copy(ω_head)
        logger.Fc[j] = copy(F .* α)


        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end


    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.vhead, logger.dwhead, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end

#! Benchmark only!!!!
function Simulation_green(sysPara, part::Particle3D, logset, Noise::Vector{SVector{3,Float64}})
    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    logger = initLogger(part::Particle3D, sysPara)


    if logset.savedata == true
        dir = savedir(part, sysPara, logset)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        @show sysPara.dir
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

    logger.pos[1] = copy(pos)
    logger.v[1] = copy(v_head)
    # logger.v[1] = copy(ω_head)
    logger.Fc[1] = copy(SA[0.0, 0.0, 0.0])

    # all_pos = [pos for _ in 1:Nstep]
    # all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    # flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    flow_field = logger.flow
    bound_vec = genBoundVec2(part)

    force_cache = copy(bound_vec[1])

    if logset.dump_field || logset.dump_flow
        dumping(logger, 1, part, sysPara, logset)
    end

    if ω0 != 0
        T = 3 * 2π / ω0
        NT = minimum([floor(Int64, T / dt), Nstep])
    else
        NT = 1000
    end
    # Noise = [SA[ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr)] for _ in 1:Nstep]
    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep

        #* calucalte chemotactic force
        # flowField!(flow_field, sysPara, part)
        # chem_field, dchem_field = diffusion!(chem_field, dchem_field, flow_field, sysPara, part, logger)
        # diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        # F = getChemForce2(dchem_field, sysPara, part, bound_vec)
        F = getChemForceGreen(gradient_kernal3D, sysPara, part, logger, bound_vec, j, force_cache)
        # chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        # F = SA[0.0, 0.0, 0.0]
        vel = v_head * v0

        part.vel = vel + α * F
        dpos = part.pos + part.vel * dt


        # noise = RotationVec(ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr))

        # ωx, ωy, ωz = ω0 * dt * ω_head .+ SA[ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr)]
        ωx, ωy, ωz = ω0 * dt * ω_head .+ Noise[j]
        torque = RotationVec(ωx, ωy, ωz)

        rot = torque

        dv_head = rot * v_head
        dω_head = rot * ω_head

        dv_head = @fastmath dv_head ./ norm(dv_head)
        dω_head = @fastmath dω_head ./ norm(dω_head)

        v_head, dv_head = dv_head, v_head
        ω_head, dω_head = dω_head, ω_head

        part.pos, dpos = dpos, part.pos

        part.v = v_head


        # chem_field, dchem_field = checkbound(chem_field, dchem_field, sysPara, part, logger)


        logger.pos[j] = copy(part.pos)
        logger.vhead[j] = copy(v_head)
        logger.dwhead[j] = copy(ω_head)
        # logger.v[j] = copy(ω_head)
        logger.Fc[j] = copy(F .* α)


        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end


    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.vhead, logger.dwhead, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end


"""
#! developing! 
Run simulation for a multiple particles in 2D

# Arguments
- 'sysPara': a struct of systerm parameters
- 'Part:Particle': a 1D vector of single 2D particle struct
- 'logset': a struct of logger setting 

"""
function Simulation(sysPara, partSet::Vector{Particle}, logset)
    @unpack v0, ω0, α, Dr = partSet[1]
    @unpack dt, Nstep, ny, nx = sysPara

    logger = initLogger(partSet, sysPara)

    if logset.savedata == true
        dir = savedir(partSet, sysPara)
        if ispath(dir)

            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        inputs = [sysPara, partSet[1]]
        dumpTxt(inputs, sysPara.dir)
    end
    # logger.field = sinPlaneWave2D(sysPara)

    chem_field = logger.field
    # chem_field = copy(sinPlaneWave2D(sysPara))
    static_field = zeros(nx, ny)
    dchem_field = copy(chem_field)

    all_pos = [part.pos for part in partSet]
    all_dpos = copy(all_pos)

    all_ϕ = [part.ϕ for part in partSet]
    all_dϕ = copy(all_ϕ)

    all_F = [SA[0.0, 0.0] for _ in 1:length(partSet)]


    dumping(static_field, logger, 1, partSet[1], sysPara, logset)
    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep
        t = dt * j
        diffusion!(chem_field, dchem_field, sysPara, partSet)
        Threads.@threads for i in eachindex(partSet)
            all_F[i] = SA[0.0, 0.0]
            # flowField!(flow_field, sysPara, part)
            # all_F[i] += α * getChemForce_periodic(dchem_field, sysPara, partSet[i])
            # sinPlaneWave2D!(static_field, sysPara, t)
            # all_F[i] += α * getChemForce_periodic_static(dchem_field, static_field, sysPara, partSet[i])
            bound_vec = genBoundVec2(partSet[i])
            all_F[i] += α * getChemForce2(dchem_field, sysPara, partSet[i], bound_vec)
            # all_F[i] += α * getChemForce2(dchem_field,sysPara, partSet[i])
            
            all_F[i] += getBoundForce(all_pos[i],  sysPara, partSet[i], i)
            # all_F[i] += WCA_force(all_pos, sysPara, partSet[i], i)
            # all_F[i] += soft_force(all_pos, sysPara, partSet[i], i)
            


            vel_i = SA[cos(all_ϕ[i]), sin(all_ϕ[i])] * v0

            partSet[i].vel = vel_i + all_F[i]
            all_dpos[i] = all_pos[i] + partSet[i].vel * dt
            # all_pos[j] = copy(dpos)
            all_dϕ[i] = all_ϕ[i] + ω0 * dt + sqrt(2 * Dr * dt) * randn()

            partSet[i].pos = all_pos[i]
            partSet[i].ϕ = all_ϕ[i]

            logger.all_pos[i][j] = copy(all_pos[i])
            # logger.Fc[j] = copy(F .* α)
        end
        chem_field, dchem_field = dchem_field, chem_field
        all_pos, all_dpos = all_dpos, all_pos
        all_ϕ, all_dϕ, = all_dϕ, all_ϕ

        if j % logset.every == 0
            # dumping(logger, j, partSet[1], sysPara, logset)
            dumping(static_field, logger, j, partSet[1], sysPara, logset)
        end

        next!(prog) #* progress bar
    end

    if logset.savedata == true
        savedata!(logger, partSet, sysPara)
    end
    logger.field = logger.field + static_field
    # return chem_field, all_pos, all_F, flow_field, logger
    return logger
end



function Simulation(sysPara, part::Particle3D, logset, old_data::Dict)
    @unpack dt, Nstep = sysPara

    @unpack pos, ϕ, θ, ϕ_ω, θ_ω, v0, ω0, α, Dr = part

    logger = initLogger(sysPara, part::Particle3D, old_data)
    part.pos = copy(logger.pos[1])

    if logset.savedata == true
        dir = savedir(part, sysPara, logset)
        if ispath(dir)
            nothing
        else
            mkpath(dir)
        end
        sysPara.dir = dir
        @show sysPara.dir
        inputs = [sysPara, part]
        dumpTxt(inputs, sysPara.dir)
    end

    pos = copy(logger.pos[1])
    ω_head = copy(logger.dwhead[1])
    v_head = copy(logger.vhead[1])

    chem_field = logger.field
    dchem_field = copy(chem_field)

    
    dpos = copy(pos)
    dω_head = copy(ω_head)
    dv_head = copy(v_head)


    flow_field = logger.flow
    bound_vec = genBoundVec2(part) #* discretise surface of sphere

    if logset.dump_field || logset.dump_flow
        dumping(logger, 1, part, sysPara, logset)
    end

    if ω0 != 0
        T = 3 * 2π / ω0
        NT = minimum([floor(Int64, T / dt), Nstep])
    else
        NT = 1000
    end

    prog = Progress(Nstep - 1, 5) #* progress bar
    for j in 2:Nstep

        #* calucalte chemotactic force
        flowField!(flow_field, sysPara, part)
        # chem_field, dchem_field = diffusion!(chem_field, dchem_field, flow_field, sysPara, part, logger)
        diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        F = getChemForce2(dchem_field, sysPara, part, bound_vec)
        chem_field, dchem_field = dchem_field, chem_field

        # all_F[j] = copy(F .* α)
        # F = SA[0.0, 0.0, 0.0]
        vel = v_head * v0

        part.vel = vel + α * F
        dpos = part.pos + part.vel * dt


        # noise = RotationVec(ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr))
        if sysPara.perturbation == true && j <= NT
            ωx, ωy, ωz = ω0*dt * ω_head .+  SA[ξ(dt, 1/500), ξ(dt, 1/500), ξ(dt, 1/500)]
        else
            ωx, ωy, ωz = ω0 * dt * ω_head .+ SA[ξ(dt, Dr), ξ(dt, Dr), ξ(dt, Dr)]
        end

        torque = RotationVec(ωx, ωy, ωz)
        rot = torque

        dv_head = rot * v_head
        dω_head = rot * ω_head

        dv_head = @fastmath dv_head ./ norm(dv_head)
        dω_head = @fastmath dω_head ./ norm(dω_head)
        
        v_head, dv_head = dv_head, v_head
        ω_head, dω_head = dω_head, ω_head

        part.pos, dpos = dpos, part.pos

        part.v = v_head

        #* extand the simulation box
        # chem_field, dchem_field = checkbound(chem_field, dchem_field, sysPara, part, logger)


        logger.pos[j] = copy(part.pos)
        logger.vhead[j] = copy(v_head)
        logger.dwhead[j] = copy(ω_head)
        # logger.v[j] = copy(ω_head)
        logger.Fc[j] = copy(F .* α)


        if j % logset.every == 0
            dumping(logger, j, part, sysPara, logset)
        end

        next!(prog) #* progress bar
    end


    if logset.savedata == true
        savedata!(logger.field, logger.pos, logger.vhead, logger.dwhead, logger.Fc, logger.flow, part, sysPara)
    end
    # return chem_field, all_pos, all_F, flow_field, logger

    return logger
end
    



# function randVec()
#     v = randn(3)
#     # return @fastmath v ./ norm(v)
#     return v
# end 

