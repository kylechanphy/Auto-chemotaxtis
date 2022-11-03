using Parameters
using Plots
using BenchmarkTools
using StaticArrays
using LinearAlgebra
using Interpolations
using JLD2

include("parameters.jl")
include("tools.jl")
include("intercation.jl")
include("hydrodynamic.jl")
include("diffusion.jl")
include("visualize.jl")



function Simulation(sysPara, part)
    @unpack pos, ϕ, v0, ω0, α, Dr = part
    @unpack dt, Nstep = sysPara

    chem_field = zeros(sysPara.nx, sysPara.ny)
    dchem_field = copy(chem_field)

    dpos = copy(pos)
    dϕ = copy(ϕ)
    
    all_pos = [pos for _ in 1:Nstep]
    all_F = [SA[0.0, 0.0] for _ in 1:Nstep] #* Chemical force
    flow_field = [[SA[0.0, 0.0] for _ in 1:sysPara.nx] for _ in 1:sysPara.ny]
    
    for j in 2:Nstep
        flowField!(flow_field, sysPara, part) 
        diffusion!(chem_field, dchem_field, flow_field, sysPara, part)
        F = getChemForce(dchem_field, sysPara, part)
        chem_field, dchem_field = dchem_field, chem_field

        all_F[j] = copy(F .* α)
        vel = SA[cos(ϕ), sin(ϕ)] * v0

        part.vel = vel + α * F
        dpos = pos + part.vel * dt
        all_pos[j] = copy(dpos)
        dϕ = ϕ + ω0 * dt + sqrt(2 * Dr * dt) * randn()
        pos, dpos = dpos, pos
        ϕ, dϕ, = dϕ, ϕ

        part.pos = pos
        part.ϕ = ϕ

    end
    return chem_field, all_pos, all_F, flow_field
end

