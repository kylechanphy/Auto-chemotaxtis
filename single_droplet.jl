using Pkg
Pkg.activate(".")
include("src/AutoChemo.jl")

using JLD2


part = Particle()
sysPara = SysPara()


sysPara.nx = 500
sysPara.ny = 500
sysPara.dx = 0.2
sysPara.dy = 0.2

sysPara.dt = 0.001

part.pos = (SA[sysPara.nx/2, sysPara.ny/2]) .* SA[sysPara.dx, sysPara.dy]
# part.ϕ = rand()*2π
# part.R = sysPara.dx*6
part.R = 1
# part.v0 = copy(part.R)
part.v0 = 1
part.ω0 = 1
part.Dr = 0.0
part.D = 0.1
# part.α =5 * 10^2
part.α = -5
# part.α = 0

sysPara.npoly = 180 #! Must be a multiple of 4
sysPara.Nstep = 10_000

sysPara.flow = [source, rotlet]
# sysPara.flow = ["off"]

@time u, all_pos, all_F, flow = Simulation(sysPara, part)

savename = "raw/v$(part.v0)_w0$(part.ω0)_flow/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep).jld2"
# savename="raw/v$(part.v0)_w0$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep).jld2"
data = Dict("field"=>u, "pos"=>all_pos, "force"=>all_F, "flow"=>flow)