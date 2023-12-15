using Pkg
Pkg.activate(".")
include("../src/AutoChemo.jl")

using JLD2

part = Particle3D()
sysPara = SysPara()
logset = LoggerSetting()

sysPara.nx = 400
sysPara.ny = 400
sysPara.nz = 400

sysPara.dx = 0.1
sysPara.dy = 0.1
sysPara.dz = 0.1

sysPara.dt = 0.001

part.pos = (SA[sysPara.nx/2, sysPara.ny/2, sysPara.nz/2]) .* SA[sysPara.dx, sysPara.dy, sysPara.dz]

part.ϕ = 0.0 #* polar angle
part.θ = π/2
part.ϕ_ω = 0.0
part.θ_ω = 0.0 

part.R = 1
# part.v0 = copy(part.R)
part.v0 = 1
part.ω0 = 0
part.Dr = 0

part.Pe = 1
part.D = 1 / part.Pe
# part.α =5 * 10^2
# part.α = -1
part.α = -0

sysPara.npoly = 8 #! Must be a multiple of 4
sysPara.Nstep = 1_000


# sysPara.flow = true

logset.savedata = false

# logset.dump_field = true
# logset.every = 500


r_prime(t_prime, v) = t_prime * v * SA[1.0, 0.0]

function kernal3D(r, r_prime, t, t_prime, D)
    dr = norm(r - r_prime)
    return exp(-(dr)^2 / (4 * D * (t - t_prime))) / (4π * D * (t - t_prime))^(3/2)
end

function integal(f, r, t, dt, r_prime, v)
    Nstep = length(0:dt:t-dt)
    s = 0
    t_prime = 0
    for i in 1:Nstep
        s += f(r, r_prime(t_prime, v), t, t_prime, 1) * dt
        # @show s
        t_prime += dt
    end
    return s
end

function integal2(f, r, t, dt, traj, v)
    Nstep = length(0:dt:t-dt)
    s = 0
    t_prime = 0
    for i in 1:Nstep
        s += f(r, traj[i], t, t_prime, 1) * dt
        # @show s
        t_prime += dt
    end
    return s
end

function diffuseTest()
    part.pos = (SA[sysPara.nx/2, sysPara.ny/2, sysPara.nz/2]) .* SA[sysPara.dx, sysPara.dy, sysPara.dz]
    p0 = copy(part.pos)

    logger = Simulation(sysPara, part, logset)

    # testPoints = [SA[xi, 2] for xi in -10:sysPara.dx:20]
    testPoints = [SA[0, 0, zi] for zi in -5:sysPara.dx:5]
    len = length(testPoints)
    result = zeros(len)
    sim = copy(result)
    traj = copy(logger.pos .- [p0])
    for i in 1:len
        # if i % 10 == 0
        x, y, z = round.(Int, (p0 .+ testPoints[i]) ./ sysPara.dx) .+ 1

        sim[i] = logger.field[x, y, z]
        # end
        # result[i] = integal(f, testPoints[i], sysPara.dt * sysPara.Nstep, sysPara.dt, r_prime, part.v0)
        result[i] = integal2(kernal3D, testPoints[i], sysPara.dt * sysPara.Nstep, sysPara.dt, traj, part.v0)
        # x, y = Int.((p0 .+ testPoints[i])./ sysPara.dx) .+ 1
        # comp[i] = u[x, y]
    end

    return result, sim, logger, testPoints
end
result, sim, logger, testPoints = diffuseTest()



z = [v[3] for v in testPoints]

fig = plot(z, result, label="Theary", framestyle=:box)
# plot!(fig, y, sim, label="Simulation")
scatter!(fig, z[1:2:end], sim[1:2:end], markershape=:circle, c=2, label="Simulation")
plot!(fig, xlabel=L"z", ylabel=L"C",
        # legendfontsize=12,
        # guidefontsize=18,
        # tickfontsize=16,
        fg_legend=:transparent,
        grid=false
        )

fig3D = fig
# savefig(fig, "paper/diffus_test/diff_test3D.svg")
