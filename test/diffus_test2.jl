using Pkg
Pkg.activate(".")
include("../src/AutoChemo.jl")

using JLD2
part = Particle()
sysPara = SysPara()
logset = LoggerSetting()

sysPara.nx = 500
sysPara.ny = 500


sysPara.dx = 0.2
sysPara.dy = 0.2


sysPara.dt = 0.01

part.pos = (SA[sysPara.nx/2, sysPara.ny/2]) .* SA[sysPara.dx, sysPara.dy]
p0 = copy(part.pos)
# part.pos = (SA[2*sysPara.nx/3, 2*sysPara.ny/3]) .* SA[sysPara.dx, sysPara.dy]
# part.ϕ = rand()*2π
# part.R = sysPara.dx*6
part.R = 1
# part.v0 = copy(part.R)
part.v0 = 1
part.ω0 = 0
part.Dr = 0

part.Pe = 1
part.D = 1 / part.Pe
# part.α =5 * 10^2
part.α = -0
# part.α = 0

sysPara.npoly = 180 #! Must be a multiple of 4
sysPara.Nstep = 1000



logset.savedata = false

r_prime(t_prime, v) = t_prime * v * SA[1.0, 0.0]

function f(r, r_prime, t, t_prime, D)
    dr = norm(r - r_prime)
    return exp(-(dr)^2 / (4 * D * (t - t_prime))) / (4π * D * (t - t_prime))
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
    part.pos = (SA[sysPara.nx/2, sysPara.ny/2]) .* SA[sysPara.dx, sysPara.dy]
    p0 = copy(part.pos)

    logger = Simulation(sysPara, part, logset)

    # testPoints = [SA[xi, 2] for xi in -10:sysPara.dx:20]
    testPoints = [SA[0, xi] for xi in -10:sysPara.dx:10]
    len = length(testPoints)
    result = zeros(len)
    sim = copy(result)
    traj = copy(logger.pos .- [p0])
    for i in 1:len
        # if i % 10 == 0
        x, y = round.(Int, (p0 .+ testPoints[i]) ./ sysPara.dx) .+ 1
        sim[i] = logger.field[x, y]
        # end
        # result[i] = integal(f, testPoints[i], sysPara.dt * sysPara.Nstep, sysPara.dt, r_prime, part.v0)
        result[i] = integal2(f, testPoints[i], sysPara.dt * sysPara.Nstep, sysPara.dt, traj, part.v0)
        # x, y = Int.((p0 .+ testPoints[i])./ sysPara.dx) .+ 1
        # comp[i] = u[x, y]
    end
    # for i in 1:10*sysPara.dx:len
    #     x, y = Int.((p0 .+ testPoints[Int(i)])./ sysPara.dx) .+ 1
    #     comp[i] = u[x, y]
    # end
    return result, sim, testPoints
end
result, sim, testPoints = diffuseTest()


y = [v[2] for v in testPoints]

fig = plot(y, result, label="Theary", framestyle=:box)
# plot!(fig, y, sim, label="Simulation")
scatter!(fig, y[1:2:end], sim[1:2:end], markershape=:circle, c=2, label="Simulation")
plot!(fig, xlabel=L"y", ylabel=L"C",
        legendfontsize=14,
        guidefontsize=18,
        tickfontsize=18,
        fg_legend=:transparent,
        grid=false
        )

# savefig("paper/diffus_test/diff_test2D.svg")
