@with_kw mutable struct SysPara
    dx::Float64 = 1
    dy::Float64 = 1
    dz::Float64 = 1
    dt::Float64 = 0.05
    nx::Int64 = 100
    ny::Int64 = 100
    nz::Int64 = 100

    Nstep::Int64 = 10

    npoly::Int64 = 64

    value = 0

    flow::Bool = false

    dir::String = ""
end


@with_kw mutable struct Particle
    v0::Float16 = 10
    ω0::Float64 = 1
    Dr::Float64 = 0
    Pe::Float64 = 10
    D::Float64 = 1 / Pe
    R::Float64 = 10
    src::Float64 = 1.0
    α::Float64 = -1.0

    pos::SVector = SA[0.0, 0.0]
    vel::SVector = SA[0.0, 0.0]
    ϕ::Float64 = 0.0            #* polar angle [0, 2π]
end


@with_kw mutable struct Particle3D
    v0::Float16 = 10
    ω0::Float64 = 1
    Dr::Float64 = 0
    Pe::Float64 = 10
    D::Float64 = 1 / Pe
    R::Float64 = 10
    src::Float64 = 1.0
    α::Float64 = -1.0

    pos::SVector = SA[0.0, 0.0, 0.0]
    vel::SVector = SA[0.0, 0.0, 0.0]
    v::SVector = SA[0.0, 0.0, 0.0]
    ϕ0::Float64 = 0.0  
    θ0::Float64 = π/2   
    ϕ_ω::Float64 = 0.0          #* polar angle [0, 2π]
    θ_ω::Float64 = 0.0          #* azimuth angle [0, π]
end


@with_kw mutable struct LoggerSetting
    savedata::Bool = false
    dump_field::Bool = false
    dump_flow::Bool = false
    every::Int = 10

end
