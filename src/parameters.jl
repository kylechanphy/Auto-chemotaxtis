@with_kw mutable struct SysPara
    dx::Float64 = 1
    dy::Float64 = 1
    dt::Float64 = 0.05
    nx::Int64 = 1000
    ny::Int64 = 1000
    Nstep::Int64 = 10

    npoly::Int64 = 64

    value = 0

    flow::Bool = false
end


@with_kw mutable struct Particle
    v0::Float16 = 10
    ω0::Float64 = 1
    Dr::Float64 = 0
    D::Float64 = 1
    R::Float64 = 10
    src::Float64 = 1.0
    α::Float64 = 1.0

    pos::SVector = SA[0.0, 0.0]
    vel::SVector = SA[0.0, 0.0]
    ϕ::Float64 = 0.0
end

@with_kw mutable struct loggerSetting
    savedata::Bool = false
    dump_field::Bool = false
    dump_flow::Bool = false
    every::Int = 10

end
