function savedata(field, all_pos, all_F, flow, part, sysPara)
    if sysPara.flow
        savedir = "raw/v$(part.v0)_w0$(part.ω0)_flow/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    else
        savedir = "raw/v$(part.v0)_w0$(part.ω0)/D$(part.D)_a$(part.α)_dx$(sysPara.dx)_nx$(sysPara.nx)_N$(sysPara.Nstep)"
    end
    data = Dict("field" => field, "pos" => all_pos, "force" => all_F, "flow" => flow, "part" => part, "sysPara" => sysPara)
    save(string(savedir * "\\data.jld2"), data)

    return savedir
end