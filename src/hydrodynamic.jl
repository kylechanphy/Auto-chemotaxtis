

function flowField!(field, sysPara, part::Particle)
    @unpack pos, vel, ω0, R = part
    @unpack nx, ny, dx, dy = sysPara

    if sysPara.flow
        Threads.@threads for i in 1:nx
            # @show vel, i,j
            for j in 1:ny
                flow = SA[0.,0.]
                r = norm(pos - (SA[i-1, j-1]) .* SA[dx, dy])
                if r == 0
                    field[i][j] = SA[0.0, 0.0]
                elseif r < 1
                    # pos0 = (SA[i, j] .- 1.0) .* SA[dx, dy]
                    # for flow_tpye in sysPara.flow
                    #     f = Symbol(flow_tpye)
                    #     flow += @eval $f(vel, ω0, pos0, pos, part)
                    # end
                    # field[i][j] = flow * r
                    field[i][j] = SA[0.,0.]

                else
                    pos0 = (SA[i, j] .- 1.0) .* SA[dx, dy]
                    # for flow_tpye in sysPara.flow
                    #     # f = Symbol(flow_tpye)
                    #     # flow += @eval $f(part.vel, part.ω0, pos0, part.pos, part)
                    #     flow = flow + flow_tpye(vel, ω0, pos0, pos, part)
                    # end
                    # field[i][j] = flow
                    # field[i][j] += dipole2D(v, pos, p, para) + rotlet(ω, pos, p, para) + f_dipole(v, pos, p, para)
                    field[i][j] = source(vel, ω0, pos0, pos, part) + rotlet(vel, ω0, pos0, pos, part)
                    # field[i][j] = dipole2D(vel, pos0, pos, part)
                    # field[i][j] = myflow(vel, pos0, pos, part)

                    # field[i][j] = SA[0.,0.]
                end
            end
        end
    end

end


function flowField!(field, sysPara, part::Particle3D)
    @unpack pos, vel, ω0, R = part
    @unpack nx, ny, dx, dy = sysPara
    #* cut off of flow field 
    # cutoff = 20*R
    # field = field

    # if sysPara.flow

        
    #     xlimlo = floor(Int, (x - 1.2R) * _dx + 1)
    #     xlimup = ceil(Int, (x + 1.2R) * _dx + 1)
    #     ylimlo = floor(Int, (y - 1.2R) * _dy + 1)
    #     ylimup = ceil(Int, (y + 1.2R) * _dy + 1)
    #     zlimlo = floor(Int, (z - 1.2R) * _dz + 1)
    #     zlimup = ceil(Int, (z + 1.2R) * _dz + 1)


    #     Threads.@threads for i in 1:nx
    #         # @show vel, i,j
    #         for j in 1:ny
    #             flow = SA[0.,0.]
    #             r = norm(pos - (SA[i-1, j-1]) .* SA[dx, dy])
    #             if r == 0
    #                 field[i][j] = SA[0.0, 0.0]
    #             elseif r < 1
    #                 # pos0 = (SA[i, j] .- 1.0) .* SA[dx, dy]
    #                 # for flow_tpye in sysPara.flow
    #                 #     f = Symbol(flow_tpye)
    #                 #     flow += @eval $f(vel, ω0, pos0, pos, part)
    #                 # end
    #                 # field[i][j] = flow * r
    #                 field[i][j] = SA[0.,0.]

    #             else
    #                 pos0 = (SA[i, j] .- 1.0) .* SA[dx, dy]
    #                 # for flow_tpye in sysPara.flow
    #                 #     # f = Symbol(flow_tpye)
    #                 #     # flow += @eval $f(part.vel, part.ω0, pos0, part.pos, part)
    #                 #     flow = flow + flow_tpye(vel, ω0, pos0, pos, part)
    #                 # end
    #                 # field[i][j] = flow
    #                 # field[i][j] += dipole2D(v, pos, p, para) + rotlet(ω, pos, p, para) + f_dipole(v, pos, p, para)
    #                 field[i][j] = source(vel, ω0, pos0, pos, part) + rotlet(vel, ω0, pos0, pos, part)
    #                 # field[i][j] = dipole2D(vel, pos0, pos, part)
    #                 # field[i][j] = myflow(vel, pos0, pos, part)

    #                 # field[i][j] = SA[0.,0.]
    #             end
    #         end
    #     end
    # end

end


function myflow(vel, w,  pos0, pos, part)
    x, y = pos0 - pos
    e0 = atand(vel[2], vel[1])
    r = norm(SA[x, y])

    return -1vel / (2π * r^2)

end

function source(vel,ω0, pos0, pos, part)
    # nx, ny = para.nx, para.ny
    # dx, dy = para.dx, para.dy

    stokesflow = SA[0.0, 0.0]
    D = norm(vel)
    # D = 1
    # @show D
    # D = 0
    e0 = atand(vel[2], vel[1])
    # x, y = nearestImage(pos, pos0, sysPara) ./ part.R
    x, y = pos0 - pos
    # @show x,y
    r = norm(SA[x, y])
    # @show r
    e = atand(y, x)

    vr = D * cosd(e - e0) / (2π * r^2)
    ve = D * sind(e - e0) / (2π * r^2)

    vx = vr * cosd(e) - ve * sind(e)
    vy = vr * sind(e) + ve * cosd(e)
    stokesflow = SA[vx, vy]
    # @show stokesflow
    return stokesflow
end

function rotlet(vel, ω, pos0, pos, part)

    # x, y = nearestImage(pos, pos0, sysPara) ./ part.R
    x, y = pos0 - pos
    r = norm(SA[x, y])

    rvec = SA[x, y, 0.0]
    flow = cross(SA[0.0, 0.0, 1.0], rvec) ./ r^2
    flow = ω .* flow ./ (4π)
    flow = SA[flow[1], flow[2]]
    return flow
end

