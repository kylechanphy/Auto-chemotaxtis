using Pkg
Pkg.activate(".")
include("../src/AutoChemo.jl")

using StaticArrays
using LinearAlgebra
using Plots


function angle_between(v1::AbstractVector, v2::AbstractVector)
    # Ensure the input vectors are 2D
    if length(v1) != 2 || length(v2) != 2
        error("Both input vectors must be 2D.")
    end

    # Calculate the dot product and the magnitudes of the vectors
    dot_product = dot(v1, v2)
    magnitude_v1 = norm(v1)
    magnitude_v2 = norm(v2)

    # Calculate the cosine of the angle between the vectors
    cosine_angle = dot_product / (magnitude_v1 * magnitude_v2)

    # Calculate the angle in radians and convert to degrees
    angle_radians = acos(cosine_angle)
    # angle_degrees = rad2deg(angle_radians)

    return angle_radians
end

function neutral_flow(pos0, pos, e)
    # flow = SA[0., 0.]
    r = pos - pos0
    if norm(r) > 3
        theta = angle_between(r, e)
        return SA[cos(theta)*-sin(theta)*norm(r)^-3, (cos(theta)^2)*norm(r)^-3,]
    else
        return SA[0., 0.]
    end

end

function source(pos0, pos, e)
    r = pos - pos0
    theta = atan(r[2], r[1])
    return SA[cos(theta), sin(theta)]
end

# function meshgrid(x, y)
#     X = [x for _ in y, x in x]
#     Y = [y for y in y, _ in x]
#     X, Y
# end

x_dim = 100
y_dim = 100

pos0 = SA[50.0, 50.0]
e = SA[0, 1]

flow = Array{SVector{2,Float64},2}(undef, x_dim, y_dim)
for i in 1:x_dim
    for j in 1:y_dim
        pos = SA[i, j]
        flow[i, j] = neutral_flow(pos0, pos, e)
        # flow[i, j] = source(pos0, pos, e)
    end
end

flow_nrom = flow ./ norm.(flow)
x_range = 1:10:x_dim
y_range = 1:10:y_dim

# Extract the u and v components of the vector field
# u_field = map(v -> v[1], flow[1:10:end, 1:10:end]).*10^5
# v_field = map(v -> v[2], flow[1:10:end, 1:10:end]).*10^5

u_field = map(v -> v[1], flow_nrom[1:10:end, 1:10:end]).*10
v_field = map(v -> v[2], flow_nrom[1:10:end, 1:10:end]).*10

# Create the meshgrid for the x and y coordinates
# X = x_range' .* ones(length(y_dim))
# Y = zeros(size(x_range)) .+ y_range
X, Y = meshgrid(x_range, y_range)

# Plot the vector field
quiver(X, Y, quiver=(u_field, v_field), color=:blue, label="", aspect_ratio=1)
# quiver(quiver=(u_field, v_field), color=:blue, label="", aspect_ratio=1)
scatter!([pos0[1]], [pos0[2]])
xlabel!("X")
ylabel!("Y")
title!("2D Vector Field")
