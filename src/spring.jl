"""
    Use the spring/repulsion model of Fruchterman and Reingold (1991):
        Attractive force:  f_a(d) =  d^2 / k
        Repulsive force:  f_r(d) = -k^2 / d
    where d is distance between two vertices and the optimal distance
    between vertices k is defined as C * sqrt( area / num_vertices )
    where C is a parameter we can adjust

    Arguments:
    adj_matrix Adjacency matrix of some type. Non-zero of the eltype
               of the matrix is used to determine if a link exists,
               but currently no sense of magnitude
    C          Constant to fiddle with density of resulting layout
    MAXITER    Number of iterations we apply the forces
    INITTEMP   Initial "temperature", controls movement per iteration
"""
function layout_spring_adj{T}(adj_matrix::Array{T,2}; C=2.0, MAXITER=100, INITTEMP=2.0)

    size(adj_matrix, 1) != size(adj_matrix, 2) && error("Adj. matrix must be square.")
    const N = size(adj_matrix, 1)

    # Initial layout is random on the square [-1,+1]^2
    locs_x = 2*rand(N) .- 1.0
    locs_y = 2*rand(N) .- 1.0

    # The optimal distance bewteen vertices
    const K = C * sqrt(4.0 / N)

    # Store forces and apply at end of iteration all at once
    force_x = zeros(N)
    force_y = zeros(N)

    # Iterate MAXITER times
    @inbounds for iter = 1:MAXITER
        # Calculate forces
        for i = 1:N
            force_vec_x = 0.0
            force_vec_y = 0.0
            for j = 1:N
                i == j && continue
                d_x = locs_x[j] - locs_x[i]
                d_y = locs_y[j] - locs_y[i]
                d   = sqrt(d_x^2 + d_y^2)
                if adj_matrix[i,j] != zero(eltype(adj_matrix)) || adj_matrix[j,i] != zero(eltype(adj_matrix))
                    # F = d^2 / K - K^2 / d
                    F_d = d / K - K^2 / d^2
                else
                    # Just repulsive
                    # F = -K^2 / d^
                    F_d = -K^2 / d^2
                end
                # d  /          sin θ = d_y/d = fy/F
                # F /| dy fy    -> fy = F*d_y/d
                #  / |          cos θ = d_x/d = fx/F
                # /---          -> fx = F*d_x/d
                # dx fx
                force_vec_x += F_d*d_x
                force_vec_y += F_d*d_y
            end
            force_x[i] = force_vec_x
            force_y[i] = force_vec_y
        end
        # Cool down
        TEMP = INITTEMP / iter
        # Now apply them, but limit to temperature
        for i = 1:N
            force_mag  = sqrt(force_x[i]^2 + force_y[i]^2)
            scale      = min(force_mag, TEMP)/force_mag
            locs_x[i] += force_x[i] * scale
            #locs_x[i]  = max(-1.0, min(locs_x[i], +1.0))
            locs_y[i] += force_y[i] * scale
            #locs_y[i]  = max(-1.0, min(locs_y[i], +1.0))
        end
    end

    # Scale to unit square
    min_x, max_x = minimum(locs_x), maximum(locs_x)
    min_y, max_y = minimum(locs_y), maximum(locs_y)
    function scaler(z, a, b)
        2.0*((z - a)/(b - a)) - 1.0
    end
    map!(z -> scaler(z, min_x, max_x), locs_x)
    map!(z -> scaler(z, min_y, max_y), locs_y)

    return locs_x,locs_y
end
