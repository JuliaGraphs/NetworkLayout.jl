module Spring
using GeometryTypes
export Layout, layout

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

immutable Layout{A<:AbstractMatrix, P<:AbstractVector, T}
  adj_matrix::A
  positions::P
  C::T
  MAXITER::Int
  INITTEMP::T
end

function layout{T}(
        adj_matrix::T, dim::Int,
        startpositions = (2*rand(Point{dim, Float64}, size(adj_matrix,1)) .- 1);
        C=2.0, MAXITER=100, INITTEMP=2.0
    )
    size(adj_matrix, 1) != size(adj_matrix, 2) && error("Adj. matrix must be square.")
    # Layout object for the graph
    network = Layout(adj_matrix,startpositions,C,MAXITER,INITTEMP)
    state = start(network)
    while !done(network,state)
        network,state = next(network,state)
    end
    return network.positions
end

Base.start(::Layout) = 1

function Base.next(network::Layout, state)
    iter = layout_iterator!(network, state)
    return network, iter
end

Base.done(network::Layout,state) = (state == network.MAXITER)

function layout_iterator!{A,P,T}(network::Layout{A,P,T}, iter)
    # The optimal distance bewteen vertices
    adj_matrix = network.adj_matrix
    N = size(adj_matrix,1)
    force = zeros(eltype(P),N)
    locs = network.positions
    C = network.C
    MAXITER = network.MAXITER
    INITTEMP = network.INITTEMP
    N = size(adj_matrix,1)
    Ftype = eltype(force)
    K = C * sqrt(4.0 / N)

    # Calculate forces
    for i = 1:N
        force_vec = Ftype(0)
        for j = 1:N
            i == j && continue
            d   = norm(locs[j]-locs[i])
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
            force_vec += Ftype(F_d*(locs[j]-locs[i]))
        end
        force[i] = force_vec
    end
    # Cool down
    TEMP = INITTEMP / iter
    # Now apply them, but limit to temperature
    for i = 1:N
        force_mag  = norm(force[i])
        scale      = min(force_mag, TEMP)/force_mag
        locs[i]   += force[i] * scale
    end
    return (iter+1)
end
end #end of module
