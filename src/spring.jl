"""
    Use the spring/repulsion model of Fruchterman and Reingold (1991):
        Attractive force:  f_a(d) =  d^2 / k
        Repulsive force:  f_r(d) = -k^2 / d
    where d is distance between two vertices and the optimal distance
    between vertices k is defined as C * sqrt( area / num_vertices )
    where C is a parameter we can adjust

    Arguments:
    adj_matrix    Adjacency matrix of some type. Non-zero of the eltype
                  of the matrix is used to determine if a link exists,
                  but currently no sense of magnitude
    C             Constant to fiddle with density of resulting layout
    iterations    Number of iterations we apply the forces
    initialtemp   Initial "temperature", controls movement per iteration
"""
module Spring

using GeometryTypes
using LinearAlgebra: norm

struct Layout{M<:AbstractMatrix, P<:AbstractVector, T<:AbstractFloat}
  adj_matrix::M
  positions::P
  C::T
  iterations::Int
  initialtemp::T
end

function Layout(
        adj_matrix,
        PT::Type{Point{N,T}}=Point{2, Float64};
        startpositions=(2*rand(typ, size(adj_matrix,1)) .- 1),
        C=2.0, iterations=100, initialtemp=2.0
    ) where {N, T}
    Layout(adj_matrix, startpositions, T(C), Int(iterations), T(initialtemp))
end

layout(adj_matrix, dim::Int; kw_args...) = layout(adj_matrix, Point{dim,Float64}; kw_args...)

function layout(
        adj_matrix, typ::Type{Point{N,T}}=Point{2, Float64};
        startpositions = (2*rand(typ, size(adj_matrix,1)) .- 1),
        kw_args...
    ) where {N, T}
    layout!(adj_matrix,startpositions;kw_args...)
end

function layout!(
         adj_matrix,
         startpositions::AbstractVector{Point{N,T}};
         kw_args...
    ) where {N, T}
    size(adj_matrix, 1) != size(adj_matrix, 2) && error("Adj. matrix must be square.")
    # Layout object for the graph
    network = Layout(adj_matrix, Point{N,T}; startpositions=startpositions, kw_args...)
    next = iterate(network)
    while next != nothing
        (i, state) = next
        next = iterate(network, state)
    end
    return network.positions
end

function iterate(network::Layout)
   network.iterations == 1 && return nothing
   return network, 1
end

function iterate(network::Layout{M,P,T}, state) where {M, P, T}
    # The optimal distance bewteen vertices
    adj_matrix = network.adj_matrix
    N = size(adj_matrix,1)
    force = zeros(eltype(P), N)
    locs = network.positions
    C = network.C
    iterations = network.iterations
    initialtemp = network.initialtemp
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
    temp = initialtemp / state
    # Now apply them, but limit to temperature
    for i = 1:N
        force_mag  = norm(force[i])
        scale      = min(force_mag, temp)/force_mag
        locs[i]   += force[i] * scale
    end

    network.iterations == state && return nothing
    return network, (state+1)
end

end #end of module
