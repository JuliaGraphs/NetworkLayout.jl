export Spring, spring

"""
    Spring(; kwargs...)(adj_matrix)
    spring(adj_matrix; kwargs...)

Use the spring/repulsion model of Fruchterman and Reingold (1991) with

- Attractive force:  `f_a(d) =  d^2 / k`
- Repulsive force:  `f_r(d) = -k^2 / d`

where `d` is distance between two vertices and the optimal distance between
vertices `k` is defined as `C * sqrt( area / num_vertices )` where `C` is a parameter
we can adjust

Takes adjacency matrix representation of a network and returns coordinates of
the nodes.

## Keyword Arguments
- `dim=2`, `Ptype=Float64`: Determines dimension and output type `Point{dim,Ptype}`.
- `C=2.0`: Constant to fiddle with density of resulting layout
- `iterations=100`: maximum number of iterations
- `initialtemp=2.0`: Initial "temperature", controls movement per iteration
- `initialpos=Point{dim,Ptype}[]`

  Provide list of initial positions. If length does not match Network size the initial
  positions will be truncated or filled up with random values between [-1,1] in every coordinate.

- `seed=1`: Seed for random initial positions.
"""
@addcall struct Spring{Dim,Ptype} <: IterativeLayout{Dim,Ptype}
    C::Float64
    iterations::Int
    initialtemp::Float64
    initialpos::Vector{Point{Dim,Ptype}}
    seed::UInt
end

function Spring(; dim=2, Ptype=Float64, C=2.0, iterations=100, initialtemp=2.0, initialpos=Point{dim,Ptype}[],
                seed=1)
    if !isempty(initialpos)
        initialpos = Point.(initialpos)
        Ptype = eltype(eltype(initialpos))
        # TODO fix initial pos if list has points of multiple types
        Ptype == Any && error("Please provide list of Point{N,T} with same T")
        dim = length(eltype(initialpos))
    end
    return Spring{dim,Ptype}(C, iterations, initialtemp, initialpos, seed)
end

function Base.iterate(iter::LayoutIterator{<:Spring{Dim,Ptype}}) where {Dim,Ptype}
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    N = size(adj_matrix, 1)
    M = length(algo.initialpos)
    rng = MersenneTwister(algo.seed)
    startpos = Vector{Point{Dim,Ptype}}(undef, N)
    # take the first
    for i in 1:min(N, M)
        startpos[i] = algo.initialpos[i]
    end
    # fill the rest with random points
    for i in (M + 1):N
        startpos[i] = 2 .* rand(rng, Point{Dim,Ptype}) .- 1
    end
    # iteratorstate: #iter nr, old pos
    return (startpos, (1, startpos))
end

function Base.iterate(iter::LayoutIterator{<:Spring}, state)
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    iteration, old_pos = state
    iteration >= algo.iterations && return nothing

    # The optimal distance bewteen vertices
    N = size(adj_matrix, 1)
    force = similar(old_pos)
    Ftype = eltype(force)
    K = algo.C * sqrt(4.0 / N)

    locs = copy(old_pos)
    # Calculate forces
    for i in 1:N
        force_vec = zero(Ftype)
        for j in 1:N
            i == j && continue
            d = norm(locs[j] .- locs[i])
            if adj_matrix[i, j] != zero(eltype(adj_matrix)) || adj_matrix[j, i] != zero(eltype(adj_matrix))
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
            force_vec += Ftype(F_d .* (locs[j] .- locs[i]))
        end
        force[i] = force_vec
    end
    # Cool down
    temp = algo.initialtemp / iteration
    # Now apply them, but limit to temperature
    for i in 1:N
        force_mag = norm(force[i])
        scale = min(force_mag, temp) ./ force_mag
        locs[i] += force[i] .* scale
    end

    return locs, (iteration + 1, locs)
end
