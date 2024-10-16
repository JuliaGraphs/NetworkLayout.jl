export Spring, spring

"""
    Spring(; kwargs...)(adj_matrix)
    spring(adj_matrix; kwargs...)

Use the spring/repulsion model of Fruchterman and Reingold (1991,
[doi 10.1002/spe.4380211102](https://doi.org/10.1002/spe.4380211102)) with

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

  Provide `Vector` or `Dict` of initial positions. All positions will be initialized
  using random coordinates between [-1,1]. Random positions will be overwritten using
  the key-val-pairs provided by this argument.

- `pin=[]`: Pin node positions (won't be updated). Can be given as `Vector` or `Dict`
   of node index -> value pairings. Values can be either
    - `(12, 4.0)` : overwrite initial position and pin
    - `true/false` : pin this position
    - `(true, false, false)` : only pin certain coordinates

- `seed=1`: Seed for random initial positions.
- `rng=DEFAULT_RNG[](seed)`

  Create rng based on seed. Defaults to `MersenneTwister`, can be specified
  by overwriting `DEFAULT_RNG[]`
"""
@addcall struct Spring{Dim,Ptype,RNG} <: IterativeLayout{Dim,Ptype}
    C::Float64
    iterations::Int
    initialtemp::Float64
    initialpos::Dict{Int,Point{Dim,Ptype}}
    pin::Dict{Int,SVector{Dim,Bool}}
    rng::RNG
end

function Spring(; dim=2, Ptype=Float64,
                C=2.0, iterations=100, initialtemp=2.0,
                initialpos=[], pin=[],
                seed=1, rng=DEFAULT_RNG[](seed))
    if !isempty(initialpos)
        dim, Ptype = infer_pointtype(initialpos)
        Ptype = promote_type(Float32, Ptype) # make sure to get at least f32 if given as int
    end
    _initialpos, _pin = _sanitize_initialpos_pin(dim, Ptype, initialpos, pin)

    return Spring{dim,Ptype,typeof(rng)}(C, iterations, initialtemp, _initialpos, _pin, rng)
end

function Base.iterate(iter::LayoutIterator{<:Spring{Dim,Ptype}}) where {Dim,Ptype}
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    N = size(adj_matrix, 1)
    rng = copy(algo.rng)
    startpos = [2 .* rand(rng, Point{Dim,Ptype}) .- 1 for _ in 1:N]

    for (k, v) in algo.initialpos
        startpos[k] = v
    end

    pin = [get(algo.pin, i, SVector{Dim,Bool}(false for _ in 1:Dim)) for i in 1:N]

    # iteratorstate: #iter nr, old pos, pin
    return (startpos, (1, startpos, pin, rng))
end

function Base.iterate(iter::LayoutIterator{<:Spring}, state)
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    iteration, old_pos, pin, rng = state
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
            if !iszero(d)
                force_vec += Ftype(F_d .* (locs[j] .- locs[i]))
            else
                # if two points are at the exact same location
                # use random force in any direction
                force_vec += randn(rng, Ftype)
            end

        end
        force[i] = force_vec
    end
    # Cool down
    temp = algo.initialtemp / iteration
    # Now apply them, but limit to temperature
    for i in 1:N
        force_mag = norm(force[i])
        iszero(force_mag) && continue
        scale = min(force_mag, temp) ./ force_mag

        mask = (!).(pin[i]) # where pin=true mask will multiply with 0
        locs[i] += force[i] .* scale .* mask
    end

    return locs, (iteration + 1, locs, pin, rng)
end
