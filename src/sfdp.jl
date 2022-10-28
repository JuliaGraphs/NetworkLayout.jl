export SFDP, sfdp

"""
    SFDP(; kwargs...)(adj_matrix)
    sfdp(adj_matrix; kwargs...)

Using the Spring-Electric [model suggested by Yifan Hu](http://yifanhu.net/PUB/graph_draw_small.pdf).
Forces are calculated as:

        f_attr(i,j) = ‖xi - xj‖ ² / K ,    i<->j
        f_repln(i,j) = -CK² / ‖xi - xj‖ ,  i!=j

Takes adjacency matrix representation of a network and returns coordinates of
the nodes.

## Keyword Arguments
- `dim=2`, `Ptype=Float64`: Determines dimension and output type `Point{dim,Ptype}`.
- `tol=1.0`: Stop if position changes of last step `Δp <= tol*K` for all nodes
- `C=0.2`, `K=1.0`: Parameters to tweak forces.
- `iterations=100`: maximum number of iterations
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
"""
@addcall struct SFDP{Dim,Ptype,T<:AbstractFloat} <: IterativeLayout{Dim,Ptype}
    tol::T
    C::T
    K::T
    iterations::Int
    initialpos::Dict{Int,Point{Dim,Ptype}}
    pin::Dict{Int,SVector{Dim,Bool}}
    seed::UInt
end

# TODO: check SFDP default parameters
function SFDP(; dim=2, Ptype=Float64,
              tol=1.0, C=0.2, K=1.0,
              iterations=100,
              initialpos=[], pin=[],
              seed=1)
    if !isempty(initialpos)
        dim, Ptype = infer_pointtype(initialpos)
        Ptype = promote_type(Float32, Ptype) # make sure to get at least f32 if given as int
    end
    _initialpos, _pin = _sanitize_initialpos_pin(dim, Ptype, initialpos, pin)

    return SFDP{dim,Ptype,typeof(tol)}(tol, C, K, iterations, _initialpos, _pin, seed)
end

function Base.iterate(iter::LayoutIterator{SFDP{Dim,Ptype,T}}) where {Dim,Ptype,T}
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    N = size(adj_matrix, 1)
    rng = MersenneTwister(algo.seed)
    startpos = [2 .* rand(rng, Point{Dim,Ptype}) .- 1 for _ in 1:N]

    for (k, v) in algo.initialpos
        startpos[k] = v
    end

    pin = [get(algo.pin, i, SVector{Dim,Bool}(false for _ in 1:Dim)) for i in 1:N]

    # iteratorstate: (#iter, energy, step, progress, old pos, pin, stopflag)
    return startpos, (1, typemax(T), one(T), 0, startpos, pin, false)
end

function Base.iterate(iter::LayoutIterator{<:SFDP}, state)
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    iter, energy0, step, progress, locs0, pin, stopflag = state
    K, C, tol = algo.K, algo.C, algo.tol

    # stop if stopflag (tol reached) or nr of iterations reached
    if iter >= algo.iterations || stopflag
        return nothing
    end

    locs = copy(locs0)
    energy = zero(energy0)
    Ftype = eltype(locs)
    N = size(adj_matrix, 1)
    for i in 1:N
        force = zero(Ftype)
        for j in 1:N
            i == j && continue
            if adj_matrix[i, j] == 1
                # Attractive forces for adjacent nodes
                force += Ftype(f_attr(locs[i], locs[j], K) .*
                               ((locs[j] .- locs[i]) / norm(locs[j] .- locs[i])))
            else
                # Repulsive forces
                force += Ftype(f_repln(locs[i], locs[j], C, K) .*
                               ((locs[j] .- locs[i]) / norm(locs[j] .- locs[i])))
            end
        end
        if any(isnan, force)
            # if two points are at the exact same location
            # use random force in any direction
            rng = MersenneTwister(algo.seed + i)
            force = randn(rng, Ftype)
        end
        mask = (!).(pin[i]) # where pin=true mask will multiply with 0
        locs[i] = locs[i] .+ (step .* (force ./ norm(force))) .* mask
        energy = energy + norm(force)^2
    end
    step, progress = update_step(step, energy, energy0, progress)

    # if the tolerance is reached set stopflag to keep claculated point but stop next iteration
    if dist_tolerance(locs, locs0, K, tol)
        stopflag = true
    end

    return locs, (iter + 1, energy, step, progress, locs, pin, stopflag)
end

# Calculate Attractive force
f_attr(a, b, K) = (norm(a .- b) .^ 2) ./ K
# Calculate Repulsive force
f_repln(a, b, C, K) = -C .* (K^2) / norm(a .- b)

function update_step(step, energy::T, energy0, progress) where {T}
    # cooldown step
    t = T(0.9)
    if energy < energy0
        progress = progress + 1
        if progress >= 5
            progress = 0
            step = step / t
        end
    else
        progress = 0
        step = t * step
    end
    return step, progress
end

function dist_tolerance(locs, locs0, K, tol)
    # check whether the layout is optimal
    for i in 1:size(locs, 1)
        if norm(locs[i] .- locs0[i]) >= K * tol
            return false
        end
    end
    return true
end
