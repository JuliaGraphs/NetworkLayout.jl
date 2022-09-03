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

  Provide list of initial positions. If length does not match Network size the initial
  positions will be truncated or filled up with random values between [-1,1] in every coordinate.

- `seed=1`: Seed for random initial positions.
- `pin=Dict{Int,Bool}()` : Anchors positions of nodes to initial position
"""
@addcall struct SFDP{Dim,Ptype,T<:AbstractFloat} <: IterativeLayout{Dim,Ptype}
    tol::T
    C::T
    K::T
    iterations::Int
    initialpos::Dict{Int, Point{Dim,Ptype}}
    seed::UInt
    pin::Dict{Int, Bool}
    function SFDP(tol::T, C::T, K::T, iterations::Int, initialpos::Dict, seed::UInt, pin::Dict) where T<:AbstractFloat
        for (ix, p) in pin
            if !haskey(initialpos, ix) && p
                @warn "No coordinate provided for pinned position $ix"
            end
        end
        dim = get_pt_dim(initialpos)
        ptype = get_pt_ptype(initialpos)
        new{dim, ptype, typeof(tol)}(tol, C, K, iterations, initialpos, seed, pin)
    end
end

# TODO: check SFDP default parameters
function SFDP(; dim=2, Ptype=Float64, tol=1.0, C=0.2, K=1.0, iterations=100, initialpos=Dict{Int, Point{dim,Ptype}}(),
              seed::UInt=UInt(1), pin = Dict{Int, Bool}())
    return SFDP(tol, C, K, iterations, initialpos, seed, pin)
end

function SFDP(tol::T, C::T, K::T, iterations::Int, initialpos::Vector, seed::UInt, pin) where T<:AbstractFloat
    initialpos = Dict(zip(1:length(initialpos), Point.(initialpos)))
    # TODO fix initial pos if list has points of multiple types
    return SFDP(tol, C, K, iterations, initialpos, seed, pin)
end

function SFDP(tol::T, C::T, K::T, iterations::Int, initialpos::Dict{Int, <:Point}, seed::UInt, pin::Vector{Bool}) where T<:AbstractFloat
    dim = get_pt_dim(initialpos)
    fixed = Dict(zip(1:length(pin), pin))
    return SFDP(tol, C, K, iterations, initialpos, seed, fixed)
end

function get_pt_ptype(ip::Dict{Int, <:Point})
    ptype = eltype(eltype(values(ip)))
    ptype == Any && error("Please provide list of Point{N,T} with same T")
    return ptype
end

function get_pt_dim(ip::Vector)
    dims = length.(ip)
    length(unique(dims)) != 1 && error("Please provide list of Point{N,T} with same N")
    return first(dims)
end

function get_pt_dim(ip::Dict{Int, <:Point}) 
    return typeof(ip).parameters[2].parameters[1] #based on type of point
end

function Base.iterate(iter::LayoutIterator{SFDP{Dim,Ptype,T}}) where {Dim,Ptype,T}
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    N = size(adj_matrix, 1)
    M = length(algo.initialpos)
    rng = MersenneTwister(algo.seed)
    startpos = Vector{Point{Dim,Ptype}}(undef, N)
    # take the first
    for (i, p) in algo.initialpos
        startpos[i] = p
    end

    # create bounds for random initial positions
    if isempty(algo.initialpos)
        startposbounds = (min = zero(Point{Dim, Ptype}), max = zero((Point{Dim, Ptype})) .+ one(Ptype))
    else
        startposbounds = (
            min = Point{Dim, Ptype}([minimum((p[d] for p in values(algo.initialpos))) for d in 1:Dim]),
            max = Point{Dim, Ptype}([maximum((p[d] for p in values(algo.initialpos))) for d in 1:Dim])
        )
    end

    # fill the rest with random points
    for i in setdiff(1:N, keys(algo.initialpos))
        startpos[i] = (startposbounds.max .- startposbounds.min) .* (2 .* rand(rng, Point{Dim,Ptype}) .- 1) .+ startposbounds.min
    end
    # iteratorstate: (#iter, energy, step, progress, old pos, stopflag)
    return startpos, (1, typemax(T), one(T), 0, startpos, false)
end

function Base.iterate(iter::LayoutIterator{<:SFDP}, state)
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    iter, energy0, step, progress, locs0, stopflag = state
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
        if !get(algo.pin, i, false)
            locs[i] = locs[i] .+ step .* (force ./ norm(force))
        end
        energy = energy + norm(force)^2
    end
    step, progress = update_step(step, energy, energy0, progress)

    # if the tolerance is reached set stopflag to keep claculated point but stop next iteration
    if dist_tolerance(locs, locs0, K, tol)
        stopflag = true
    end

    return locs, (iter + 1, energy, step, progress, locs, stopflag)
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
