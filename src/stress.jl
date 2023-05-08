using LinearAlgebra: checksquare, norm, pinv, mul!
using SparseArrays: SparseMatrixCSC

export Stress, stress

"""
    Stress(; kwargs...)(adj_matrix)
    stress(adj_matrix; kwargs...)

Compute graph layout using stress majorization. Takes adjacency matrix
representation of a network and returns coordinates of the nodes.

The main equation to solve is (8) in Gansner, Koren and North (2005,
[doi 10.1007/978-3-540-31843-9_25](https://doi.org/10.1007/978-3-540-31843-9_25)).

## Inputs:
- `adj_matrix`: Matrix of pairwise distances.

## Keyword Arguments
- `dim=2`, `Ptype=Float64`: Determines dimension and output type `Point{dim,Ptype}`.
- `iterations=:auto`: maximum number of iterations (`:auto` means `400*N^2` where `N` are the number of vertices)
- `abstols=0`

   Absolute tolerance for convergence of stress. The iterations terminate if the
   difference between two successive stresses is less than abstol.

- `reltols=10e-6`

  Relative tolerance for convergence of stress. The iterations terminate if the
  difference between two successive stresses relative to the current stress is
  less than reltol.

- `abstolx=10e-6`

  Absolute tolerance for convergence of layout. The iterations terminate if the
  Frobenius norm of two successive layouts is less than abstolx.

- `weights=Array{Float64}(undef, 0, 0)`

  Matrix of weights. If empty (i.e. not specified), defaults to `weights[i,j] = δ[i,j]^-2` if
  `δ[i,j]` is nonzero, or `0` otherwise.

- `initialpos=Point{dim,Ptype}[]`

  Provide `Vector` or `Dict` of initial positions. All positions will be
  initialized using random coordinates from normal distribution. Random
  positions will be overwritten using the key-val-pairs provided by this
  argument.

- `pin=[]`: Pin node positions (won't be updated). Can be given as `Vector` or `Dict`
   of node index -> value pairings. Values can be either
    - `(12, 4.0)` : overwrite initial position and pin
    - `true/false` : pin this position
    - `(true, false, false)` : only pin certain coordinates

- `seed=1`: Seed for random initial positions.
"""
@addcall struct Stress{Dim,Ptype,IT<:Union{Symbol,Int},FT<:AbstractFloat,M<:AbstractMatrix} <:
                IterativeLayout{Dim,Ptype}
    iterations::IT
    abstols::FT
    reltols::FT
    abstolx::FT
    weights::M
    initialpos::Dict{Int,Point{Dim,Ptype}}
    pin::Dict{Int,SVector{Dim,Bool}}
    seed::UInt
end

function Stress(; dim=2,
                Ptype=Float64,
                iterations=:auto,
                abstols=0.0,
                reltols=10e-6,
                abstolx=10e-6,
                weights=Array{Float64}(undef, 0, 0),
                initialpos=[], pin=[],
                seed=1)
    if !isempty(initialpos)
        dim, Ptype = infer_pointtype(initialpos)
        Ptype = promote_type(Float32, Ptype) # make sure to get at least f32 if given as int
    end

    _initialpos, _pin = _sanitize_initialpos_pin(dim, Ptype, initialpos, pin)

    IT, FT, WT = typeof(iterations), typeof(abstols), typeof(weights)
    Stress{dim,Ptype,IT,FT,WT}(iterations, abstols, reltols, abstolx, weights, _initialpos, _pin, seed)
end

function Base.iterate(iter::LayoutIterator{<:Stress{Dim,Ptype,IT,FT}}) where {Dim,Ptype,IT,FT}
    algo, δ = iter.algorithm, iter.adj_matrix
    N = size(δ, 1)
    M = length(algo.initialpos)
    rng = MersenneTwister(algo.seed)
    startpos = randn(rng, Point{Dim,Ptype}, N)

    for (k, v) in algo.initialpos
        startpos[k] = v
    end

    if isempty(algo.pin)
        pin = nothing
    else
        isbitstype(Ptype) || error("Pin position only available for isbitstype (got $Ptype)!")
        pin = [get(algo.pin, i, SVector{Dim,Bool}(false for _ in 1:Dim)) for i in 1:N]
        pin = reinterpret(reshape, Bool, pin)
    end

    # calculate iteration if :auto
    maxiter = algo.iterations === :auto ? 400 * N^2 : algo.iterations
    @assert maxiter > 0 "Iterations need to be > 0"

    # if user provided weights not empty try those
    make_symmetric!(δ)
    distances = pairwise_distance(δ, FT)
    weights = isempty(algo.weights) ? distances .^ (-2) : algo.weights

    @assert length(startpos) == size(δ, 1) == size(δ, 2) == size(weights, 1) == size(weights, 2) "Wrong size of weights?"

    Lw = weightedlaplacian(weights)
    pinvLw = pinv(Lw)
    oldstress = stress(startpos, distances, weights)

    # the `state` of the iterator is (#iter, old stress, old pos, weights, distances pinvLw, pin, stopflag)
    return startpos, (1, oldstress, startpos, weights, distances, pinvLw, maxiter, pin, false)
end

function Base.iterate(iter::LayoutIterator{<:Stress{Dim,Ptype}}, state) where {Dim,Ptype}
    algo, δ = iter.algorithm, iter.adj_matrix
    i, oldstress, oldpos, weights, distances, pinvLw, maxiter, pin, stopflag = state

    if i >= maxiter || stopflag
        return nothing
    end

    # TODO the faster way is to drop the first row and col from the iteration
    t = LZ(oldpos, distances, weights)
    positions = similar(oldpos)
    mul!(positions, pinvLw, (t * oldpos))

    if !isnothing(pin)
        # on pin positions multiply newpos with zero and add oldpos
        _pos = reinterpret(reshape, Ptype, positions)
        _oldpos = reinterpret(reshape, Ptype, oldpos)
        _pos .= ((!).(pin) .* _pos) + (pin .* _oldpos)
    end

    @assert all(x -> all(map(isfinite, x)), positions)
    newstress = stress(positions, distances, weights)

    if abs(newstress - oldstress) < algo.reltols * newstress ||
       abs(newstress - oldstress) < algo.abstols ||
       norm(positions - oldpos) < algo.abstolx
        stopflag = true
    end

    return positions, (i + 1, newstress, positions, weights, distances, pinvLw, maxiter, pin, stopflag)
end

"""
Stress function to majorize

Input:
    positions: A particular layout (coordinates in rows)
    d: Matrix of pairwise distances
    weights: Weights for each pairwise distance

See (1) of Reference
"""
function stress(positions::AbstractArray{Point{T,N}}, d, weights) where {T,N}
    s = zero(T)
    n = length(positions)
    @assert n == size(d, 1) == size(d, 2) == size(weights, 1) == size(weights, 2)
    for j in 1:n, i in 1:(j - 1)
        s += weights[i, j] * (norm(positions[i] - positions[j]) - d[i, j])^2
    end
    @assert isfinite(s)
    return s
end

"""
Compute weighted Laplacian given ideal weights weights

Lʷ defined in (4) of the Reference
"""
function weightedlaplacian(weights::AbstractMatrix{T}) where {T}
    n = checksquare(weights)
    Lw = zeros(T, n, n)
    for i in 1:n
        D = zero(T)
        for j in 1:n
            i == j && continue
            Lw[i, j] = -weights[i, j]
            D += weights[i, j]
        end
        Lw[i, i] = D
    end
    return Lw
end

"""
Computes L^Z defined in (5) of the Reference

Input: Z: current layout (coordinates)
       d: Ideal distances
       weights: weights
"""
function LZ(Z::AbstractVector{Point{N,T}}, d, weights) where {N,T}
    n = length(Z)
    L = zeros(T, n, n)
    for i in 1:n
        D = zero(T)
        for j in 1:n
            i == j && continue
            nrmz = norm(Z[i] - Z[j])
            nrmz == 0 && continue
            δ = weights[i, j] * d[i, j]
            L[i, j] = -δ / nrmz
            D -= -δ / nrmz
        end
        @assert isfinite(D)
        L[i, i] = D
    end
    return L
end

"""
   pairwise_distance(δ)

Calculate the pairwise distances of a the graph from a adjacency matrix
using the Floyd-Warshall algorithm.

https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
"""
function pairwise_distance(δ, ::Type{T}=Float64) where {T}
    N = size(δ, 1)
    d = Matrix{T}(undef, N, N)
    @inbounds for j in 1:N, i in 1:N
        if i == j
            d[i, j] = zero(eltype(d))
        elseif iszero(δ[i, j])
            d[i, j] = typemax(eltype(d))
        else
            d[i, j] = δ[i, j]
        end
    end

    @inbounds for k in 1:N, i in 1:N, j in 1:N
        if d[i, k] + d[k, j] < d[i, j]
            d[i, j] = d[i, k] + d[k, j]
        end
    end
    return d
end
