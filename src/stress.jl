using LinearAlgebra: checksquare, norm, pinv, mul!
using SparseArrays: SparseMatrixCSC

export Stress, stress

"""
    Stress(; kwargs...)(adj_matrix)
    stress(adj_matrix; kwargs...)

Compute graph layout using stress majorization. Takes adjacency matrix
representation of a network and returns coordinates of the nodes.

## Inputs:
- `adj_matrix`: Matrix of pairwise distances.

## Keyword Arguments
- `dim=2`, `Ptype=Float64`: Determines dimension and output type `Point{dim,Ptype}`.
- `iterations=:auto`: maximum number of iterations (`:auto` means `400*N^2` where `N` are the number of vertices)
- `abstols=(√(eps(Float64)))`

   Absolute tolerance for convergence of stress. The iterations terminate if the
   difference between two successive stresses is less than abstol.

- `reltols=(√(eps(Float64)))`

  Relative tolerance for convergence of stress. The iterations terminate if the
  difference between two successive stresses relative to the current stress is
  less than reltol.

- `abstolx=(√(eps(Float64)))`

  Absolute tolerance for convergence of layout. The iterations terminate if the
  Frobenius norm of two successive layouts is less than abstolx.

- `weights=Array{Float64}(undef, 0, 0)`

  Matrix of weights. If empty (i.e. not specified), defaults to `weights[i,j] = δ[i,j]^-2` if
  `δ[i,j]` is nonzero, or `0` otherwise.

- `initialpos=Point{dim,Ptype}[]`

  Provide list of initial positions. If length does not match Network size the initial
  positions will be truncated or filled up with random normal distributed values in every coordinate.

- `seed=1`: Seed for random initial positions.

## Reference:
The main equation to solve is (8) of:

    @incollection{
        author = {Emden R Gansner and Yehuda Koren and Stephen North},
        title = {Graph Drawing by Stress Majorization}
        year={2005},
        isbn={978-3-540-24528-5},
        booktitle={Graph Drawing},
        seriesvolume={3383},
        series={Lecture Notes in Computer Science},
        editor={Pach, J\'anos},
        doi={10.1007/978-3-540-31843-9_25},
        publisher={Springer Berlin Heidelberg},
        pages={239--250},
    }
"""
@addcall struct Stress{Dim,Ptype,IT<:Union{Symbol,Int},FT<:AbstractFloat,M<:AbstractMatrix} <:
                IterativeLayout{Dim,Ptype}
    iterations::IT
    abstols::FT
    reltols::FT
    abstolx::FT
    weights::M
    initialpos::Vector{Point{Dim,Ptype}}
    seed::UInt
end

function Stress(; dim=2,
                Ptype=Float64,
                iterations=:auto,
                abstols=0.0,
                reltols=10e-6,
                abstolx=10e-6,
                weights=Array{Float64}(undef, 0, 0),
                initialpos=Point{dim,Ptype}[],
                seed=1)
    if !isempty(initialpos)
        initialpos = Point.(initialpos)
        Ptype = eltype(eltype(initialpos))
        # TODO fix initial pos if list has points of multiple types
        Ptype == Any && error("Please provide list of Point{N,T} with same T")
        dim = length(eltype(initialpos))
    end
    IT, FT, WT = typeof(iterations), typeof(abstols), typeof(weights)
    Stress{dim,Ptype,IT,FT,WT}(iterations, abstols, reltols, abstolx, weights, initialpos, seed)
end

function Base.iterate(iter::LayoutIterator{<:Stress{Dim,Ptype,IT,FT}}) where {Dim,Ptype,IT,FT}
    algo, δ = iter.algorithm, iter.adj_matrix
    N = size(δ, 1)
    M = length(algo.initialpos)
    rng = MersenneTwister(algo.seed)
    startpos = Vector{Point{Dim,Ptype}}(undef, N)
    # take the first
    for i in 1:min(N, M)
        startpos[i] = algo.initialpos[i]
    end
    # fill the rest with random points
    for i in (M + 1):N
        startpos[i] = randn(rng, Point{Dim,Ptype})
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

    # the `state` of the iterator is (#iter, old stress, old pos, weights, distances pinvLw, stopflag)
    return startpos, (1, oldstress, startpos, weights, distances, pinvLw, maxiter, false)
end

function Base.iterate(iter::LayoutIterator{<:Stress{Dim,Ptype}}, state) where {Dim,Ptype}
    algo, δ = iter.algorithm, iter.adj_matrix
    i, oldstress, oldpos, weights, distances, pinvLw, maxiter, stopflag = state

    if i >= maxiter || stopflag
        return nothing
    end

    # TODO the faster way is to drop the first row and col from the iteration
    t = LZ(oldpos, distances, weights)
    positions = similar(oldpos) # allocate new array but keep type of oldpos
    mul!(positions, pinvLw, (t * oldpos))
    @assert all(x -> all(map(isfinite, x)), positions)
    newstress = stress(positions, distances, weights)

    if abs(newstress - oldstress) < algo.reltols * newstress ||
       abs(newstress - oldstress) < algo.abstols ||
       norm(positions - oldpos) < algo.abstolx
        stopflag = true
    end

    return positions, (i + 1, newstress, positions, weights, distances, pinvLw, maxiter, stopflag)
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
