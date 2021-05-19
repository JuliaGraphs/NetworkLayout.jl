using LinearAlgebra: checksquare, norm, pinv, mul!
using SparseArrays: SparseMatrixCSC

"""
Compute graph layout using stress majorization

Inputs:

    δ: Matrix of pairwise distances
    p: Dimension of embedding (default: 2)
    weights: Matrix of weights. If not specified, defaults to
           weights[i,j] = δ[i,j]^-2 if δ[i,j] is nonzero, or 0 otherwise
    X0: Initial guess for the layout. Coordinates are given in rows.
        If not specified, default to random matrix of Gaussians

Additional optional keyword arguments control the convergence of the algorithm
and the additional output as requested:

    iterations:   Maximum number of iterations. Default: 400size(X0, 1)^2
    abstols:      Absolute tolerance for convergence of stress.
                  The iterations terminate if the difference between two
                  successive stresses is less than abstol.
                  Default: √(eps(eltype(X0))
    reltols:      Relative tolerance for convergence of stress.
                  The iterations terminate if the difference between two
                  successive stresses relative to the current stress is less than
                  reltol. Default: √(eps(eltype(X0))
    abstolx:      Absolute tolerance for convergence of layout.
                  The iterations terminate if the Frobenius norm of two successive
                  layouts is less than abstolx. Default: √(eps(eltype(X0))

Output:

    The final layout positions.

Reference:

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
struct Stress{Dim,Ptype,IT<:Union{Symbol,Int},FT<:AbstractFloat,M<:AbstractMatrix} <:
       IterativeLayout{Dim,Ptype}
    iterations::IT
    abstols::FT
    reltols::FT
    abstolx::FT
    weights::M
    initialpos::Vector{Point{Dim,Ptype}}
end

function Stress(; dim=2, Ptype=Float64, iterations=:auto, abstols=(√(eps(Float64))),
                reltols=(√(eps(Float64))), abstolx=(√(eps(Float64))), weights=Array{Float64}(undef, 0, 0),
                initialpos=Point{dim,Ptype}[])
    if !isempty(initialpos)
        initialpos = Point.(initialpos)
        Ptype = eltype(eltype(initialpos))
        # TODO fix initial pos if list has points of multiple types
        Ptype == Any && error("Please provide list of Point{N,T} with same T")
        dim = length(eltype(initialpos))
    end
    IT, FT, WT = typeof(iterations), typeof(abstols), typeof(weights)
    Stress{dim,Ptype,IT,FT,WT}(iterations, abstols, reltols, abstolx, weights, initialpos)
end

function initialweights(D, T)::SparseMatrixCSC{T,Int64}
    map(D) do d
        x = T(d^(-2.0))
        return isfinite(x) ? x : zero(T)
    end
end

function Base.iterate(iter::LayoutIterator{<:Stress{Dim,Ptype,IT,FT}}) where {Dim,Ptype,IT,FT}
    algo, δ = iter.algorithm, iter.adj_matrix
    N = size(δ, 1)
    M = length(algo.initialpos)
    startpos = Vector{Point{Dim,Ptype}}(undef, N)
    # take the first
    for i in 1:min(N, M)
        startpos[i] = algo.initialpos[i]
    end
    # fill the rest with random points
    for i in (M + 1):N
        startpos[i] = 2 .* rand(Point{Dim,Ptype}) .- 1
    end

    # calculate iteration if :auto
    maxiter = algo.iterations === :auto ? 400 * size(δ, 1)^2 : algo.iterations
    @assert maxiter > 0 "Iterations need to be > 0"

    # if user provided weights not empty try those
    weights = isempty(algo.weights) ? initialweights(δ, FT) : algo.weights

    @assert length(startpos) == size(δ, 1) == size(δ, 2) == size(weights, 1) == size(weights, 2) "Wrong size of weights?"

    Lw = weightedlaplacian(weights)
    pinvLw = pinv(Lw)
    s = stress(startpos, δ, weights)

    # the `state` of the iterator is (#iter, old stress, old pos, weights, pinvLw, stopflag)
    return startpos, (1, s, startpos, weights, pinvLw, maxiter, false)
end

function Base.iterate(iter::LayoutIterator{<:Stress{Dim,Ptype}}, state) where {Dim,Ptype}
    algo, δ = iter.algorithm, iter.adj_matrix
    i, oldstress, oldpos, weights, pinvLw, maxiter, stopflag = state
    # newstress, oldstress, X0, i = state

    if i >= maxiter || stopflag
        return nothing
    end

    # TODO the faster way is to drop the first row and col from the iteration
    t = LZ(oldpos, δ, weights)
    positions = similar(oldpos) # allocate new array but keep type of oldpos
    mul!(positions, pinvLw, (t * oldpos))
    @assert all(x -> all(map(isfinite, x)), positions)
    newstress = stress(positions, δ, weights)

    if abs(newstress - oldstress) < algo.reltols * newstress ||
       abs(newstress - oldstress) < algo.abstols ||
       norm(positions - oldpos) < algo.abstolx
        stopflag = true
    end

    return positions, (i + 1, newstress, positions, weights, pinvLw, maxiter, stopflag)
end

"""
Stress function to majorize

Input:
    positions: A particular layout (coordinates in rows)
    d: Matrix of pairwise distances
    weights: Weights for each pairwise distance

See (1) of Reference
"""
function stress(positions::AbstractArray{Point{T,N}}, d=ones(T, length(positions), length(positions)),
                weights=initialweights(d, T)) where {T,N}
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
       d: Ideal distances (default: all 1)
       weights: weights (default: d.^-2)
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
