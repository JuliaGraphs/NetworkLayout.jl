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
module Stress

using GeometryTypes
using LinearAlgebra: checksquare, norm, pinv
using SparseArrays: SparseMatrixCSC

struct Layout{M1<:AbstractMatrix, M2<:AbstractMatrix, VP<:AbstractVector, FT<:AbstractFloat}
    δ::M1
    weights::M2
    positions::VP
    pinvLw::Matrix{FT}
    iterations::Int
    abstols::FT
    reltols::FT
    abstolx::FT
end


function initialweights(D, T=Float64)::SparseMatrixCSC{T,Int64}
    map(D) do d
        x = T(d^(-2.0))
        isfinite(x) ? x : zero(T)
    end
end

function Layout(
        δ, PT::Type{Point{N, T}}=Point{2, Float64};
        startpositions=rand(PT, size(δ,1)), weights=initialweights(δ,T),
        iterations=400*size(δ,1)^2, abstols=√(eps(T)),
        reltols=√(eps(T)), abstolx=√(eps(T))
    ) where {N, T}
    @assert size(startpositions, 1)==size(δ, 1)==size(δ, 2)==size(weights, 1)==size(weights, 2)
    Lw = weightedlaplacian(weights)
    pinvLw = pinv(Lw)
    return Layout(δ, weights, startpositions, pinvLw, iterations, abstols, reltols, abstolx)
end

layout(δ, dim::Int; kw_args...) = layout(δ, Point{dim,Float64}; kw_args...)

function layout(
        δ, PT::Type{Point{N, T}}=Point{2, Float64};
        startpositions=rand(PT, size(δ,1)), kw_args...
    ) where {N, T}
    layout!(δ, startpositions; kw_args...)
end

function layout!(
        δ, startpositions::AbstractVector{Point{N, T}};
        iterations=400*size(δ,1)^2, kw_args...
    ) where {N, T}
    iter = Layout(δ, Point{N,T}; startpositions=startpositions, kw_args...)
    num_iterations = 0
    next = iterate(iter)
    while next != nothing
        i, state = next
        next = iterate(iter, state)
        num_iterations += 1
    end
    num_iterations > iterations && @warn("Maximum number of iterations reached without convergence")
    iter.positions
end

function iterate(network::Layout)
    network.iterations == 0 && return nothing
    s = stress(network.positions, network.δ, network.weights)
    return network, (s, s, network.positions, 0)
end

function iterate(network::Layout, state)
    newstress, oldstress, X0, i = state
    δ, weights, pinvLw, positions, X0 = network.δ, network.weights, network.pinvLw, network.positions, copy(network.positions)
    #TODO the faster way is to drop the first row and col from the iteration
    t = LZ(X0, δ, weights)
    positions = pinvLw * (t*X0)
    @assert all(x->all(map(isfinite, x)), positions)
    newstress, oldstress = stress(positions, δ, weights), newstress
    network.positions[:] = positions

    if i > network.iterations ||
            abs(newstress - oldstress) < network.reltols * newstress ||
            abs(newstress - oldstress) < network.abstols ||
            norm(positions - X0) < network.abstolx
        return nothing
    end

    return network, (newstress, oldstress, X0, (i+1))

end

"""
Stress function to majorize

Input:
    positions: A particular layout (coordinates in rows)
    d: Matrix of pairwise distances
    weights: Weights for each pairwise distance

See (1) of Reference
"""
function stress(
        positions::AbstractArray{Point{T, N}},
        d=ones(T, length(positions), length(positions)), weights=initialweights(d, T)
    ) where {T, N}
    s = zero(T); n = length(positions)
    @assert n==size(d, 1)==size(d, 2)==size(weights, 1)==size(weights, 2)
    for j=1:n, i=1:j-1
        s += weights[i, j] * (norm(positions[i] - positions[j]) - d[i,j])^2
    end
    @assert isfinite(s)
    s
end



"""
Compute weighted Laplacian given ideal weights weights

Lʷ defined in (4) of the Reference
"""
function weightedlaplacian(weights::AbstractMatrix{T}) where {T}
    n = checksquare(weights)
    Lw = zeros(T, n, n)
    for i=1:n
        D = zero(T)
        for j=1:n
            i==j && continue
            Lw[i, j] = -weights[i, j]
            D += weights[i, j]
        end
        Lw[i, i] = D
    end
    Lw
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
    for i=1:n
        D = zero(T)
        for j=1:n
            i==j && continue
            nrmz = norm(Z[i] - Z[j])
            nrmz==0 && continue
            δ = weights[i, j] * d[i, j]
            L[i, j] = -δ/nrmz
            D -= -δ/nrmz
        end
        @assert isfinite(D)
        L[i, i] = D
    end
    L
end

end # end of module
