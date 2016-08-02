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

    maxiter:   Maximum number of iterations. Default: 400size(X0, 1)^2
    abstols:   Absolute tolerance for convergence of stress.
               The iterations terminate if the difference between two
               successive stresses is less than abstol.
               Default: √(eps(eltype(X0))
    reltols:   Relative tolerance for convergence of stress.
               The iterations terminate if the difference between two
               successive stresses relative to the current stress is less than
               reltol. Default: √(eps(eltype(X0))
    abstolx:   Absolute tolerance for convergence of layout.
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

using GeometryTypes, Compat, FixedSizeArrays
import Base: start, next, done, *

function (*){T<:LinAlg.BlasFloat,S<:FixedArray}(A::StridedMatrix{T}, x::StridedVector{S})
    A_mul_B!(similar(x, S, size(A,1)), A, x)
end
immutable Layout{M1<:AbstractMatrix, M2<:AbstractMatrix, VP<:AbstractVector,FT<:AbstractFloat}
    δ::M1
    weights::M2
    positions::VP
    pinvLw::Matrix{FT}
    maxiter::Int
    abstols::FT
    reltols::FT
    abstolx::FT
end


function initialweights(D, T=eltype(D))
    map(D) do d
        x = T(d^(-2.0))
        isfinite(x) ? x : zero(T)
    end
end

function Layout{N, T}(
        δ, PT::Type{Point{N, T}}=Point{2, Float64};
        startpositions=rand(PT, size(δ,1)), weights=initialweights(δ, T),
        iterations=400*size(δ,1)^2, abstols=√(eps(T)),
        reltols=√(eps(T)), abstolx=√(eps(T))
    )
    @assert size(startpositions, 1)==size(δ, 1)==size(δ, 2)==size(weights, 1)==size(weights, 2)
    Lw = weightedlaplacian(weights)
    pinvLw = pinv(Lw)
    return Layout(δ, weights, startpositions, pinvLw, iterations, abstols, reltols, abstolx)
end


function layout{N, T}(
        δ, PT::Type{Point{N, T}}=Point{2, Float64};
        startpositions=rand(PT, size(δ,1)), kw_args...
    )
    layout!(δ, startpositions; kw_args...)
end

function layout!{N, T}(
        δ, startpositions::AbstractVector{Point{N, T}};
        iterations=400*size(δ,1)^2, kw_args...
    )
    iter = Layout(δ, Point{N,T}; startpositions=startpositions, kw_args...)
    state = start(iter)
    while !done(iter, state)
        iter, state = next(iter, state)
    end
    state[end] > iterations && warn("Maximum number of iterations reached without convergence")
    iter.positions
end

function start(network::Layout)
    s = stress(network.positions, network.δ, network.weights)
    return (s, s, network.positions, 0)
end

function next(iter::Layout, state)
    newstress, oldstress, X0, i = state
    δ, weights, pinvLw, positions = iter.δ, iter.weights, iter.pinvLw, iter.positions
    #TODO the faster way is to drop the first row and col from the iteration
    t = LZ(X0, δ, weights)
    positions = pinvLw * (t*X0)
    @assert all(x->all(map(isfinite, x)), positions)
    newstress, oldstress = stress(positions, δ, weights), newstress
    iter.positions[:] = positions
    return iter, (newstress, oldstress, positions, (i+1))
end

function done(iter::Layout, state)
    newstress, oldstress, X0, i = state
    maxiter, reltols = iter.maxiter, iter.reltols
    positions, abstols, abstolx = iter.positions, iter.abstols, iter.abstolx
    (i == 0) && (i<maxiter) && return false #special case 0th iteration

    return (
        i > maxiter ||
        abs(newstress - oldstress) < reltols * newstress ||
        abs(newstress - oldstress) < abstols ||
        vecnorm(positions - X0) < abstolx
    )
end


"""
Stress function to majorize

Input:
    positions: A particular layout (coordinates in rows)
    d: Matrix of pairwise distances
    weights: Weights for each pairwise distance

See (1) of Reference
"""
function stress{T, N}(
        positions::AbstractArray{Point{T, N}},
        d=ones(T, length(positions), length(positions)), weights=initialweights(d, T)
    )
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
function weightedlaplacian{T}(weights::AbstractMatrix{T})
    n = Compat.LinAlg.checksquare(weights)
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
function LZ{N,T}(Z::AbstractVector{Point{N,T}}, d, weights)
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
