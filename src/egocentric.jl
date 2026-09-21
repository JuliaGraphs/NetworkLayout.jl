using LinearAlgebra: norm, eigen, Symmetric

export Egocentric, egocentric

"""
    Egocentric(; kwargs...)(adj_matrix)
    egocentric(adj_matrix; kwargs...)

Compute an egocentric ("focus") graph layout using stress majorization, centered
on a single focal vertex. Takes an adjacency matrix representation of a network
and returns coordinates of the nodes, translated such that the focal vertex sits
at the origin.

The layout interpolates between two objectives, following Brandes and Pich,
"More Flexible Radial Layout", Journal of Graph Algorithms and Applications
15(1):157-173 (2011,
[doi 10.7155/jgaa.00221](https://doi.org/10.7155/jgaa.00221)):

- a plain stress objective, which tries to match *all* pairwise euclidean
  distances to the corresponding graph distances (this is [`Stress`](@ref)), and
- a *focus* objective, which only weights the pairs involving the focal vertex.

Minimizing the focus objective alone places every vertex at a radius equal to
its graph distance from the focal vertex, producing concentric rings of
constant geodesic distance. Optimization proceeds along a schedule `tseq` of
mixing parameters `t ∈ [0,1]`, where the weights used in iteration stage `t` are
`(1-t)*W + t*Z` with `W[i,j] = d[i,j]^-2` and `Z` equal to `W` on the row and
column of the focal vertex and zero elsewhere. Each stage is majorized to
convergence before moving on to the next. Starting at `t=0` and ending at `t=1`
therefore uses the unconstrained stress layout to pick sensible *angles*, then
gradually enforces the radii.

## Inputs:
- `adj_matrix`: Matrix of pairwise distances.

## Keyword Arguments
- `focus=1`: Index of the focal vertex. It is held at the origin, so all returned
  positions are relative to it. Pinning it has no effect.
- `dim=2`, `Ptype=Float64`: Determines dimension and output type `Point{dim,Ptype}`.
- `tseq=0.0:0.1:1.0`

  Schedule of mixing parameters, from pure stress (`t=0`) to pure focus (`t=1`).
  Must be non-empty with entries in `[0,1]`.

- `iterations=100`: maximum number of majorization steps *per* entry of `tseq`.
- `abstols=0.0`

  Absolute tolerance for convergence of stress. A stage terminates if the
  difference between two successive stresses is less than abstol.

- `reltols=10e-5`

  Relative tolerance for convergence of stress. A stage terminates if the
  improvement in stress relative to the current stress is less than reltol.

- `abstolx=10e-6`

  Absolute tolerance for convergence of layout. A stage terminates if the
  largest movement of any single node is less than abstolx.

- `maxdist=nothing`

  If given, radially compress every vertex further than `maxdist` from the focal
  vertex onto a narrow band outside `maxdist`, keeping its angle. This bounds the
  extent of the layout so that the interesting, near part of the network fills
  the frame. Vertices that are unconnected to the focal vertex are affected by
  this too (see `uncon_dist`). Without `maxdist` no compression is applied.

- `compress=log1p`

  How far past `maxdist` a vertex at excess radius `Δ = r - maxdist` is drawn: its
  new radius is `maxdist + compress(Δ)`. May be a function, a real number (all
  distant vertices land on a single ring at `maxdist + compress`), or `nothing`
  (equivalent to `0`, i.e. clamp onto the `maxdist` ring). Ignored when `maxdist`
  is `nothing`.

- `uncon_dist=(maxdist, Ncomps)->maxdist*Ncomps^(1/3)`

  Per default, unconnected vertices in the graph get a pairwise "ideal" distance
  which scales with the number of connected components and the maximum distance
  within the components.

- `initialpos=Point{dim,Ptype}[]`

  Provide `Vector` or `Dict` of initial positions. By default all positions are
  initialized using classical multidimensional scaling of the graph distances
  plus a small random jitter, which makes the result largely deterministic.
  Those positions will be overwritten using the key-val-pairs provided by this
  argument.

- `pin=[]`: Pin node positions (won't be updated). Can be given as `Vector` or `Dict`
   of node index -> value pairings. Values can be either
    - `(12, 4.0)` : overwrite initial position and pin
    - `true/false` : pin this position
    - `(true, false, false)` : only pin certain coordinates

   Pinned positions are given relative to the focal vertex (which sits at the
   origin), and are exempt from the `maxdist` compression.

- `seed=1`: Seed for the random jitter on the initial positions.
- `rng=DEFAULT_RNG[](seed)`

  Create rng based on seed. Defaults to `MersenneTwister`, can be specified
  by overwriting `DEFAULT_RNG[]`
"""
@addcall struct Egocentric{Dim,Ptype,FT<:AbstractFloat,TS,MD,CF,UF,RNG} <:
                IterativeLayout{Dim,Ptype}
    focus::Int
    tseq::TS
    iterations::Int
    abstols::FT
    reltols::FT
    abstolx::FT
    maxdist::MD
    compress::CF
    uncon_dist::UF
    initialpos::Dict{Int,Point{Dim,Ptype}}
    pin::Dict{Int,SVector{Dim,Bool}}
    rng::RNG
end

function Egocentric(; focus=1,
                    dim=2,
                    Ptype=Float64,
                    tseq=0.0:0.1:1.0,
                    iterations=100,
                    abstols=0.0,
                    reltols=10e-5,
                    abstolx=10e-6,
                    maxdist=nothing,
                    compress=log1p,
                    uncon_dist=(maxd, N) -> maxd * N^(1 / 3),
                    initialpos=[], pin=[],
                    seed=1, rng=DEFAULT_RNG[](seed))
    if !isempty(initialpos)
        dim, Ptype = infer_pointtype(initialpos)
        Ptype = promote_type(Float32, Ptype) # make sure to get at least f32 if given as int
    end

    focus > 0 || throw(ArgumentError("focus needs to be a valid vertex index, got $focus"))
    isempty(tseq) && throw(ArgumentError("tseq needs to be non-empty!"))
    all(t -> 0 ≤ t ≤ 1, tseq) || throw(ArgumentError("All entries of tseq need to be in [0,1]!"))
    iterations > 0 || throw(ArgumentError("Iterations need to be > 0"))

    _initialpos, _pin = _sanitize_initialpos_pin(dim, Ptype, initialpos, pin)

    _compress = _compressfun(compress)
    _maxdist = maxdist === nothing ? nothing : float(maxdist)

    FT = promote_type(Float64, typeof(abstols), typeof(reltols), typeof(abstolx))
    TS, MD, CF, UF, RNG = typeof(tseq), typeof(_maxdist), typeof(_compress), typeof(uncon_dist),
                          typeof(rng)
    return Egocentric{dim,Ptype,FT,TS,MD,CF,UF,RNG}(focus, tseq, iterations, FT(abstols),
                                                    FT(reltols), FT(abstolx),
                                                    _maxdist, _compress, uncon_dist,
                                                    _initialpos, _pin, rng)
end

# `Returns` is only available from 1.7 on, but compat allows 1.6
@static if !isdefined(Base, :Returns)
    struct Returns{V} <: Function
        value::V
    end
    (obj::Returns)(args...; kw...) = obj.value
end

"""
Normalize the `compress` keyword into a function `Δ -> extra radius`.
"""
_compressfun(f) = f
_compressfun(::Nothing) = Returns(0.0)
_compressfun(b::Bool) = b ? Returns(1.0) : Returns(0.0)
_compressfun(c::Real) = Returns(float(c))

"""
Iteration state of the [`Egocentric`](@ref) layout. `positions` are the positions
of the majorization itself, already focus-centered but not yet radially
compressed; the compression is applied on the way out in `_egoemit`. `pintarget`
holds the fixed positions of the pinned vertices.
"""
mutable struct EgocentricState{PT,FT}
    positions::Vector{PT}
    D::Matrix{FT}
    W::Matrix{FT}
    Z::Matrix{FT}
    wsum::Vector{FT}
    zsum::Vector{FT}
    pin::Union{Nothing,Vector{<:SVector}}
    pintarget::Vector{PT}
    stage::Int      # index into tseq
    iter::Int       # majorization steps taken within the current stage
    laststress::FT
    finished::Bool
end

function Base.iterate(iter::LayoutIterator{<:Egocentric{Dim,Ptype,FT}}) where {Dim,Ptype,FT}
    algo, δ = iter.algorithm, iter.adj_matrix
    N = assertsquare(δ)
    algo.focus ≤ N ||
        throw(ArgumentError("focus=$(algo.focus) is out of bounds for a graph with $N vertices!"))

    make_symmetric!(δ)
    D = pairwise_distance(δ, FT)
    _replace_unconnected!(D, algo.uncon_dist)

    W = zeros(FT, N, N)
    for j in 1:N, i in 1:N
        i == j && continue
        W[i, j] = D[i, j]^-2
    end

    # focus weights: W restricted to the row and column of the focal vertex
    Z = zeros(FT, N, N)
    Z[algo.focus, :] .= @view W[algo.focus, :]
    Z[:, algo.focus] .= @view W[:, algo.focus]

    positions = _egoinitialpos(algo, D, N)

    if isempty(algo.pin)
        pin = nothing
    else
        pin = [get(algo.pin, i, SVector{Dim,Bool}(false for _ in 1:Dim)) for i in 1:N]
    end

    state = EgocentricState(positions, D, W, Z, vec(sum(W; dims=2)), vec(sum(Z; dims=2)),
                            pin, copy(positions), 1, 0, zero(FT), false)
    state.laststress = _mixedstress(state, first(algo.tseq))

    return _egoemit(algo, state), state
end

function Base.iterate(iter::LayoutIterator{<:Egocentric}, state)
    algo = iter.algorithm

    state.finished && return nothing

    if state.stage > length(algo.tseq)
        # emit the final layout a second time: `layout` returns the second to
        # last item of the iterator
        state.finished = true
        return _egoemit(algo, state), state
    end

    t = algo.tseq[state.stage]
    oldpos = copy(state.positions)
    _majorize!(state.positions, algo.focus, state, t)

    moved = maximum(i -> norm(state.positions[i] - oldpos[i]), eachindex(oldpos))
    state.iter += 1

    newstress = _mixedstress(state, t)
    converged = (state.laststress - newstress) ≤ algo.reltols * state.laststress ||
                abs(state.laststress - newstress) ≤ algo.abstols ||
                moved ≤ algo.abstolx
    state.laststress = newstress

    if converged || state.iter ≥ algo.iterations
        state.stage += 1
        state.iter = 0
        if state.stage ≤ length(algo.tseq)
            state.laststress = _mixedstress(state, algo.tseq[state.stage])
        end
    end

    return _egoemit(algo, state), state
end

"""
One SMACOF majorization sweep with the mixed weights `(1-t)*W + t*Z`.

`pos` is updated in place, i.e. the update of node `i` already sees the new
positions of nodes `1:i-1`. The simultaneous (Jacobi) variant of this update can
oscillate instead of converging.

The focal vertex is never moved. Stress is translation invariant, so holding one
vertex fixed only fixes the gauge; keeping it at the origin means that user
supplied positions (`initialpos`, `pin`) and the returned layout live in the same,
focus-centered frame.
"""
function _majorize!(pos::Vector{PT}, focus::Int, state::EgocentricState, t) where {PT}
    D, W, Z = state.D, state.W, state.Z
    pin, pintarget = state.pin, state.pintarget
    for i in eachindex(pos)
        i == focus && continue # holds the gauge, see above
        pinned = pin === nothing ? nothing : pin[i]
        pinned !== nothing && all(pinned) && continue

        acc = zero(PT)
        for j in eachindex(pos)
            i == j && continue
            w = (1 - t) * W[i, j] + t * Z[i, j]
            iszero(w) && continue
            offset = pos[i] - pos[j]
            nrm = norm(offset)
            # coincident nodes carry no direction information, only pull towards `pos[j]`
            inv_nrm = nrm > 1e-5 ? inv(nrm) : zero(nrm)
            acc += w * (pos[j] + D[i, j] * offset * inv_nrm)
        end

        denom = (1 - t) * state.wsum[i] + t * state.zsum[i]
        iszero(denom) && continue
        new = acc / denom
        if pinned !== nothing && any(pinned)
            new = PT(ntuple(k -> pinned[k] ? pintarget[i][k] : new[k], length(new)))
        end
        pos[i] = new
    end
    return pos
end

"""
Stress of the current layout under the mixed weights of stage `t`.
"""
function _mixedstress(state::EgocentricState, t)
    (1 - t) * stress(state.positions, state.D, state.W) +
    t * stress(state.positions, state.D, state.Z)
end

"""
Copy out the current layout, applying the radial compression beyond `maxdist`.
The focal vertex already sits at the origin, see [`_majorize!`](@ref).
"""
function _egoemit(algo::Egocentric{Dim,Ptype}, state::EgocentricState) where {Dim,Ptype}
    pos = copy(state.positions)

    if algo.maxdist !== nothing
        maxdist = algo.maxdist
        for i in eachindex(pos)
            state.pin !== nothing && any(state.pin[i]) && continue
            r = norm(pos[i])
            r > maxdist || continue
            pos[i] = pos[i] * Ptype((maxdist + algo.compress(r - maxdist)) / r)
        end
    end
    return pos
end

"""
Vertices which can't reach each other end up with a distance of `typemax` or
`Inf` (the latter if the Floyd-Warshall sum overflowed). Replace both with a
finite "ideal" distance based on the number of connected components.
"""
function _replace_unconnected!(D::Matrix{T}, uncon_dist) where {T}
    unreachable(x) = !isfinite(x) || x ≥ typemax(T)
    any(unreachable, D) || return D

    maxd = zero(T)
    for d in D
        unreachable(d) || (maxd = max(maxd, d))
    end
    iszero(maxd) && (maxd = one(T))

    dist = T(uncon_dist(maxd, _count_components(D, unreachable)))
    for i in eachindex(D)
        unreachable(D[i]) && (D[i] = dist)
    end
    return D
end

function _count_components(D, unreachable)
    N = size(D, 1)
    seen = falses(N)
    ncomps = 0
    for i in 1:N
        seen[i] && continue
        ncomps += 1
        for j in i:N
            unreachable(D[i, j]) || (seen[j] = true)
        end
    end
    return ncomps
end

"""
Initial positions from classical multidimensional scaling of the graph distances
plus a small jitter, overwritten by the user provided `initialpos`.
"""
function _egoinitialpos(algo::Egocentric{Dim,Ptype}, D, N) where {Dim,Ptype}
    coords = _classical_mds(D, Dim)
    rng = copy(algo.rng)
    startpos = Vector{Point{Dim,Ptype}}(undef, N)
    for i in 1:N
        startpos[i] = Point{Dim,Ptype}(ntuple(k -> coords[k, i] + Ptype(0.2 * (rand(rng) - 0.5)),
                                              Dim))
    end

    # the MDS solution is centered on the centroid, move it onto the focal vertex
    # so that `initialpos` and `pin` are interpreted in the focus-centered frame
    origin = startpos[algo.focus]
    startpos .= startpos .- Ref(origin)

    for (k, v) in algo.initialpos
        startpos[k] = v
    end
    startpos[algo.focus] = zero(Point{Dim,Ptype})
    return startpos
end

"""
Classical (Torgerson) multidimensional scaling of the distance matrix `D` into
`dim` dimensions. Returns a `dim × N` matrix of coordinates. Falls back to zeros
if the eigendecomposition fails, in which case the jitter alone seeds the
majorization.
"""
function _classical_mds(D::Matrix{T}, dim::Int) where {T}
    N = size(D, 1)
    coords = zeros(T, dim, N)
    N < 2 && return coords

    # double centering of the squared distances
    D2 = D .^ 2
    rowmean = vec(sum(D2; dims=2)) ./ N
    grandmean = sum(rowmean) / N
    B = Matrix{T}(undef, N, N)
    for j in 1:N, i in 1:N
        B[i, j] = -(D2[i, j] - rowmean[i] - rowmean[j] + grandmean) / 2
    end

    local F
    try
        F = eigen(Symmetric(B))
    catch
        return coords
    end

    # eigen returns ascending eigenvalues, the largest ones are at the end
    for k in 1:min(dim, N)
        λ = F.values[end - k + 1]
        λ > 0 || continue
        coords[k, :] .= @views F.vectors[:, end - k + 1] .* sqrt(λ)
    end
    return coords
end
