export Sugiyama, sugiyama

include("sugiyama/graph.jl")
include("sugiyama/cycles.jl")
include("sugiyama/ranking.jl")
include("sugiyama/ordering.jl")
include("sugiyama/coordinates.jl")

"""
    Sugiyama(; kwargs...)(adj_matrix)
    sugiyama(adj_matrix; kwargs...)

Layered ("hierarchical") layout for directed graphs: vertices are grouped
into ranks, edge crossings between ranks are heuristically minimized, and
coordinates are assigned within each rank. Ranking and crossing
minimization follow Gansner, Koutsofios, North & Vo (1993,
[doi 10.1109/32.221135](https://doi.org/10.1109/32.221135)); coordinate
assignment follows Brandes & Köpf (2002,
[doi 10.1007/3-540-45848-4_3](https://doi.org/10.1007/3-540-45848-4_3)).

Takes the adjacency matrix of a directed graph and returns coordinates of
the nodes. Cycles are broken by implicitly reversing edges; disconnected
components are laid out independently and placed side by side.

## Keyword Arguments
- `Ptype=Float64`: Determines the output type `Point{2,Ptype}`.
- `nodesize=Float64[]`: Size of each node. Filled up with `ones` or
  truncated to match the number of nodes.
- `nodespacing=1.0`: Minimum gap between neighboring nodes, both within a
  rank and between ranks.
- `dummysize=0.0`: Width of the invisible "dummy" nodes used internally to
  route edges spanning more than one rank.
- `minimum_length=1`: Minimum number of ranks every edge must span.
- `ranking_type=:networksimplex`: `:networksimplex` minimizes total edge
  length (as in the papers above); `:longestpath`, `:up` and `:down` are
  cheaper longest-path schedules (from both ends, from sources only, and
  from sinks only, respectively).
- `crossing_minimization=:barycenter`: `:barycenter` or `:median` heuristic
  used to reorder each rank.
- `transpose=true`: Follow up with a greedy pairwise-swap pass that further
  reduces crossings, at the cost of runtime.
- `direction=:down`: Which way ranks flow: `:down`, `:up`, `:left` or
  `:right`.

This implementation is a Julia port of
[rust-sugiyama](https://github.com/paddison/rust-sugiyama).
"""
@addcall struct Sugiyama{Ptype,T} <: AbstractLayout{2,Ptype}
    nodesize::Vector{T}
    nodespacing::Float64
    dummysize::Float64
    minimum_length::Int
    ranking_type::Symbol
    crossing_minimization::Symbol
    transpose::Bool
    direction::Symbol
end

function Sugiyama(; Ptype=Float64,
                   nodesize=Float64[],
                   nodespacing=1.0,
                   dummysize=0.0,
                   minimum_length=1,
                   ranking_type=:networksimplex,
                   crossing_minimization=:barycenter,
                   transpose=true,
                   direction=:down)
    minimum_length >= 1 || throw(ArgumentError("minimum_length must be >= 1"))
    direction in (:down, :up, :left, :right) ||
        throw(ArgumentError("direction must be one of :down, :up, :left, :right"))
    ranking_type in (:networksimplex, :longestpath, :up, :down) ||
        throw(ArgumentError("ranking_type must be one of :networksimplex, :longestpath, :up, :down"))
    crossing_minimization in (:barycenter, :median) ||
        throw(ArgumentError("crossing_minimization must be :barycenter or :median"))
    Sugiyama{Ptype,eltype(nodesize)}(nodesize, Float64(nodespacing), Float64(dummysize),
                                      Int(minimum_length), ranking_type, crossing_minimization,
                                      transpose, direction)
end

function layout(algo::Sugiyama{Ptype,T}, adj_matrix::AbstractMatrix) where {Ptype,T}
    n = assertsquare(adj_matrix)
    positions = Vector{Point{2,Ptype}}(undef, n)
    n == 0 && return positions

    nodesize = ones(Float64, n)
    for i in 1:min(n, length(algo.nodesize))
        nodesize[i] = Float64(algo.nodesize[i])
    end

    edgelist = Tuple{Int,Int}[]
    for j in 1:n, i in 1:n
        i != j && !iszero(adj_matrix[i, j]) && push!(edgelist, (i, j))
    end

    xoffset = 0.0
    for comp in _weakly_connected_components(n, edgelist)
        m = length(comp)
        localid = Dict(v => k for (k, v) in enumerate(comp))

        sub = SugiGraph()
        for v in comp
            _add_vertex!(sub; width=nodesize[v] + algo.nodespacing, height=nodesize[v] + algo.nodespacing)
        end
        for (i, j) in edgelist
            (haskey(localid, i) && haskey(localid, j)) || continue
            _add_edge!(sub, localid[i], localid[j])
        end

        xs, ys = layout_component!(sub, algo)

        minx, maxx = extrema(view(xs, 1:m))
        for (k, v) in enumerate(comp)
            positions[v] = to_point(Ptype, xs[k] - minx + xoffset, ys[k], algo.direction)
        end
        xoffset += (maxx - minx) + algo.nodespacing
    end
    return positions
end

function to_point(::Type{Ptype}, x::Float64, y::Float64, direction::Symbol) where {Ptype}
    if direction === :down
        Point{2,Ptype}(x, -y)
    elseif direction === :up
        Point{2,Ptype}(x, y)
    elseif direction === :right
        Point{2,Ptype}(y, x)
    else # :left
        Point{2,Ptype}(-y, x)
    end
end

"""Weakly connected components of the graph given by `n` vertices `1:n` and
directed edge list `edgelist`, as sorted vectors of (global) vertex ids."""
function _weakly_connected_components(n::Int, edgelist::Vector{Tuple{Int,Int}})
    adj = [Int[] for _ in 1:n]
    for (i, j) in edgelist
        push!(adj[i], j)
        push!(adj[j], i)
    end
    visited = falses(n)
    comps = Vector{Int}[]
    for start in 1:n
        visited[start] && continue
        comp = Int[]
        stack = [start]
        visited[start] = true
        while !isempty(stack)
            v = pop!(stack)
            push!(comp, v)
            for w in adj[v]
                if !visited[w]
                    visited[w] = true
                    push!(stack, w)
                end
            end
        end
        push!(comps, sort!(comp))
    end
    return comps
end

"""Run all four layout phases on a single weakly-connected component `g`
(vertices `1:m` are the "real" input vertices; the pipeline may append
dummy vertices after that). Returns `(xs, ys)` covering every vertex of the
(possibly grown) graph."""
function layout_component!(g::SugiGraph, algo::Sugiyama)
    remove_cycles!(g)
    rank!(g, algo.minimum_length, algo.ranking_type)
    insert_dummy_vertices!(g, algo.minimum_length, algo.dummysize)
    layers = ordering(g, algo.crossing_minimization, algo.transpose)

    layouts = create_layouts(g, layers)
    align_to_smallest_width_layout!(layouts)
    xcoords = calculate_relative_coords(layouts)

    n = _nv(g)
    xs = [xcoords[v] for v in 1:n]
    xs .-= minimum(xs)

    rank_to_max_height = Dict{Int,Float64}()
    for v in 1:n
        r = g.verts[v].rank
        rank_to_max_height[r] = max(get(rank_to_max_height, r, 0.0), g.verts[v].height)
    end
    ranks_sorted = sort!(collect(keys(rank_to_max_height)))
    rank_to_y = Dict{Int,Float64}()
    top = -rank_to_max_height[ranks_sorted[1]] * 0.5
    for r in ranks_sorted
        mh = rank_to_max_height[r]
        rank_to_y[r] = top + mh * 0.5
        top += mh
    end
    ys = [rank_to_y[g.verts[v].rank] for v in 1:n]

    return xs, ys
end
