# Internal mutable directed graph used to implement the Sugiyama-style
# layered layout (phases 0-3 below). Not a general purpose graph type: it
# only supports what the algorithm needs (add/remove edges and vertices,
# look up edges by endpoints, walk incident edges). Edge ids stay valid
# after removal (removed edges are tombstoned to `nothing` instead of
# compacted), since e.g. network simplex keeps edge ids around across
# mutations.

mutable struct SugiVertex
    rank::Int
    pos::Int                  # position within its layer
    low::Int
    lim::Int
    parent::Union{Nothing,Int} # `nothing` at the root of the feasible tree
    is_tree_vertex::Bool
    is_dummy::Bool
    root::Int
    align::Int
    shift::Float64
    sink::Int
    block_width::Float64
    width::Float64
    height::Float64
end

function SugiVertex(; width=1.0, height=1.0, is_dummy=false)
    # root/align/sink are placeholders, overwritten with the vertex's own id
    # right after construction (see `_add_vertex!`).
    SugiVertex(0, 0, 0, 0, nothing, false, is_dummy, 0, 0, Inf, 0, 0.0, width, height)
end

mutable struct SugiEdge
    tail::Int
    head::Int
    weight::Int
    cut_value::Union{Nothing,Int}
    is_tree_edge::Bool
    has_type1_conflict::Bool
end

SugiEdge(tail, head; weight=1) = SugiEdge(tail, head, weight, nothing, false, false)

struct SugiGraph
    verts::Vector{SugiVertex}
    edges::Vector{Union{Nothing,SugiEdge}}
    out::Vector{Vector{Int}}   # out[v] = ids of edges with tail == v
    inn::Vector{Vector{Int}}   # inn[v] = ids of edges with head == v
end

SugiGraph() = SugiGraph(SugiVertex[], Union{Nothing,SugiEdge}[], Vector{Int}[], Vector{Int}[])

_nv(g::SugiGraph) = length(g.verts)

function _add_vertex!(g::SugiGraph; kwargs...)
    push!(g.verts, SugiVertex(; kwargs...))
    push!(g.out, Int[])
    push!(g.inn, Int[])
    id = length(g.verts)
    v = g.verts[id]
    v.root = id
    v.align = id
    v.sink = id
    return id
end

function _add_edge!(g::SugiGraph, tail::Int, head::Int; weight=1)
    push!(g.edges, SugiEdge(tail, head; weight))
    eid = length(g.edges)
    push!(g.out[tail], eid)
    push!(g.inn[head], eid)
    return eid
end

function _rem_edge!(g::SugiGraph, eid::Int)
    e = g.edges[eid]
    g.edges[eid] = nothing
    deleteat!(g.out[e.tail], findfirst(==(eid), g.out[e.tail]))
    deleteat!(g.inn[e.head], findfirst(==(eid), g.inn[e.head]))
    return e
end

"""Swap the direction of every edge in the graph (in place)."""
function reverse_graph!(g::SugiGraph)
    for e in g.edges
        e === nothing && continue
        e.tail, e.head = e.head, e.tail
    end
    for v in 1:_nv(g)
        g.out[v], g.inn[v] = g.inn[v], g.out[v]
    end
    return g
end

out_neighbors(g::SugiGraph, v::Int) = (g.edges[eid].head for eid in g.out[v])
in_neighbors(g::SugiGraph, v::Int) = (g.edges[eid].tail for eid in g.inn[v])

"""Edge ids and the vertex at the other end, considering both directions."""
function incident(g::SugiGraph, v::Int)
    Iterators.flatten((((eid, g.edges[eid].head) for eid in g.out[v]),
                        ((eid, g.edges[eid].tail) for eid in g.inn[v])))
end

function find_edge(g::SugiGraph, tail::Int, head::Int)
    for eid in g.out[tail]
        g.edges[eid].head == head && return eid
    end
    return nothing
end

function find_edge_undirected(g::SugiGraph, a::Int, b::Int)
    e = find_edge(g, a, b)
    e !== nothing && return e
    return find_edge(g, b, a)
end

function tree_degree(g::SugiGraph, v::Int)
    c = 0
    for eid in g.out[v]
        g.edges[eid].is_tree_edge && (c += 1)
    end
    for eid in g.inn[v]
        g.edges[eid].is_tree_edge && (c += 1)
    end
    return c
end

"""slack of edge `eid`: how much longer it is than the minimum rank length."""
function slack(g::SugiGraph, eid::Int, minimum_length::Int)
    e = g.edges[eid]
    return g.verts[e.head].rank - g.verts[e.tail].rank - minimum_length
end
