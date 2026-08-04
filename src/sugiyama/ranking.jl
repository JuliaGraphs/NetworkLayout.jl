# Phase 1: assign each vertex a rank (layer).
#
# Follows Gansner, Koutsofios, North & Vo (1993), "A Technique for Drawing
# Directed Graphs": build a feasible spanning tree that is *tight* (every
# tree edge has zero slack), then repeatedly find a tree edge with negative
# cut value and swap it for a non-tree edge until all cut values are
# non-negative, which is optimal (minimizes the total weighted edge length).
#
# `low`/`lim` are a postorder DFS numbering of the feasible tree that lets
# `enter_edge` decide in O(1) whether a candidate edge goes from the tail
# component to the head component after removing a tree edge, without
# re-walking the tree.

function toposort(g::SugiGraph)
    n = _nv(g)
    indeg = zeros(Int, n)
    for e in g.edges
        e === nothing && continue
        indeg[e.head] += 1
    end
    queue = [v for v in 1:n if indeg[v] == 0]
    order = Int[]
    sizehint!(order, n)
    while !isempty(queue)
        v = popfirst!(queue)
        push!(order, v)
        for eid in g.out[v]
            e = g.edges[eid]
            e === nothing && continue
            indeg[e.head] -= 1
            indeg[e.head] == 0 && push!(queue, e.head)
        end
    end
    length(order) == n || throw(ArgumentError("graph must be acyclic to rank it"))
    return order
end

function init_rank!(g::SugiGraph, minimum_length::Int)
    for v in toposort(g)
        maxr = nothing
        for eid in g.inn[v]
            cand = g.verts[g.edges[eid].tail].rank + minimum_length
            maxr = maxr === nothing ? cand : max(maxr, cand)
        end
        maxr !== nothing && (g.verts[v].rank = maxr)
    end
end

function move_vertices_up!(g::SugiGraph, minimum_length::Int)
    for v in 1:_nv(g)
        maxr = nothing
        for eid in g.inn[v]
            cand = g.verts[g.edges[eid].tail].rank + minimum_length
            maxr = maxr === nothing ? cand : max(maxr, cand)
        end
        g.verts[v].rank = maxr === nothing ? 0 : maxr
    end
end

function move_vertices_down!(g::SugiGraph, minimum_length::Int)
    _nv(g) == 0 && return
    max_rank = maximum(v.rank for v in g.verts)
    for v in 1:_nv(g)
        minr = nothing
        for eid in g.out[v]
            cand = g.verts[g.edges[eid].head].rank - minimum_length
            minr = minr === nothing ? cand : min(minr, cand)
        end
        g.verts[v].rank = minr === nothing ? max_rank : minr
    end
end

# ---- feasible tree construction -------------------------------------------

function tight_tree!(g::SugiGraph, v::Int, visited_edges::Set{Int}, minimum_length::Int)
    node_count = 1
    g.verts[v].is_tree_vertex = true
    for (eid, other) in collect(incident(g, v))
        eid in visited_edges && continue
        push!(visited_edges, eid)
        e = g.edges[eid]
        if e.is_tree_edge
            node_count += tight_tree!(g, other, visited_edges, minimum_length)
        elseif slack(g, eid, minimum_length) == 0 && !g.verts[other].is_tree_vertex
            e.is_tree_edge = true
            node_count += tight_tree!(g, other, visited_edges, minimum_length)
        end
    end
    return node_count
end

function is_incident_edge(g::SugiGraph, eid::Int)
    e = g.edges[eid]
    g.verts[e.tail].is_tree_vertex ⊻ g.verts[e.head].is_tree_vertex
end

function find_non_tight_edge(g::SugiGraph, minimum_length::Int)
    best, bestslack = nothing, nothing
    for eid in eachindex(g.edges)
        e = g.edges[eid]
        (e === nothing || e.is_tree_edge) && continue
        is_incident_edge(g, eid) || continue
        s = slack(g, eid, minimum_length)
        if bestslack === nothing || s < bestslack
            bestslack, best = s, eid
        end
    end
    best === nothing && throw(ArgumentError("no non-tight edge found while building feasible tree"))
    return best
end

function tighten_edge!(g::SugiGraph, delta::Int)
    for v in g.verts
        v.is_tree_vertex && (v.rank += delta)
    end
end

function feasible_tree!(g::SugiGraph, minimum_length::Int)
    root = 1
    while tight_tree!(g, root, Set{Int}(), minimum_length) < _nv(g)
        eid = find_non_tight_edge(g, minimum_length)
        e = g.edges[eid]
        delta = slack(g, eid, minimum_length)
        g.verts[e.head].is_tree_vertex && (delta = -delta)
        tighten_edge!(g, delta)
    end
    init_cutvalues!(g)
    init_low_lim!(g)
end

# ---- cut values -------------------------------------------------------

struct NeighborhoodInfo
    cut_value_sum::Int
    tree_edge_weight_sum::Int
    non_tree_edge_weight_sum::Int
    missing_v::Union{Nothing,Int}
end

function get_neighborhood_info(g::SugiGraph, v::Int, dir::Symbol)
    cut_value_sum = 0
    tree_edge_weight_sum = 0
    non_tree_edge_weight_sum = 0
    missing_v = nothing
    edge_ids = dir === :in ? g.inn[v] : g.out[v]
    for eid in edge_ids
        e = g.edges[eid]
        if !e.is_tree_edge
            non_tree_edge_weight_sum += e.weight
        elseif e.cut_value !== nothing
            cut_value_sum += e.cut_value
            tree_edge_weight_sum += e.weight
        elseif missing_v === nothing
            missing_v = dir === :in ? e.tail : e.head
        else
            return nothing
        end
    end
    return NeighborhoodInfo(cut_value_sum, tree_edge_weight_sum, non_tree_edge_weight_sum, missing_v)
end

function calculate_cut_value(edge_weight::Int, incoming::NeighborhoodInfo, outgoing::NeighborhoodInfo)
    edge_weight + incoming.non_tree_edge_weight_sum - incoming.cut_value_sum +
        incoming.tree_edge_weight_sum - outgoing.non_tree_edge_weight_sum +
        outgoing.cut_value_sum - outgoing.tree_edge_weight_sum
end

function calculate_cut_values!(g::SugiGraph, queue::Vector{Int})
    while !isempty(queue)
        v = popfirst!(queue)
        incoming = get_neighborhood_info(g, v, :in)
        outgoing = get_neighborhood_info(g, v, :out)
        (incoming === nothing || outgoing === nothing) && continue

        missing_v = if incoming.missing_v !== nothing && outgoing.missing_v === nothing
            incoming.missing_v
        elseif incoming.missing_v === nothing && outgoing.missing_v !== nothing
            outgoing.missing_v
        else
            continue
        end

        eid = find_edge(g, v, missing_v)
        if eid !== nothing
            incoming, outgoing = outgoing, incoming
        else
            eid = find_edge(g, missing_v, v)
        end
        g.edges[eid].cut_value = calculate_cut_value(g.edges[eid].weight, incoming, outgoing)
        push!(queue, missing_v)
    end
end

function leaves(g::SugiGraph)
    return [v for v in 1:_nv(g) if tree_degree(g, v) == 1]
end

init_cutvalues!(g::SugiGraph) = calculate_cut_values!(g, leaves(g))

function update_cutvalues!(g::SugiGraph, removed_edge::Int, swap_edge::Int)
    lca = remove_outdated_cut_values!(g, swap_edge, removed_edge)
    calculate_cut_values!(g, [g.edges[removed_edge].tail])
    return lca
end

function remove_outdated_cut_values!(g::SugiGraph, swap_edge::Int, removed_edge::Int)
    g.edges[removed_edge].cut_value = nothing
    w, x = g.edges[swap_edge].tail, g.edges[swap_edge].head
    g.verts[w].lim > g.verts[x].lim && ((w, x) = (x, w))

    lca = clear_cut_values_to_lca!(g, w, x)

    l = x
    while l != lca
        parent = g.verts[l].parent
        eid = find_edge_undirected(g, l, parent)
        g.edges[eid].cut_value = nothing
        l = parent
    end
    return lca
end

"""Clear cut values along the tree path from `w` towards the root, stopping
at (and returning) the least common ancestor of `w` and `x`."""
function clear_cut_values_to_lca!(g::SugiGraph, w::Int, x::Int)
    parent = g.verts[w].parent
    parent === nothing && return w
    l = w
    while true
        eid = find_edge_undirected(g, l, parent)
        g.edges[eid].cut_value = nothing
        l = parent
        if (g.verts[l].low <= g.verts[w].lim && g.verts[x].lim <= g.verts[l].lim) || g.verts[l].parent === nothing
            return l
        end
        parent = g.verts[l].parent
    end
end

# ---- low/lim (postorder DFS numbering of the feasible tree) -----------

function dfs_low_lim!(g::SugiGraph, v::Int, parent::Union{Nothing,Int}, max_lim::Base.RefValue{Int},
                       visited::BitVector)
    visited[v] = true
    g.verts[v].lim = max_lim[]
    g.verts[v].parent = parent
    for (eid, other) in collect(incident(g, v))
        if !visited[other] && g.edges[eid].is_tree_edge
            max_lim[] -= 1
            dfs_low_lim!(g, other, v, max_lim, visited)
        end
    end
    g.verts[v].low = max_lim[]
end

function init_low_lim!(g::SugiGraph)
    root = 1
    dfs_low_lim!(g, root, nothing, Ref(_nv(g)), falses(_nv(g)))
end

function update_low_lim!(g::SugiGraph, lca::Int)
    parent = g.verts[lca].parent
    visited = falses(_nv(g))
    parent !== nothing && (visited[parent] = true)
    dfs_low_lim!(g, lca, parent, Ref(g.verts[lca].lim), visited)
end

# ---- network simplex iteration -----------------------------------------

function leave_edge(g::SugiGraph)
    for eid in eachindex(g.edges)
        e = g.edges[eid]
        e === nothing && continue
        e.is_tree_edge && e.cut_value !== nothing && e.cut_value < 0 && return eid
    end
    return nothing
end

function is_head_to_tail(g::SugiGraph, eid::Int, u::Int, is_root_in_head::Bool)
    e = g.edges[eid]
    head_in_tail_component = g.verts[u].low <= g.verts[e.head].lim <= g.verts[u].lim
    tail_in_head_component = g.verts[u].low <= g.verts[e.tail].lim <= g.verts[u].lim
    return (is_root_in_head == head_in_tail_component) && (is_root_in_head != tail_in_head_component)
end

function enter_edge(g::SugiGraph, eid::Int, minimum_length::Int)
    e = g.edges[eid]
    u, v = e.tail, e.head
    is_root_in_head = g.verts[u].lim < g.verts[v].lim
    is_root_in_head || ((u, v) = (v, u))

    best, bestslack = nothing, nothing
    for cid in eachindex(g.edges)
        ce = g.edges[cid]
        (ce === nothing || ce.is_tree_edge) && continue
        is_head_to_tail(g, cid, u, is_root_in_head) || continue
        s = slack(g, cid, minimum_length)
        if bestslack === nothing || s < bestslack
            bestslack, best = s, cid
        end
    end
    best === nothing && throw(ArgumentError("no entering edge found (this should not happen)"))
    return best
end

function exchange!(g::SugiGraph, removed_edge::Int, swap_edge::Int, minimum_length::Int)
    g.edges[removed_edge].is_tree_edge = false
    g.edges[swap_edge].is_tree_edge = true
    lca = update_cutvalues!(g, removed_edge, swap_edge)
    update_low_lim!(g, lca)
    update_ranks!(g, minimum_length)
end

function update_neighbor_ranks!(g::SugiGraph, parent::Int, dir::Symbol, coeff::Int,
                                 queue::Vector{Int}, visited::BitVector, minimum_length::Int)
    edge_ids = dir === :out ? g.out[parent] : g.inn[parent]
    for eid in edge_ids
        e = g.edges[eid]
        e.is_tree_edge || continue
        other = dir === :out ? e.head : e.tail
        visited[other] && continue
        g.verts[other].rank = g.verts[parent].rank + minimum_length * coeff
        push!(queue, other)
        visited[other] = true
    end
end

function update_ranks!(g::SugiGraph, minimum_length::Int)
    node = 1
    visited = falses(_nv(g))
    visited[node] = true
    g.verts[node].rank = 0
    queue = [node]
    while !isempty(queue)
        parent = popfirst!(queue)
        update_neighbor_ranks!(g, parent, :out, 1, queue, visited, minimum_length)
        update_neighbor_ranks!(g, parent, :in, -1, queue, visited, minimum_length)
    end
end

"""Shift ranks so the smallest one is `1` (ranks are later used directly as
`layers` row indices, see `coordinates.jl`)."""
function normalize_ranks!(g::SugiGraph)
    _nv(g) == 0 && return
    shift = 1 - minimum(v.rank for v in g.verts)
    for v in g.verts
        v.rank += shift
    end
end

function minimize_edge_length!(g::SugiGraph, minimum_length::Int)
    feasible_tree!(g, minimum_length)
    while (e = leave_edge(g)) !== nothing
        swap_edge = enter_edge(g, e, minimum_length)
        exchange!(g, e, swap_edge, minimum_length)
    end
end

"""
    rank!(g, minimum_length, ranking_type)

Assign a `rank` to every vertex of `g`. `ranking_type` is one of:
- `:networksimplex`: optimal (minimizes total weighted edge length), via the
  network simplex method of Gansner et al. (1993).
- `:longestpath`: move every vertex as far up as possible, then pull sinks
  as far down as possible without violating edge lengths. Cheap but tends
  to produce wide, uneven layers.
- `:up`: longest path from the sources only.
- `:down`: longest path from the sinks only.
"""
function rank!(g::SugiGraph, minimum_length::Int, ranking_type::Symbol)
    init_rank!(g, minimum_length)
    if ranking_type === :networksimplex
        minimize_edge_length!(g, minimum_length)
    elseif ranking_type === :longestpath
        move_vertices_up!(g, minimum_length)
        move_vertices_down!(g, minimum_length)
    elseif ranking_type === :up
        move_vertices_up!(g, minimum_length)
    elseif ranking_type === :down
        move_vertices_down!(g, minimum_length)
    else
        throw(ArgumentError("Unknown ranking_type $(repr(ranking_type)), must be one of " *
                             ":networksimplex, :longestpath, :up, :down"))
    end
    normalize_ranks!(g)
end
