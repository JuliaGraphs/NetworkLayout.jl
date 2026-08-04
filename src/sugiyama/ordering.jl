# Phase 2: insert dummy vertices for edges spanning more than one rank, then
# reorder vertices within each rank to reduce edge crossings.
#
# Crossing minimization uses the classic median/barycenter bilayer-sweep
# heuristic (Sugiyama, Tagawa & Toda 1981 / Gansner et al. 1993): repeatedly
# sweep down then up through the ranks, reordering each rank by the
# median/average position of its already-fixed neighbors in the previous
# sweep direction, optionally followed by a greedy pairwise `transpose` pass
# that swaps adjacent vertices whenever doing so reduces the local crossing
# count. We keep the best ordering found and stop once a few sweeps in a row
# fail to improve on it (finding the true minimum is NP-hard).

function insert_dummy_vertices!(g::SugiGraph, minimum_length::Int, dummy_width::Float64)
    for eid in [eid for eid in eachindex(g.edges) if g.edges[eid] !== nothing]
        slack(g, eid, minimum_length) <= 0 && continue
        e = g.edges[eid]
        tail, head = e.tail, e.head
        _rem_edge!(g, eid)
        t = tail
        for r in (g.verts[tail].rank + 1):(g.verts[head].rank - 1)
            d = _add_vertex!(g; width=dummy_width, height=0.0, is_dummy=true)
            g.verts[d].rank = r
            _add_edge!(g, t, d)
            t = d
        end
        _add_edge!(g, t, head)
    end
end

# ---- vertex ordering within ranks --------------------------------------

mutable struct SugiOrder
    layers::Vector{Vector{Int}}
    positions::Vector{Int}     # position within layer, indexed by vertex id
end

function SugiOrder(layers::Vector{Vector{Int}}, n::Int)
    positions = zeros(Int, n)
    for layer in layers, (p, v) in enumerate(layer)
        positions[v] = p
    end
    return SugiOrder(layers, positions)
end

Base.copy(o::SugiOrder) = SugiOrder([copy(l) for l in o.layers], copy(o.positions))

function init_order(g::SugiGraph)
    n = _nv(g)
    n == 0 && return SugiOrder(Vector{Int}[], 0)
    max_rank = maximum(v.rank for v in g.verts)
    layers = [Int[] for _ in 1:max_rank]
    visited = falses(n)
    function dfs(v::Int)
        visited[v] && return
        visited[v] = true
        push!(layers[g.verts[v].rank], v)
        for n in collect(out_neighbors(g, v))
            dfs(n)
        end
    end
    for v in 1:n
        dfs(v)
    end
    return SugiOrder(layers, n)
end

function exchange!(order::SugiOrder, r::Int, i::Int, j::Int)
    layer = order.layers[r]
    order.positions[layer[i]], order.positions[layer[j]] = j, i
    layer[i], layer[j] = layer[j], layer[i]
    return order
end

function barycenter(g::SugiGraph, v::Int, move_down::Bool, positions::Vector{Int})
    neighbors = collect(move_down ? in_neighbors(g, v) : out_neighbors(g, v))
    isempty(neighbors) && return Float64(positions[v])
    return sum(positions[n] for n in neighbors) / length(neighbors)
end

function median(g::SugiGraph, v::Int, move_down::Bool, positions::Vector{Int})
    neighbors = move_down ? in_neighbors(g, v) : out_neighbors(g, v)
    adj = sort!([positions[n] for n in neighbors if abs(g.verts[v].rank - g.verts[n].rank) == 1])
    p = length(adj)
    p == 0 && return Inf
    m = p ÷ 2
    if isodd(p)
        return Float64(adj[m + 1])
    elseif p == 2
        return (adj[1] + adj[2]) / 2.0
    else
        left = adj[m] - adj[1]
        right = adj[p] - adj[m + 1]
        return (adj[m] * right + adj[m + 1] * left) / (left + right)
    end
end

function order_layer(g::SugiGraph, move_down::Bool, cur::SugiOrder, cm::F) where {F}
    nlayers = length(cur.layers)
    new_layers = [copy(l) for l in cur.layers]
    positions = copy(cur.positions)
    ranks = move_down ? (2:nlayers) : (nlayers - 1):-1:1

    for r in ranks
        layer = new_layers[r]
        scores = Dict(v => cm(g, v, move_down, positions) for v in layer)
        sort!(layer; by=v -> scores[v])
        for (p, v) in enumerate(layer)
            positions[v] = p
        end
    end
    return SugiOrder(new_layers, positions)
end

function bilayer_cross_count(g::SugiGraph, order::SugiOrder, rank::Int)
    north = order.layers[rank]
    south = order.layers[rank + 1]
    endpoints = Int[]
    for v in north
        for n in out_neighbors(g, v)
            abs(g.verts[v].rank - g.verts[n].rank) == 1 || continue
            push!(endpoints, order.positions[n])
        end
    end
    return count_crossings(endpoints, length(south))
end

"""Count inversions of `endpoints` (positions in `1:south_len`) using a
Fenwick tree; equivalent to (but simpler than) the accumulator-tree method
in the reference implementation."""
function count_crossings(endpoints::Vector{Int}, south_len::Int)
    south_len == 0 && return 0
    bit = zeros(Int, south_len)
    function bit_update!(i::Int)
        while i <= south_len
            bit[i] += 1
            i += i & (-i)
        end
    end
    function bit_query(i::Int)
        s = 0
        while i > 0
            s += bit[i]
            i -= i & (-i)
        end
        return s
    end
    cross = 0
    inserted = 0
    for pos in endpoints
        cross += inserted - bit_query(pos)
        bit_update!(pos)
        inserted += 1
    end
    return cross
end

function total_crossings(g::SugiGraph, order::SugiOrder)
    isempty(order.layers) && return 0
    return sum(bilayer_cross_count(g, order, r) for r in 1:(length(order.layers) - 1); init=0)
end

function cross_count_two_vertices(g::SugiGraph, order::SugiOrder, v::Int, w::Int)
    crossings = 0
    for dir in (:in, :out)
        v_adj = [order.positions[n] for n in (dir === :in ? in_neighbors(g, v) : out_neighbors(g, v))]
        w_adj = [order.positions[n] for n in (dir === :in ? in_neighbors(g, w) : out_neighbors(g, w))]
        for i in v_adj, j in w_adj
            i > j && (crossings += 1)
        end
    end
    return crossings
end

function transpose!(g::SugiGraph, order::SugiOrder, move_down::Bool)
    nlayers = length(order.layers)
    improved = true
    ranks = move_down ? (1:nlayers) : (nlayers:-1:1)
    while improved
        improved = false
        for r in ranks
            layer = order.layers[r]
            for i in 1:(length(layer) - 1)
                v, w = layer[i], layer[i + 1]
                vw = cross_count_two_vertices(g, order, v, w)
                wv = cross_count_two_vertices(g, order, w, v)
                if vw > wv
                    improved = true
                    exchange!(order, r, i, i + 1)
                end
            end
        end
    end
end

function reduce_crossings_bilayer_sweep(g::SugiGraph, order::SugiOrder, cm::F,
                                         do_transpose::Bool) where {F}
    length(order.layers) <= 1 && return order
    best_crossings = total_crossings(g, order)
    best = copy(order)
    cur = order
    last_best = 0
    i = 0
    while true
        cur = order_layer(g, iseven(i), cur, cm)
        do_transpose && transpose!(g, cur, iseven(i))
        crossings = total_crossings(g, cur)
        if crossings < best_crossings
            best_crossings = crossings
            best = copy(cur)
            last_best = 0
        else
            last_best += 1
        end
        last_best == 4 && return best
        i += 1
    end
end

"""
    ordering(g, crossing_minimization, transpose)

Return a `Vector{Vector{Int}}` giving, for every rank, the vertices in that
rank ordered left-to-right so as to (heuristically) minimize edge crossings.
`crossing_minimization` is `:barycenter` or `:median`.
"""
function ordering(g::SugiGraph, crossing_minimization::Symbol, do_transpose::Bool)
    order = init_order(g)
    cm = crossing_minimization === :barycenter ? barycenter :
         crossing_minimization === :median ? median :
         throw(ArgumentError("Unknown crossing_minimization $(repr(crossing_minimization)), " *
                              "must be :barycenter or :median"))
    order = reduce_crossings_bilayer_sweep(g, order, cm, do_transpose)
    return order.layers
end
