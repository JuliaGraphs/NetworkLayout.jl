# Phase 3: assign x/y coordinates given the ranks and per-rank ordering,
# using the horizontal-alignment method of Brandes & Köpf (2002), "Fast and
# Simple Horizontal Coordinate Assignment".
#
# The idea: align each vertex with (at most) one "median" neighbor per rank
# to build vertical chains ("blocks"), then compact each block as far left
# as it can go without overlapping its neighbors. Doing this once produces
# skewed results, so the algorithm is run 4 times (top-down/bottom-up ×
# left-biased/right-biased, via reversing the graph's edges and/or each
# rank's order) and the 4 resulting x-coordinates per vertex are combined by
# averaging the two middle (median) values, which Brandes & Köpf show is
# both order- and separation-preserving.
#
function reset_alignment!(g::SugiGraph, layers::Vector{Vector{Int}})
    for (r, layer) in enumerate(layers)
        for (p, v) in enumerate(layer)
            vert = g.verts[v]
            vert.rank = r
            vert.pos = p
            vert.shift = Inf
            vert.align = v
            vert.root = v
            vert.sink = v
        end
    end
end

function is_incident_to_inner_segment(g::SugiGraph, id::Int)
    g.verts[id].is_dummy || return false
    return any(g.verts[t].is_dummy for t in in_neighbors(g, id))
end

function get_inner_segment_upper_neighbor(g::SugiGraph, id::Int)
    is_incident_to_inner_segment(g, id) || return nothing
    for t in in_neighbors(g, id)
        return t
    end
    return nothing
end

# Must run after `reset_alignment!` has populated `.pos` for `layers`
# (rust-port #26). The inner loop includes `l_1` itself, not just
# everything strictly before it, or the last vertex of each boundary
# segment never gets checked (rust-port #27).
function mark_type1_conflicts!(g::SugiGraph, layers::Vector{Vector{Int}})
    for r in 1:(length(layers) - 1)
        level = layers[r]
        next_level = layers[r + 1]
        left_dummy_index = 0
        l = 0
        for l_1 in eachindex(next_level)
            dummy_candidate = next_level[l_1]
            upn = get_inner_segment_upper_neighbor(g, dummy_candidate)
            if upn !== nothing
                right_dummy_index = g.verts[upn].pos
            elseif l_1 == length(next_level)
                right_dummy_index = length(level)
            else
                continue
            end
            while l < l_1
                vertex = next_level[l + 1]
                upper_neighbors = sort!(collect(in_neighbors(g, vertex)); by=u -> g.verts[u].pos)
                for un in upper_neighbors
                    vertex_index = g.verts[un].pos
                    if vertex_index < left_dummy_index || vertex_index > right_dummy_index
                        eid = find_edge(g, un, vertex)
                        g.edges[eid].has_type1_conflict = true
                    end
                end
                l += 1
            end
            left_dummy_index = right_dummy_index
        end
    end
end

function create_vertical_alignments!(g::SugiGraph, layers::Vector{Vector{Int}})
    for layer in layers
        r = 0 # sentinel: no neighbor aligned yet
        for v in layer
            edges = [(eid, g.edges[eid].tail) for eid in g.inn[v] if slack(g, eid, 1) == 0]
            isempty(edges) && continue
            sort!(edges; by=x -> g.verts[x[2]].pos)

            d = (length(edges) + 1) / 2 - 1
            for m in (floor(Int, d), ceil(Int, d))
                g.verts[v].align == v || continue
                eid, median_neighbor = edges[m + 1]
                if !g.edges[eid].has_type1_conflict && r < g.verts[median_neighbor].pos
                    g.verts[median_neighbor].align = v
                    g.verts[v].root = g.verts[median_neighbor].root
                    g.verts[v].align = g.verts[v].root
                    r = g.verts[median_neighbor].pos
                end
            end
        end
    end
end

function compute_block_max_vertex_widths!(g::SugiGraph)
    for root in 1:_nv(g)
        g.verts[root].root == root || continue
        maxw = g.verts[root].width
        cur = g.verts[root].align
        while cur != root
            maxw = max(maxw, g.verts[cur].width)
            cur = g.verts[cur].align
        end
        g.verts[root].block_width = maxw
        cur = g.verts[root].align
        while cur != root
            g.verts[cur].block_width = maxw
            cur = g.verts[cur].align
        end
    end
end

"""The vertex immediately to the left of `v` within its own rank."""
pred(g::SugiGraph, v::Int, layers::Vector{Vector{Int}}) = layers[g.verts[v].rank][g.verts[v].pos - 1]

function place_block!(g::SugiGraph, layers::Vector{Vector{Int}}, root::Int, x::Dict{Int,Float64})
    haskey(x, root) && return
    x[root] = 0.0
    w = root
    while true
        if g.verts[w].pos > 1
            u = g.verts[pred(g, w, layers)].root
            place_block!(g, layers, u, x)
            g.verts[root].sink == root && (g.verts[root].sink = g.verts[u].sink)
            if g.verts[root].sink == g.verts[u].sink
                gap = (g.verts[root].block_width + g.verts[u].block_width) * 0.5
                x[root] = max(x[root], x[u] + gap)
            end
        end
        w = g.verts[w].align
        w == root && break
    end
    while g.verts[w].align != root
        w = g.verts[w].align
        x[w] = x[root]
        g.verts[w].sink = g.verts[root].sink
    end
end

function place_blocks(g::SugiGraph, layers::Vector{Vector{Int}})
    x = Dict{Int,Float64}()
    for root in 1:_nv(g)
        g.verts[root].root == root && place_block!(g, layers, root, x)
    end
    return x
end

function do_horizontal_compaction!(g::SugiGraph, layers::Vector{Vector{Int}})
    compute_block_max_vertex_widths!(g)
    x = place_blocks(g, layers)

    nlayers = length(layers)
    for i in 1:nlayers
        v = layers[i][1]
        g.verts[v].sink == v || continue
        vsink = g.verts[v].sink
        g.verts[vsink].shift == Inf && (g.verts[vsink].shift = 0.0)

        j, k = i, 1
        while true
            v = layers[j][k]
            while g.verts[v].align != g.verts[v].root
                v = g.verts[v].align
                j += 1
                if g.verts[v].pos > 1
                    u = pred(g, v, layers)
                    gap = (g.verts[v].block_width + g.verts[u].block_width) * 0.5
                    distance_v_u = x[v] - (x[u] + gap)
                    u_sink = g.verts[u].sink
                    g.verts[u_sink].shift = min(g.verts[u_sink].shift,
                                                 g.verts[g.verts[v].sink].shift + distance_v_u)
                end
            end
            k = g.verts[v].pos + 1
            (k > length(layers[j]) || g.verts[v].sink != g.verts[layers[j][k]].sink) && break
        end
    end

    for v in 1:_nv(g)
        x[v] = x[v] + g.verts[g.verts[v].sink].shift
    end
    return x
end

"""
    create_layouts(g, layers)

Run the 4-direction Brandes & Köpf alignment and return the resulting list
of 4 coordinate maps (`Dict{Int,Float64}`, vertex id => x coordinate).
"""
function create_layouts(g::SugiGraph, layers::Vector{Vector{Int}})
    reset_alignment!(g, layers)
    mark_type1_conflicts!(g, layers)

    layouts = Dict{Int,Float64}[]
    cur_layers = [copy(l) for l in layers]
    for _vdir in 1:2
        for hdir in 1:2
            reset_alignment!(g, cur_layers)
            create_vertical_alignments!(g, cur_layers)
            layout = do_horizontal_compaction!(g, cur_layers)
            if hdir == 2 # :left
                for k in keys(layout)
                    layout[k] = -layout[k]
                end
            end
            push!(layouts, layout)
            for row in cur_layers
                reverse!(row)
            end
        end
        reverse_graph!(g)
        reverse!(cur_layers)
    end
    reset_alignment!(g, layers)
    return layouts
end

function align_to_smallest_width_layout!(layouts::Vector{Dict{Int,Float64}})
    isempty(layouts) && return layouts
    min_max = map(layouts) do c
        lo, hi = extrema(values(c))
        (lo, hi, hi - lo)
    end
    _, min_width = findmin(x -> x[3], min_max)

    for (i, layout) in enumerate(layouts)
        shift = isodd(i) ? min_max[i][1] - min_max[min_width][1] :
                min_max[min_width][2] - min_max[i][2]
        for k in keys(layout)
            layout[k] += shift
        end
    end
    return layouts
end

"""Average the two median x-coordinates (of the 4 alignment directions) for
each vertex — "the average median is both order and separation preserving"
(Brandes & Köpf, 2002)."""
function calculate_relative_coords(layouts::Vector{Dict{Int,Float64}})
    coords = Dict{Int,Float64}()
    for k in keys(layouts[1])
        v = sort!([layouts[1][k], layouts[2][k], layouts[3][k], layouts[4][k]])
        coords[k] = (v[2] + v[3]) / 2.0
    end
    return coords
end
