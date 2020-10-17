"""
Using the algorithm proposed in the paper,
"Improving Walker's Algorithm to Run in Linear Time"
by Christoph Buchheim, Michael Junger, Sebastian Leipert
(http://dirk.jivas.de/papers/buchheim02improving.pdf)

Arguments
tree    Adjacency List that represents the given tree

Returns
positions     co-ordinates of the layout
"""
module Buchheim

using GeometryBasics

struct Tree{A<:AbstractVector,P<:AbstractVector,F}
    nodes::A
    mod::F
    thread::Vector{Int}
    ancestor::Vector{Int}
    prelim::F
    shift::F
    change::F
    positions::P
    nodesize::F
end

function Tree(tree::AbstractVector, nodesize)
    len = length(tree)
    mod = zeros(len)
    thread = zeros(Int, len)
    prelim = zeros(len)
    shift = zeros(len)
    change = zeros(len)
    ancestor = collect(1:len)
    nodes = copy(tree)
    positions = zeros(Point{2,Float64}, len)
    t = Tree(nodes, mod, thread, ancestor, prelim, shift, change, positions, nodesize)
    return t
end

function layout(t::AbstractVector; nodesize=ones(length(t)))
    return layout!(t, nodesize)
end

function layout!(t::AbstractVector, nodesize)
    tree = Tree(t, nodesize)
    first_walk(1, tree)
    second_walk(1, -tree.prelim[1], 0.0, tree)
    return tree.positions
end

function parent(v, t::Tree)
    tree = t.nodes
    for i in 1:length(tree)
        y = findall(x -> (x == v), tree[i])
        if length(y) != 0
            return i
        end
    end
    return nothing
end

function first_walk(v, t::Tree)
    prelim = t.prelim
    mod = t.mod
    tree = t.nodes
    nodesize = t.nodesize
    p = parent(v, t)
    if p != nothing
        index = findall(x -> (x == v), tree[p])[1]
    else
        index = 1
    end
    if length(tree[v]) == 0
        if v != tree[p][1]
            prelim[v] = prelim[tree[p][index - 1]] + (nodesize[tree[p][index - 1]])
        else
            prelim[v] = 0
        end
    else
        defaultAncestor = tree[v][1]
        for w in tree[v]
            first_walk(w, t)
            defaultAncestor = apportion(w, defaultAncestor, t)
        end
        execute_shifts(v, t)
        midpoint = (prelim[tree[v][1]] + prelim[tree[v][end]]) / 2
        if index > 1
            w = tree[p][index - 1]
            prelim[v] = prelim[w] + (nodesize[w] + 1.0)
            mod[v] = prelim[v] - midpoint
        else
            prelim[v] = midpoint
        end
    end
end

function apportion(v::T, defaultAncestor::T, t::Tree) where {T}
    tree = t.nodes
    ancestor = t.ancestor
    prelim = t.prelim
    mod = t.mod
    thread = t.thread
    p = parent(v, t)
    nodesize = t.nodesize
    if p != nothing
        index = findall(x -> (x == v), tree[p])[1]
    else
        index = 1
    end
    if index > 1
        w = tree[p][index - 1]
        v_in_right = v_out_right = v
        v_in_left = w
        v_out_left = tree[parent(v_in_right, t)][1]
        s_in_right = mod[v_in_right]
        s_out_right = mod[v_out_right]
        s_in_left = mod[v_in_left]
        s_out_left = mod[v_out_left]
        while next_right(v_in_left, t) != 0 && next_left(v_in_right, t) != 0
            v_in_left = next_right(v_in_left, t)
            v_in_right = next_left(v_in_right, t)
            v_out_left = next_left(v_out_left, t)
            v_out_right = next_right(v_out_right, t)
            ancestor[v_out_right] = v
            shift = (prelim[v_in_left] + s_in_left) - (prelim[v_in_right] + s_in_right) +
                    (nodesize[v_in_left])
            if shift > 0
                move_subtree(find_ancestor(v_in_left, v, defaultAncestor, t), v, shift, t)
                s_in_right += shift
                s_out_right += shift
            end
            s_in_left += mod[v_in_left]
            s_in_right += mod[v_in_right]
            s_out_left += mod[v_out_left]
            s_out_right += mod[v_out_right]
        end
        if next_right(v_in_left, t) != 0 && next_right(v_out_right, t) == 0
            thread[v_out_right] = next_right(v_in_left, t)
            mod[v_out_right] += s_in_left - s_out_right
        else
            if next_left(v_in_right, t) != 0 && next_left(v_out_left, t) == 0
                thread[v_out_left] = next_left(v_in_right, t)
                mod[v_out_left] += s_in_right - s_out_left
                defaultAncestor = v
            end
        end
    end
    return defaultAncestor
end

function number(v, t::Tree)
    p = parent(v, t)
    index = findall(x -> (x == v), t.nodes[p])[1]
    return index
end

function move_subtree(w_left::T, w_right::T, shift::Float64, t::Tree) where {T}
    change = t.change
    prelim = t.prelim
    tree = t.nodes
    shifttree = t.shift
    mod = t.mod
    n_wl = number(w_left, t)
    n_wr = number(w_right, t)
    subtrees = n_wr - n_wl
    change[w_right] -= shift / subtrees
    shifttree[w_right] += shift
    change[w_left] += shift / subtrees
    prelim[w_right] += shift
    return mod[w_right] += shift
end

function second_walk(v, m::Float64, depth::Float64, t::Tree)
    prelim = t.prelim
    mod = t.mod
    positions = t.positions
    nodesize = t.nodesize
    positions[v] = Point(prelim[v] + m, -depth)
    if length(t.nodes[v]) != 0
        maxdist = maximum([nodesize[i] for i in t.nodes[v]])
    else
        maxdist = 0
    end
    for w in t.nodes[v]
        second_walk(w, m + mod[v], Float64(depth + 1 + maxdist), t)
    end
end

function find_ancestor(w::T, v::T, defaultAncestor::T, tree::Tree) where {T}
    ancestor = tree.ancestor
    if ancestor[w] in tree.nodes[parent(v, tree)]
        return ancestor[w]
    else
        return defaultAncestor
    end
end

function execute_shifts(v, t::Tree)
    tree = t.nodes
    shift = t.shift
    change = t.change
    prelim = t.prelim
    mod = t.mod
    shiftnode = 0
    changenode = 0
    for w in reverse(tree[v])
        prelim[w] += shiftnode
        mod[w] += shiftnode
        changenode += change[w]
        shiftnode += shift[w] + changenode
    end
end

function next_left(v, t::Tree)
    tree = t.nodes
    thread = t.thread
    if length(tree[v]) != 0
        return tree[v][1]
    else
        return thread[v]
    end
end

function next_right(v, t::Tree)
    tree = t.nodes
    thread = t.thread
    if length(tree[v]) != 0
        return tree[v][end]
    else
        return thread[v]
    end
end

end #end of module
