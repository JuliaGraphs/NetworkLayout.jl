"""
This function is copy from [IainNZ](https://github.com/IainNZ)'s [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl)
Position nodes in concentric circles.
**Parameters**
*adj_matrix*
a graph
*nlist*
Vector of Vector, Vector of node Vector for each shell.
**Examples**
```
julia> g = graphfamous("karate")
julia> nlist = Array(Vector{Int}, 2)
julia> nlist[1] = [1:5]
julia> nlist[2] = [6:num_vertiecs(g)]
julia> locs_x, locs_y = shell_layout(g, nlist)
```
"""
struct Shell{Ptype} <: AbstractLayout{2,Ptype}
    nlist::Vector{Vector{Int}}
end

Shell(; Ptype=Float64, nlist=Int[]) = Shell{Ptype}(nlist)

function layout(algo::Shell{Ptype}, adj_matrix) where {Ptype}
    if size(adj_matrix, 1) == 1
        return Point{2,Float64}[Point(0.0, 0.0)]
    end

    N = size(adj_matrix, 1)
    nlist = copy(algo.nlist)

    # if the list does not contain all the nodes push missing nodes to new shell
    listed_nodes = Iterators.flatten(nlist)
    @assert allunique(listed_nodes)
    diff = setdiff(1:N, listed_nodes)
    if !isempty(diff)
        push!(nlist, diff)
    end

    radius = 0.0
    if length(nlist[1]) > 1
        radius = 1.0
    end

    T = Point{2,Ptype}
    locs = Vector{T}(undef, N)
    for nodes in nlist
        # Discard the extra angle since it matches 0 radians.
        θ = range(0; stop=2pi, length=length(nodes) + 1)[1:(end - 1)]
        x = T[(radius * cos(o), radius * sin(o)) for o in θ]
        locs[nodes] = x
        radius += 1.0
    end
    return locs
end
