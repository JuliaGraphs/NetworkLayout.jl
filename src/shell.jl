"""
This function is copy from [IainNZ](https://github.com/IainNZ)'s [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl)
Position nodes in concentric circles.
**Parameters**
*G*
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
function shell_layout(G, nlist::Union{Void, Vector{Vector{Int}}} = nothing)
    if _nv(G) == 1
        return [0.0], [0.0]
    end
    if nlist == nothing
        nlist = Array(Vector{Int}, 1)
        nlist[1] = collect(1:_nv(G))
    end
    radius = 0.0
    if length(nlist[1]) > 1
        radius = 1.0
    end
    locs_x = Float64[]
    locs_y = Float64[]
    for nodes in nlist
        # Discard the extra angle since it matches 0 radians.
        θ = linspace(0, 2pi, length(nodes) + 1)[1:end-1]
        append!(locs_x, radius*cos(θ))
        append!(locs_y, radius*sin(θ))
        radius += 1.0
    end
    locs_x, locs_y
end
