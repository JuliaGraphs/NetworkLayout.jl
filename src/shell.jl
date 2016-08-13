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

module Shell

using GeometryTypes

function layout{M<:AbstractMatrix}(adj_matrix::M; nlist::Union{Void, Vector{Vector{Int}}} = nothing)
    layout!(adj_matrix,nlist)
end

function layout!{M<:AbstractMatrix}(adj_matrix::M, nlist::Union{Void, Vector{Vector{Int}}})
    if size(adj_matrix, 1) == 1
        return Point{2,Float64}[Point(0.0,0.0)]
    end
    if nlist == nothing
        nlist = Array(Vector{Int}, 1)
        nlist[1] = collect(1:size(adj_matrix,1))
    end
    radius = 0.0
    if length(nlist[1]) > 1
        radius = 1.0
    end
    T = Point{2, Float64}
    locs = T[]
    for (i, nodes) in enumerate(nlist)
        # Discard the extra angle since it matches 0 radians.
        θ = linspace(0, 2pi, length(nodes) + 1)[1:end-1]
        x = T[(radius*cos(o), radius*sin(o)) for o in θ]
        append!(locs, x)
        radius += 1.0
    end
    locs
end

end # end of module
