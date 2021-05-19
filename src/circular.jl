"""
This function wrap from [NetworkX](https://github.com/networkx/networkx)
Position nodes on a circle.
**Parameters**
*adj_matrix*
a graph
**Returns**
*locs_x, locs_y*
Locations of the nodes. Can be any units you want,
but will be normalized and centered anyway
**Examples**
```
julia> g = simple_house_graph()
julia> locs_x, locs_y = circular_layout(g)
```
"""
struct Circular{Ptype} <: AbstractLayout{2,Ptype} end

Circular(; Ptype=Float64) = Circular{Ptype}()

function layout(::Circular{Ptype}, adj_matrix) where {Ptype}
    if size(adj_matrix, 1) == 1
        return Point{2,Ptype}[Point(0.0, 0.0)]
    else
        # Discard the extra angle since it matches 0 radians.
        θ = range(0; stop=2pi, length=size(adj_matrix, 1) + 1)[1:(end - 1)]
        return Point{2,Ptype}[(cos(o), sin(o)) for o in θ]
    end
end
