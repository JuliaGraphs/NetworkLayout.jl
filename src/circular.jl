module Circular

using GeometryTypes
export layout

"""
This function wrap from [NetworkX](https://github.com/networkx/networkx)
Position nodes on a circle.
**Parameters**
*G*
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
function layout(G)
    if size(G,1) == 1
        return Point{2,Float64}[Point(0.0,0.0)]
    else
        # Discard the extra angle since it matches 0 radians.
        θ = linspace(0, 2pi, size(G,1) + 1)[1:end-1]
        return Point{2,Float64}[(cos(o), sin(o)) for o in θ]
    end
end

end # end of module
