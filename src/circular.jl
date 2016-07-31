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
function circular_layout(G)
    if _nv(G) == 1
        return [0.0], [0.0]
    else
        # Discard the extra angle since it matches 0 radians.
        θ = linspace(0, 2pi, _nv(G) + 1)[1:end-1]
        return cos(θ), sin(θ)
    end
end
