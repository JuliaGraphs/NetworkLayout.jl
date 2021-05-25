export Circular

"""
    Circular(; kwargs...)(adj_matrix)
    layout(algo::Circular, adj_matrix)

Position nodes on a circle with radius 1.

Takes adjacency matrix representation of a network and returns coordinates of
the nodes.

## Keyword Arguments
- `Ptype=Float64`: Determines the output type `Point{2,Ptype}`.
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
