module NetworkLayout

using GeometryBasics
using Requires
using LinearAlgebra: norm
using Random

export LayoutIterator, layout

"""
    AbstractLayout{Dim,Ptype}

Abstract supertype for all layouts. Each layout `Layout <: AbstractLayout` needs to
implement

    layout(::Layout, adj_matrix)::Vector{Point{Dim,Ptype}}

which takes the adjacency matrix representation of a network and returns a list of
node positions. Each `Layout` object holds all of the necessary parameters.

The type parameters specify the returntype `Vector{Point{Dim,Ptype}}`:
- `Dim`: the dimensionality of the layout (i.e. 2 or 3)
- `Ptype`: the type of the returned points (i.e. `Float32` or `Float64`)

By implementing `layout` the Layout also inherits the function-like property

    Layout(; kwargs...)(adj_matrix) -> node_positions
"""
abstract type AbstractLayout{Dim,Ptype} end

dim(::AbstractLayout{Dim,Ptype}) where {Dim,Ptype} = Dim
ptype(::AbstractLayout{Dim,Ptype}) where {Dim,Ptype} = Ptype

(lay::AbstractLayout)(adj_matrix) = layout(lay, adj_matrix)

"""
    IterativeLayout{Dim,Ptype} <: AbstractLayout{Dim,Ptype}

Abstract supertype for iterative layouts. Instead of implementing `layout` directly,
subtypes `Algorithm<:IterativeLayout` need to implement the [iterator
interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration)

    Base.iterate(iter::LayoutIterator{<:Algorithm})
    Base.iterate(iter::LayoutIterator{<:Algorithm}, state)

where the iteration _item_ is a `Vector{Point{Dim,Ptype}}` and the iteration
_state_ depends on the algorithm.

By implementing the iterator interface the `Algorithm` inherits the `layout` and
function-like call

    layout(algo::Algorithm, adj_matrix) -> node_postions
    Algorithm(; kwargs...)(adj_matrix) -> node_positions
"""
abstract type IterativeLayout{Dim,Ptype} <: AbstractLayout{Dim,Ptype} end

"""
    LayoutIterator{T<:IterativeLayout,M<:AbstractMatrix}(algorithm, adj_matrix)

This type bundles an [`IterativeLayout`](@ref) with an adjacency matrix to form an
iterable object whose items are the node positions.

## Example
```
for p in LayoutIterator(Stress(), adj_matrix)
    # do stuff with positions p
end
```
"""
struct LayoutIterator{T<:IterativeLayout,M<:AbstractMatrix}
    algorithm::T
    adj_matrix::M
end

function layout(alg::IterativeLayout, adj_matrix::AbstractMatrix)
    assertsquare(adj_matrix)
    iter = LayoutIterator(alg, adj_matrix)
    next = Base.iterate(iter)
    pos = next[1]
    while next !== nothing
        (item, state) = next
        next = Base.iterate(iter, state)
        pos = next !== nothing ? item : pos
    end
    return pos
end

function __init__()
    @require LightGraphs="093fc24a-ae57-5d10-9952-331d41423f4d" layout(l::AbstractLayout, g::LightGraphs.AbstractGraph) = layout(l, LightGraphs.adjacency_matrix(g))
end

"""
    assertsquare(M)

Throws `ArgumentArror` if matrix is not square. Returns size.
"""
function assertsquare(M::AbstractMatrix)
    (a, b) = size(M)
    a != b && throw(ArgumentError("Adjecency Matrix needs to be square!"))
    return a
end

include("sfdp.jl")
include("buchheim.jl")
include("spring.jl")
include("stress.jl")
include("spectral.jl")
include("shell.jl")
include("squaregrid.jl")

end
