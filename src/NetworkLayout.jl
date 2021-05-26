module NetworkLayout

export LayoutIterator, layout

using GeometryBasics
using LinearAlgebra: norm

"""
    AbstractLayout{Dim,Ptype}

Abstract supertype for all layouts. Each layout `Algorithm <: AbstractLayout` needs to
implement

    layout(algo::Algorithm, adj_matrix)::Vector{Point{Dim,Ptype}}

which takes the adjacency matrix representation of a network and returns a list of
node positions. Each `Algorithm` object holds all of the necessary parameters.

By implementing `layout` the algorithm also inherits the function-like property

    Algorithm(; kwargs...)(adj_matrix) -> node_positions
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

function layout(alg::IterativeLayout, adj_matrix)
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

include("sfdp.jl")
include("buchheim.jl")
include("spring.jl")
include("stress.jl")
include("spectral.jl")
include("circular.jl")
include("shell.jl")

end
