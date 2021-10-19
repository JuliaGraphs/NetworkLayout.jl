module NetworkLayout

using GeometryBasics
using Requires
using LinearAlgebra: norm
using Random

export LayoutIterator

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
    @require LightGraphs = "093fc24a-ae57-5d10-9952-331d41423f4d" begin
        function layout(l::AbstractLayout, g::LightGraphs.AbstractGraph)
            layout(l, LightGraphs.adjacency_matrix(g))
        end
        function LayoutIterator(l::IterativeLayout, g::LightGraphs.AbstractGraph)
            LayoutIterator(l, LightGraphs.adjacency_matrix(g))
        end
    end
    @require Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6" begin
        function layout(l::AbstractLayout, g::Graphs.AbstractGraph)
            layout(l, Graphs.adjacency_matrix(g))
        end
        function LayoutIterator(l::IterativeLayout, g::Graphs.AbstractGraph)
            LayoutIterator(l, Graphs.adjacency_matrix(g))
        end
    end
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

"""
    make_symmetric!(M::AbstractMatrix)

Pairwise check [i,j] and [j,i]. If one is zero, make symmetric.
If both are different and nonzero throw ArgumentError.
"""
function make_symmetric!(A::AbstractMatrix)
    indsm, indsn = axes(A)
    for i in first(indsn):last(indsn), j in (i):last(indsn)
        if A[i,j] == A[j,i]   # allready symmetric
            continue
        elseif iszero(A[i,j])
            A[i,j] = A[j,i]
        elseif iszero(A[j,i])
            A[j,i] = A[i,j]
        else
            throw(ArgumentError("Matrix can't be symmetrized!"))
        end
    end
    return A
end

"""
    @addcall

Annotate subtypes of `AbstractLayout` to create a lowercase function call for them.

    @addcall struct MyLayout{Dim, Ptype} <: AbstractLayout{Dim, Ptype}
        para
    end

will add the function

    mylayout(g; kwargs...) = layout(MyLayout(; kwargs...), g)
"""
macro addcall(expr::Expr)
    # assert subtype of abstract layout
    @assert expr.head === :struct "Macro not used on struct!"
    typedef = expr.args[2]
    @assert typedef isa Expr &&
            typedef.head === :<: &&
            typedef.args[2] isa Expr && # supertype
            typedef.args[2].args[1] ∈ [:AbstractLayout, :IterativeLayout] "Macro musst be used on subtype of AbstractLayout"

    if typedef.args[1] isa Symbol # no type parameters
        name = typedef.args[1]
    elseif typedef.args[1] isa Expr && typedef.args[1].head === :curly && typedef.args[1].args[1] isa Symbol
        name = typedef.args[1].args[1]
    else
        throw(ArgumentError("Can't find the layout name!"))
    end

    fname = Symbol(lowercase(String(name)))
    return quote
        Base.@__doc__ $(esc(expr))
        Base.@__doc__ $(esc(fname))(g; kwargs...) = layout($name(; kwargs...), g)
    end
end

include("sfdp.jl")
include("buchheim.jl")
include("spring.jl")
include("stress.jl")
include("spectral.jl")
include("shell.jl")
include("squaregrid.jl")

end
