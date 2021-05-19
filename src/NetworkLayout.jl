module NetworkLayout

export LayoutIterator, layout

using GeometryBasics
using LinearAlgebra: norm

abstract type AbstractLayout{Dim,Ptype} end

dim(::AbstractLayout{Dim,Ptype}) where {Dim,Ptype} = Dim
ptype(::AbstractLayout{Dim,Ptype}) where {Dim,Ptype} = Ptype

(lay::AbstractLayout)(adj_matrix) = layout(lay, adj_matrix)

abstract type IterativeLayout{Dim,Ptype} <: AbstractLayout{Dim,Ptype} end

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


# XXX TODO FIXME incredible hacky solution, bootstrapping problem betwen GraphMakie and NetworkLayout
# this is needed right now to crate the docs...
function Base.getproperty(x::Type, s::Symbol)
    if s===:layout && x == getfield(NetworkLayout, :Spring)
        return Spring()
    else
        return getfield(x, s)
    end
end

end
