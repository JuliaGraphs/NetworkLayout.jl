module NetworkLayout

export LayoutIterator

using GeometryBasics
using LinearAlgebra: norm

abstract type AbstractLayout{Dim,Ptype} end

(lay::AbstractLayout)(adj_matrix) = layout(lay, adj_matrix)

abstract type IterativeLayout{Dim,Ptype} <: AbstractLayout{Dim,Ptype} end

mutable struct LayoutIterator{Dim,Ptype,T<:IterativeLayout{Dim,Ptype},M<:AbstractMatrix}
    algorithm::T
    adj_matrix::M
    positions::Vector{Point{Dim,Ptype}}
    function LayoutIterator(algorithm::T, matrix::M) where {Dim,Ptype,T<:AbstractLayout{Dim,Ptype},M}
        new{Dim,Ptype,T,M}(algorithm, matrix, Point{Dim,Ptype}[])
    end
end

function Base.iterate(iter::LayoutIterator)
    iter.positions, state = init(iter.algorithm, iter.adj_matrix)
    (iter.positions, state)
end

function Base.iterate(iter::LayoutIterator, state)
    next = step(iter.algorithm, iter.adj_matrix, iter.positions, state)

    next == nothing && return nothing

    iter.positions, newstate = next
    return (iter.positions, newstate)
end

function layout(alg::IterativeLayout, adj_matrix)
    iterable = LayoutIterator(alg, adj_matrix)
    for pos in iterable
    end
    return iterable.positions
end

include("sfdp.jl")
# include("buchheim.jl")
include("spring.jl")
# include("stress.jl")
# include("spectral.jl")
# include("circular.jl")
# include("shell.jl")

end
