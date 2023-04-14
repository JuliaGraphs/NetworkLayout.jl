export SquareGrid, squaregrid

"""
    SquareGrid(; kwargs...)(adj_matrix)
    squaregrid(adj_matrix; kwargs...)

Position nodes on a 2 dimensional rectagular grid. The nodes are placed in order
from upper left to lower right. To skip positions see `skip` argument.

Takes adjacency matrix representation of a network and returns coordinates of
the nodes.

## Keyword Arguments
- `Ptype=Float64`: Determines the output type `Point{2,Ptype}`
- `cols=:auto`: Columns of the grid, the rows are determined automatic. If `:auto` the layout will be square-ish.
- `dx=Ptype(1), dy=Ptype(-1)`: Ofsets between rows/cols.
- `skip=Tuple{Int,Int}[]`: Specify positions to skip when placing nodes.
  `skip=[(i,j)]` means to keep the position in the `i`-th row and `j`-th column
  empty.
"""
@addcall struct SquareGrid{Ptype,CT} <: AbstractLayout{2,Ptype}
    cols::CT
    dx::Ptype
    dy::Ptype
    skip::Vector{Tuple{Int,Int}}
end

function SquareGrid(; Ptype=Float64, cols=:auto, dx=Ptype(1), dy=Ptype(-1), skip=Tuple{Int,Int}[])
    return SquareGrid{Ptype,typeof(cols)}(cols, dx, dy, skip)
end

function layout(algo::SquareGrid{Ptype}, adj_matrix::AbstractMatrix) where {Ptype}
    N = assertsquare(adj_matrix)
    M = N + length(algo.skip)

    if algo.cols === :auto
        cols = (M == 0) ? 0 : (isqrt(M - 1) + 1)
    else
        cols = algo.cols
    end

    positions = Vector{Point2{Ptype}}(undef, N)

    n = 1
    for j in 1:typemax(Int), i in 1:cols
        if (j, i) ∉ algo.skip
            positions[n] = Point2{Ptype}((i - 1) * algo.dx, (j - 1) * algo.dy)
            n += 1
            n > N && break
        end
    end

    return positions
end
