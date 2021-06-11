using LinearAlgebra: diag, eigen, Diagonal

export Spectral, spectral

"""
    Spectral(; kwargs...)(adj_matrix)
    spectral(adj_matrix; kwargs...)

This algorithm uses the technique of Spectral Graph Drawing, which is an
under-appreciated method of graph layouts; easier, simpler, and faster than the
more common spring-based methods. For reference see [Yehuda Koren,
2002](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf).

Takes adjacency matrix representation of a network and returns coordinates of
the nodes.

## Keyword Arguments
- `dim=3`, `Ptype=Float64`: Determines dimension and output type `Point{dim,Ptype}`.
- `nodeweights=Float64[]`: Vector of weights. If network size does not match the
  length of `nodesize` use `ones` instead.
"""
@addcall struct Spectral{Dim,Ptype,FT<:AbstractFloat} <: AbstractLayout{Dim,Ptype}
    nodeweights::Vector{FT}
end

Spectral(; dim=3, Ptype=Float64, nodeweights=Float64[]) = Spectral{dim,Ptype,eltype(nodeweights)}(nodeweights)

function make_symmetric(adj_matrix::AbstractMatrix)
    adj_matrix = copy(adj_matrix)
    for i in 1:size(adj_matrix, 1), j in (i + 1):size(adj_matrix, 2)
        adj_matrix[i, j] = adj_matrix[j, i] = adj_matrix[i, j] + adj_matrix[j, i]
    end
    return adj_matrix
end

function compute_laplacian(adj_matrix, node_weights)
    n, m = size(adj_matrix)
    # @show size(adj_matrix), size(node_weights)
    @assert n == m == length(node_weights)

    # scale the edge values by the product of node_weights, so that "heavier" nodes also form
    # stronger connections
    adj_matrix = adj_matrix .* sqrt.(node_weights * node_weights')

    # D is a diagonal matrix with the degrees (total weights for that node) on the diagonal
    deg = vec(sum(adj_matrix; dims=1)) - diag(adj_matrix)
    D = Matrix(Diagonal(deg))
    T = eltype(node_weights)
    # Laplacian (L = D - adj_matrix)
    L = T[i == j ? deg[i] : -adj_matrix[i, j] for i in 1:n, j in 1:n]
    return L, D
end

function layout(algo::Spectral{Dim,Ptype,FT}, adj_matrix::AbstractMatrix) where {Dim,Ptype,FT}
    N = assertsquare(adj_matrix)
    # try to use user provided nodeweights
    nodeweights = if length(algo.nodeweights) == N
        algo.nodeweights
    else
        ones(FT, N)
    end

    adj_matrix = make_symmetric(adj_matrix)
    L, D = compute_laplacian(adj_matrix, nodeweights)
    # get the matrix of eigenvectors
    v = eigen(L, D).vectors
    # x, y, and z are the 2nd through 4th eigenvectors of the solution to the
    # generalized eigenvalue problem Lv = λDv
    return map(x -> Point{Dim,Ptype}(x[2:(Dim + 1)]...), eachrow(v))
end
