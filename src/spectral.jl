# -----------------------------------------------------
# -----------------------------------------------------

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
# also: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf

# this recipe uses the technique of Spectral Graph Drawing, which is an
# under-appreciated method of graph layouts; easier, simpler, and faster
# than the more common spring-based methods.

module Spectral

using GeometryTypes
using LinearAlgebra: diag, eigen, Diagonal

function make_symmetric(adj_matrix::AbstractMatrix)
    adj_matrix = copy(adj_matrix)
    for i=1:size(adj_matrix,1), j=i+1:size(adj_matrix,2)
        adj_matrix[i,j] = adj_matrix[j,i] = adj_matrix[i,j]+adj_matrix[j,i]
    end
    adj_matrix
end

function compute_laplacian(adj_matrix, node_weights)
    n, m = size(adj_matrix)
    # @show size(adj_matrix), size(node_weights)
    @assert n == m == length(node_weights)

    # scale the edge values by the product of node_weights, so that "heavier" nodes also form
    # stronger connections
    adj_matrix = adj_matrix .* sqrt.(node_weights * node_weights')

    # D is a diagonal matrix with the degrees (total weights for that node) on the diagonal
    deg = vec(sum(adj_matrix, dims=1)) - diag(adj_matrix)
    D = Matrix(Diagonal(deg))
    T = eltype(node_weights)
    # Laplacian (L = D - adj_matrix)
    L = T[i == j ? deg[i] : -adj_matrix[i,j] for i=1:n,j=1:n]

    L, D
end

function layout(
            adj_matrix::M;
            node_weights=ones(eltype(M),
            size(adj_matrix,1)),
            kw_args...
        ) where {M<:AbstractMatrix}
    layout!(adj_matrix,node_weights,kw_args...)
end

function layout!(adj_matrix, node_weights, kw_args...)
    adj_matrix = make_symmetric(adj_matrix)
    L, D = compute_laplacian(adj_matrix, node_weights)

    # get the matrix of eigenvectors
    v = eigen(L, D).vectors
    # x, y, and z are the 2nd through 4th eigenvectors of the solution to the
    # generalized eigenvalue problem Lv = λDv
    Point{3, Float64}[(v[2, i], v[3, i], v[4, i]) for i in 1:size(v,2)]
end

end # end of module
