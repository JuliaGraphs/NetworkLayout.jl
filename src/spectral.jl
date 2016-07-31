
function make_symmetric(A::AMat)
    A = copy(A)
    for i=1:size(A,1), j=i+1:size(A,2)
        A[i,j] = A[j,i] = A[i,j]+A[j,i]
    end
    A
end

function compute_laplacian(adjmat::AMat, node_weights::AVec)
    n, m = size(adjmat)
    # @show size(adjmat), size(node_weights)
    @assert n == m == length(node_weights)

    # scale the edge values by the product of node_weights, so that "heavier" nodes also form
    # stronger connections
    adjmat = adjmat .* sqrt(node_weights * node_weights')

    # D is a diagonal matrix with the degrees (total weights for that node) on the diagonal
    deg = vec(sum(adjmat,1)) - diag(adjmat)
    D = diagm(deg)

    # Laplacian (L = D - adjmat)
    L = Float64[i == j ? deg[i] : -adjmat[i,j] for i=1:n,j=1:n]

    L, D
end


# -----------------------------------------------------
# -----------------------------------------------------

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
# also: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf

# this recipe uses the technique of Spectral Graph Drawing, which is an
# under-appreciated method of graph layouts; easier, simpler, and faster
# than the more common spring-based methods.
function spectral_graph(adjmat::AMat; node_weights::AVec = ones(size(adjmat,1)), kw...)
    adjmat = make_symmetric(adjmat)
    L, D = compute_laplacian(adjmat, node_weights)

    # get the matrix of eigenvectors
    v = eig(L, D)[2]

    # x, y, and z are the 2nd through 4th eigenvectors of the solution to the
    # generalized eigenvalue problem Lv = λDv
    vec(v[2,:]), vec(v[3,:]), vec(v[4,:])
end
