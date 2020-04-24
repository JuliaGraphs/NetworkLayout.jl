"""
Using the Spring-Electric model suggested by Yifan Hu
(http://yifanhu.net/PUB/graph_draw_small.pdf)
Forces are calculated as :
        f_attr(i,j) = ||xi - xj||^2 / K ,     i<->j
        f_repln(i,j) = -CK^2 / ||xi - xj|| ,  i!=j
Arguments :
  adj_matrix      Sparse/Full Adjacency matrix of the graph
  tol             Tolerance distance - Minimum distance between 2 nodes
  C, K            Constants that help scale the layout
Output :
  positions       Co-ordinates for the nodes
"""
module ParallelSFDP

using Distributed
using GeometryTypes
using LinearAlgebra: norm
using SharedArrays

struct Layout{M<:SharedArray, P<:SharedArray, T<:AbstractFloat}
    adj_matrix::M
    positions::P
    tol::T
    C::T
    K::T
    iterations::Int
end

function Layout(
        adj_matrix, PT::Type{Point{N, T}}=Point{2, Float64};
        startpositions=(2*rand(PT, size(adj_matrix,1)) .- 1),
        tol=1.0, C=0.2, K=1.0, iterations=100
    ) where {N, T}
    Layout(SharedArray(Matrix(adj_matrix)), SharedArray(startpositions), T(tol), T(C), T(K), Int(iterations))
end

layout(adj_matrix, dim::Int; kw_args...) = layout(adj_matrix, Point{dim,Float64}; kw_args...)

function layout(
        adj_matrix, PT::Type{Point{N, T}}=Point{2, Float64};
        startpositions = (2*rand(PT, size(adj_matrix,1)) .- 1),
        kw_args...
    ) where {N, T}
    layout!(adj_matrix,startpositions;kw_args...)
end

function layout!(adj_matrix,
         startpositions::AbstractVector{Point{N, T}};
         kw_args...
    ) where {N, T}
    network = Layout(adj_matrix, Point{N,T}; startpositions=startpositions, kw_args...)
    next = iterate(network)
    while next != nothing
        (i, state) = next
        next = iterate(network, state)
    end
    return network.positions
end

# TODO this iterator is a bit strange, it looks as if it network.iterations == 1
# then it could take a long time. Furthermore I'm not sure if the iterators are
# the best solution in this situation

function iterate(network::Layout{M, P, T}) where {M, P, T}
    return network, (one(T), typemax(T), 0, true, 1, copy(network.positions))
end

function iterate(network::Layout, state)
    step, energy, progress, start, iter, locs0 = state
    K, C, tol, adj_matrix = network.K, network.C, network.tol, network.adj_matrix
    locs = network.positions; locs0 = copy(locs)
    energy0 = energy; energy = zero(energy)
    energy_vec = zeros(nprocs())
    energy_vec = SharedArray(energy_vec)
    @sync begin
        for p in procs(adj_matrix)
            @async remotecall_wait(compute_locs!, p, adj_matrix, locs, energy_vec, step, C, K, p)
        end
    end
    energy = sum(energy_vec)
    step, progress = update_step(step, energy, energy0, progress)

    if iter == network.iterations && !start ||
       dist_tolerance(network.positions, locs0, network.K, network.tol)
        return nothing
    end

    return network, (step, energy, progress, false, iter+1, locs0)
end

# Calculate Attractive force
f_attr(a, b, K) = (norm(a-b)^2) / K
# Calculate Repulsive force
f_repln(a, b, C, K) = -C*(K^2) / norm(a-b)

function loop_range(adj_matrix)
    nchunks = length(procs(adj_matrix))
    idx = indexpids(adj_matrix)
    if idx == 0
        return 1, 1
    end
    splits = [round(Int, s) for s in range(0, stop=size(adj_matrix,1), length=nchunks+1)]
    return splits[idx]+1,splits[idx+1]
end

function compute_locs!(adj_matrix, locs, energy_vec, step, C, K, p)
    F = eltype(locs); N = size(adj_matrix,1)
    start_i, end_i = loop_range(adj_matrix)
    for i in start_i:end_i
        force = F(0)
        for j in 1:N
            i == j && continue
            if adj_matrix[i,j] == 1
            # Attractive forces for adjacent nodes
                force += F(f_attr(locs[i],locs[j],K) * ((locs[j] - locs[i]) / norm(locs[j] - locs[i])))
            else
            # Repulsive forces
                force += F(f_repln(locs[i],locs[j],C,K) * ((locs[j] - locs[i]) / norm(locs[j] - locs[i])))
            end
        end
        locs[i] = locs[i] + step * (force / norm(force))
        energy_vec[p-1] = energy_vec[p-1] + norm(force)^2
    end
end

function update_step(step, energy::T, energy0, progress) where {T}
    # cooldown step
    t = T(0.9)
    if energy < energy0
        progress = progress + 1
        if progress >= 5
            progress = 0
            step = step/t
        end
    else
        progress = 0
        step = t * step
    end
    return step, progress
end

function dist_tolerance(locs,locs0, K, tol)
    # check whether the layout is optimal
    for i in 1:size(locs, 1)
        if norm(locs[i]-locs0[i]) >= K*tol
            return false
        end
    end
    return true
end

end #end of module
