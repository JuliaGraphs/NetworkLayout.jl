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
struct SFDP{Dim,Ptype,T<:AbstractFloat} <: IterativeLayout{Dim,Ptype}
    tol::T
    C::T
    K::T
    iterations::Int
    initialpos::Vector{Point{Dim,Ptype}}
end

function SFDP(; dim=2, Ptype=Float64, tol=1.0, C=0.2, K=1.0, iterations=100, initialpos=Point{dim,Ptype}[])
    if !isempty(initialpos)
        initialpos = Point.(initialpos)
        Ptype = eltype(eltype(initialpos))
        # TODO fix initial pos if list has points of multiple types
        Ptype == Any && error("Please provide list of Point{N,T} with same T")
        dim = length(eltype(initialpos))
    end
    return SFDP{dim,Ptype,typeof(tol)}(tol, C, K, iterations, initialpos)
end

function Base.iterate(iter::LayoutIterator{SFDP{Dim,Ptype,T}}) where {Dim,Ptype,T}
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    N = size(adj_matrix, 1)
    M = length(algo.initialpos)
    startpos = Vector{Point{Dim,Ptype}}(undef, N)
    # take the first
    for i in 1:min(N, M)
        startpos[i] = algo.initialpos[i]
    end
    # fill the rest with random points
    for i in (M + 1):N
        startpos[i] = 2 .* rand(Point{Dim,Ptype}) .- 1
    end
    # iteratorstate: (#iter, energy, step, progress, old pos, stopflag)
    return startpos, (1, typemax(T), one(T), 0, startpos, false)
end

function Base.iterate(iter::LayoutIterator{<:SFDP}, state)
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    iter, energy0, step, progress, locs0, stopflag = state
    K, C, tol = algo.K, algo.C, algo.tol

    # stop if stopflag (tol reached) or nr of iterations reached
    if iter >= algo.iterations || stopflag
        return nothing
    end

    locs = copy(locs0)
    energy = zero(energy0)
    Ftype = eltype(locs)
    N = size(adj_matrix, 1)
    for i in 1:N
        force = zero(Ftype)
        for j in 1:N
            i == j && continue
            if adj_matrix[i, j] == 1
                # Attractive forces for adjacent nodes
                force += Ftype(f_attr(locs[i], locs[j], K) .*
                               ((locs[j] .- locs[i]) / norm(locs[j] .- locs[i])))
            else
                # Repulsive forces
                force += Ftype(f_repln(locs[i], locs[j], C, K) .*
                               ((locs[j] .- locs[i]) / norm(locs[j] .- locs[i])))
            end
        end
        locs[i] = locs[i] .+ step .* (force ./ norm(force))
        energy = energy + norm(force)^2
    end
    step, progress = update_step(step, energy, energy0, progress)

    # if the tolerance is reached set stopflag to keep claculated point but stop next iteration
    if dist_tolerance(locs, locs0, K, tol)
        stopflag = true
    end

    return locs, (iter+1, energy, step, progress, locs, stopflag)
end

# Calculate Attractive force
f_attr(a, b, K) = (norm(a .- b) .^ 2) ./ K
# Calculate Repulsive force
f_repln(a, b, C, K) = -C .* (K^2) / norm(a .- b)

function update_step(step, energy::T, energy0, progress) where {T}
    # cooldown step
    t = T(0.9)
    if energy < energy0
        progress = progress + 1
        if progress >= 5
            progress = 0
            step = step / t
        end
    else
        progress = 0
        step = t * step
    end
    return step, progress
end

function dist_tolerance(locs, locs0, K, tol)
    # check whether the layout is optimal
    for i in 1:size(locs, 1)
        if norm(locs[i] .- locs0[i]) >= K * tol
            return false
        end
    end
    return true
end
