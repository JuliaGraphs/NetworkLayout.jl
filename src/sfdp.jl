using GeometryTypes
export layout_fdp, Layout

"""
Using the Spring-Electric model suggested by Yifan Hu
(http://yifanhu.net/PUB/graph_draw_small.pdf)
Forces are calculated as :
        f_attr(i,j) = ||xi - xj||^2 / K ,     i<->j
        f_repln(i,j) = -CK^2 / ||xi - xj|| ,  i!=j
Arguments :
  g      Sparse/Full Adjacency matrix of the graph
  tol    Tolerance distance - Minimum distance between 2 nodes
  C, K   Constants that help scale the layout
Output :
  x, y   Co-ordinates for the nodes
"""

# Calculate Attractive force
f_attr(a,b,K) = (norm(a-b)^2) / K
# Calculate Repulsive force
f_repln(a,b,C,K) = -C*(K^2) / norm(a-b)

immutable Layout{A, P, F}
  adj_matrix::A
  positions::P
  force::F
  tol
  C
  K
end

function layout_fdp{T}(g::T, dim::Int, locs = (2*rand(Point{dim, Float64}, size(g,1)) .- 1); tol=1.0, C=0.2, K=1.0)
  N = size(g,1)
  P = eltype(locs)
  force = P(0)
  network = Layout(g,locs,force,tol,C,K)
  converged = false
  step = 1
  energy = typemax(Float64)
  progress = 0
  while !converged
    converged,step,energy,progress = layout_iterator!(network,converged,step,energy,progress)
  end
  return network.positions
end

function Base.next(network::Layout, converged)
  step = 1
  progress = 0
  energy = typemax(Float64)
  converged, step, energy, progress = layout_iterator!(network,converged,step,energy,progress)
  return network, converged
end

function layout_iterator!{A,P,F}(network::Layout{A,P,F},converged,step,energy,progress)
  g = network.adj_matrix
  N = size(g,1)
  locs = network.positions
  force = network.force
  locs0 = copy(locs)
  energy0 = energy
  energy = 0
  K = network.K
  C = network.C
  tol = network.tol
  for i in 1:N
    force = F(0)
    for j in 1:N
      i == j && continue
      if g[i,j] == 1
        # Attractive forces for adjacent nodes
        force = force + F(f_attr(locs[i],locs[j],K) * ((locs[j] - locs[i]) / norm(locs[j] - locs[i])))
      else
        # Repulsive forces
        force = force + F(f_repln(locs[i],locs[j],C,K) * ((locs[j] - locs[i]) / norm(locs[j] - locs[i])))
      end
    end
    locs[i] = locs[i] + step * (force / norm(force))
    energy = energy + norm(force)^2
  end
  step, progress = update_step(step, energy, energy0, progress)
  converged = dist_tolerance(locs,locs0,K,tol)
  return converged,step,energy,progress
end

function update_step(step, energy, energy0, progress)
  # cooldown step
  const t = 0.9
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

function dist_tolerance(locs,locs0,K,tol)
  # check whether the layout is optimal
  const N = size(locs,1)
  for i in 1:N
    if norm(locs[i]-locs0[i]) >= K*tol
      return false
    end
  end
  return true
end
