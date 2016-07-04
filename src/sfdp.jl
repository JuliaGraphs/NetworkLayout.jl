export layout_fdp

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

# Function for calculating dist b/w vectors
dist(xj,xi,yj,yi) = sqrt((xj-xi)^2 + (yj-yi)^2)
# Function to calculate magnitude of the vector
mag(a,b) = a^2 + b^2
# Calculate Attractive force
f_attr(xi,xj,yi,yj,K) = (dist(xi,xj,yi,yj)^2) / K
# Calculate Repulsive force
f_repln(xi,xj,yi,yj,C,K) = -C*(K^2) / dist(xi,xj,yi,yj)

function layout_fdp(g; tol=1, C=0.2, K=1)
  converged = false
  step = 1
  progress = 0
  energy = typemax(Float64)
  N = size(g,1)
  x = 2*rand(N) .- 1
  y = 2*rand(N) .- 1
  while !converged
    x0 = copy(x)
    y0 = copy(y)
    energy0 = energy
    energy = 0
    for i in 1:N
      fx = 0.0
      fy = 0.0
      for j in 1:N
        i == j && continue
        if g[i,j] == 1
          # Attractive forces for adjacent nodes
          fx = fx + f_attr(x[i],x[j],y[i],y[j],K) * ((x[j] - x[i]) / dist(x[j],x[i],y[j],y[i]))
          fy = fy + f_attr(x[i],x[j],y[i],y[j],K) * ((y[j] - y[i]) / dist(x[j],x[i],y[j],y[i]))
        else
          # Repulsive forces
          fx = fx + f_repln(x[i],x[j],y[i],y[j],C,K) * ((x[j] - x[i]) / dist(x[j],x[i],y[j],y[i]))
          fy = fy + f_repln(x[i],x[j],y[i],y[j],C,K) * ((y[j] - y[i]) / dist(x[j],x[i],y[j],y[i]))
        end
      end
      x[i] = x[i] + step * (fx / sqrt(mag(fx,fy)))
      y[i] = y[i] + step * (fy / sqrt(mag(fx,fy)))
      energy = energy + mag(fx,fy)
    end
    step, progress = update_step(step, energy, energy0, progress)
    converged = dist_tolerance(x,x0,y,y0,K,tol)
  end
  return x, y
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

function dist_tolerance(x,x0,y,y0,K,tol)
  # check whether the layout is optimal
  xt = x-x0
  yt = y-y0
  const N = size(x,1)
  for i in 1:N
    if sqrt(mag(xt[i],yt[i])) >= K*tol
      return false
    end
  end
  return true
end
