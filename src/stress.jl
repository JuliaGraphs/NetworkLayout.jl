



function compute_Lω(g,k)
  centres = kcentres(g,k)
  d = distances()
  m = inv(d)
  uv = svdfact(m)
  s = sign(u'*m*u)
  dw = diag()
  return (centres)
end

function update_anchors(g,centers,x)
  y = copy(x)
  for i in centers
    sum = 0
    for j in centers
      if i == j
        sum += x[j]
      end
    end
    y[i] = (1/centers)*sum
  end
  for i in centers
    sum1 = sum2 = 0
    for j in centers
      if i != j
        sum1 += w[i,j]*(x[j]+d[i,j]*(x[i]-x[j])*inv(norm(x[i]-x[j])))
        sum2 += w[i,j]
      end
    end
    y[i] = sum1/sum2
  end
  return y
end

function findrhs(x)
  N = size(x,1)
  rhs = copy(x)
  for i in 1:N
    sum = 0
    for j in 1:N
      if i!=j
        sum += (x[i]-x[j])/norm(x[i]-x[j])
      end
    end
    rhs[i] = sum
  end
  return rhs
end


function dist_tolerance(Xnew,X,tol)
  # check whether the layout is optimal
  const N = size(Xnew,1)
  for i in 1:N
    if norm(Xnew[i]-X[i]) >= K*tol
      return false
    end
  end
  return true
end

function approximate_layout(g,k)
  Lω,centers = compute_Lω(g,k)
  Xnew = (2*rand(Point{dim, Float64}, size(g,1)) .- 1)
  converged = false
  while !converged
    X = copy(Xnew)
    rhs = findrhs(X)
    Xnew = solve_Lω(Lω,rhs)
    Xnew = update_anchors(g,centers,Xnew)
    converged = dist_tolerance(Xnew,X,tol)
  end
  return Xnew
end
