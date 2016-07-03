export layout_fdp

dist(xj,xi,yj,yi) = sqrt((xj-xi)^2 + (yj-yi)^2)
mag(a,b) = a^2 + b^2
f_attr(xi,xj,yi,yj) = (dist(xi,xj,yi,yj)^2) / K
f_attr(xi,xj,yi,yj) = -C*(K^2) / dist(xi,xj,yi,yj)

function layout_fdp(g, tol)
  converged = false
  step = 1
  energy = typemax(Float64)
  while !converged
    x0 = x
    y0 = y
    energy0 = energy
    energy = 0
    for i in 1:N
      fx = 0.0
      fy = 0.0
      for j in 1:N
        if g[i,j] == 1
          fx = fx + f_attr(x[i],x[j],y[i],y[j]) * ((x[j] - x[i]) / dist(x[j],x[i],y[j],y[i]))
          fy = fy + f_attr(x[i],x[j],y[i],y[j]) * ((y[j] - y[i]) / dist(x[j],x[i],y[j],y[i]))
        else
          fx = fx + f_repln(x[i],x[j],y[i],y[j]) * ((x[j] - x[i]) / dist(x[j],x[i],y[j],y[i]))
          fy = fy + f_repln(x[i],x[j],y[i],y[j]) * ((y[j] - y[i]) / dist(x[j],x[i],y[j],y[i]))
        end
      end
        x[i] = x[i] + step * (fx / sqrt(mag(fx,fy)))
        y[i] = y[i] + step * (fy / sqrt(mag(fx,fy)))
        energy = energy + mag(fx,fy)
    end
    step = update_step(step, energy, energy0)
    converged = dist_tolerance(x,x0,y,y0,tol)
  end
  return x, y
end

function update_step(step, energy, energy0)
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
end

function dist_tolerance(x,x0,y,y0,K)
  xt = x-x0
  yt = y-y0
  const N = size(x,1)
  for i in 1:N
    if sqrt(mag(xt[i],yt[i])) > K
      return false
    end
  end
  return true
end
