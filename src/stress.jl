

function layout_stress()

end

function compute_approx(g,k)
  centres = kcentres(g,k)
  d = distances()
  m = inv(d)
  uv = svdfact(m)
  s = sign(u'*m*u)
  dw = diag()
  return (centres)
end

function update_anchors(g,centers,x)
  y = x
  for i in centers
    sum = 0
    for j in centers
      if i == j
        sum += x[j]
      end
    end
    yi = (1/centers)*sum
  end
  for i in centers
    sum1 = sum2 = 0
    for j in centers
      if i != j
        sum1 += w[i,j]*(x[j]+d[i,j]*(x[i]-x[j])*inv(norm(x[i]-x[j])))
        sum2 += w[i,j]
      end
    end
    yi = sum1/sum2
  end
  return y
end
