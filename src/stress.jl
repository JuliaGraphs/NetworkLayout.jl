

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
