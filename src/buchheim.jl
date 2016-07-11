export layout_tree_buchheim

typealias list{T} Vector{Vector{T}}

immutable Tree
  nodes
  mod
  thread
  ancestor
  prelim::Array{Float64}
  shift
  change
  level
  x
  y
  distance
end

function typegen(tree)
  mod = zeros(length(tree))
  thread = zeros(length(tree))
  prelim = zeros(length(tree))
  shift = zeros(length(tree))
  change = zeros(length(tree))
  level = zeros(length(tree))
  ancestor = [i for i in 1:length(tree)]
  nodes = copy(tree)
  x = rand(length(tree))
  y = rand(length(tree))
  t = Tree(nodes,mod,thread,ancestor,prelim,shift,change,level,x,y,1.0)
  return t
end

function layout_tree_buchheim(t)
  tree = typegen(t)
  first_walk(1,tree)
  second_walk(1,-tree.prelim[1],tree)
  return tree
end

function parent(v,t)
  tree = t.nodes
  for i in 1:length(tree)
    y = find(x->(x==v),tree[i])
    if length(y)!=0
      return i
    end
  end
end

function first_walk(v,t)
  prelim = t.prelim
  mod = t.mod
  tree = t.nodes
  distance = t.distance
  if length(tree[v]) == 0
    prelim[v] = 0
  else
    defaultAncestor = tree[v][1]
    for w in tree[v]
      first_walk(w,t)
      apportion(w,defaultAncestor,t)
    end
    execute_shifts(v,t)
    midpoint = (prelim[tree[v][1]] + prelim[tree[v][end]]) / 2
    p = parent(v,t)
    if p != nothing
      index = find(x->(x==v),tree[p])[1]
    else
      index = 1
    end
    if index > 1
      w = tree[p][index-1]
      prelim[v] = prelim[w] + distance
      mod[v] = prelim[v] - midpoint
    else
      prelim[v] = midpoint
    end
  end
end

function apportion(v,defaultAncestor,t)
  tree = t.nodes
  ancestor = t.ancestor
  prelim = t.prelim
  mod = t.mod
  p = parent(v,t)
  distance = t.distance
  if p != nothing
    index = find(x->(x==v),tree[p])[1]
  else
    index = 1
  end
  if index > 1
    w = tree[p][index-1]
    v_in_right = v_out_right = v
    v_in_left = w
    v_out_left = tree[parent(v_in_right,t)][1]
    s_in_right = mod[v_in_right]
    s_out_right = mod[v_out_right]
    s_in_left = mod[v_in_left]
    s_out_left = mod[v_out_left]
    while next_right(v_in_left,t)!=0 && next_left(v_in_right,t)!=0
      v_in_left = next_right(v_in_left,t)
      v_in_right = next_left(v_in_right,t)
      v_out_left = next_left(v_out_left,t)
      v_out_right = next_right(v_out_right,t)
      ancestor[v_out_right] = v
      shift = (prelim[v_in_left] + s_in_left) - (prelim[v_in_right] + s_in_right) + distance
      if shift > 0
        move_subtree(find_ancestor(v_in_left,v,defaultAncestor),v,shift,t)
        s_in_right += shift
        s_out_right += shift
      end
      s_in_left += mod[v_in_left]
      s_in_right += mod[v_in_right]
      s_out_left += mod[v_out_left]
      s_out_right += mod[v_out_right]
    end
    if next_right(v_in_left,t)!=0 && next_right(v_out_right,t)==0
      thread[v_out_right] = next_right(v_in_left,t)
      mod[v_out_right] += s_in_left - s_out_right
    end
    if next_left(v_in_right,t)!=0 next_left(v_out_left,t)==0
      thread[v_out_left] = next_left(v_in_right,t)
      mod[v_out_left] += s_in_right - s_out_left
      defaultAncestor = v
    end
  end
end

function move_subtree(w_left,w_right,shift,t)
  change = t.change
  prelim = t.prelim
  tree = t.nodes
  shifttree = t.shift
  mod = t.mod
  p_wl = parent(w_left,t)
  p_wr = parent(w_right,t)
  n_wl = find(x->(x==w_left),tree[p_wl])[1]
  n_wr = find(x->(x==w_right),tree[p_wr])[1]
  subtrees = n_wr - n_wl
  change[w_right] += shift / subtrees
  shiftree[w_right] += shift
  change[w_left] += shift / subtrees
  prelim[w_right] += shift
  mod[w_right] += shift
end

function second_walk(v,m,t)
  t = compute_depth(t)
  prelim = t.prelim
  mod = t.mod
  x = t.x
  y = t.y
  level = t.level
  x[v] = prelim[v] + m
  y[v] = level[v]
  for w in t.nodes[v]
    second_walk(w,m+mod[v],t)
  end
end

function find_ancestor(w,v,defaultAncestor,tree)
  if parent(w,t) == parent(v,t)
    return tree.ancestor[w]
  else
    return defaultAncestor
  end
end

function execute_shifts(v,t)
  tree = t.nodes
  shift = t.shift
  change = t.change
  prelim = t.prelim
  mod = t.mod
  shiftnode = 0
  changenode = 0
  for w in reverse(tree[v])
    prelim[w] += shiftnode
    mod[w] += shiftnode
    change += change[w]
    shiftnode = shiftnode + shift[w] + changenode
  end
end

function next_left(v,t)
  tree = t.nodes
  thread = t.thread
  if length(tree[v]) != 0
    return tree[v][1]
  else
    return thread[v]
  end
end

function next_right(v,t)
  tree = t.nodes
  thread = t.thread
  if length(tree[v]) != 0
    return tree[v][end]
  else
    return thread[v]
  end
end

function compute_depth(t)
  tree = t.nodes
  level = t.level
  N = length(tree)
  for v in 1:N
    for w in tree[v]
      level[w] = level[v] + 1
    end
  end
  return t
end
