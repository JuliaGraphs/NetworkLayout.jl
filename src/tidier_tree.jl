export layout_tree_tr

typealias list{T} Vector{Vector{T}}

immutable Tree
  nodes
  mod
  thread
  ancestor
  prelim
  shift
  change
  x
  y
end

function typegen(tree)
  mod = zeros(length(tree))
  thread = zeros(length(tree))
  prelim = zeros(length(tree))
  shift = zeros(length(tree))
  change = zeros(length(tree))
  ancestor = [i for i in 1:length(tree)]
  nodes = copy(tree)
  x = rand(length(tree))
  y = rand(length(tree))
  t = Tree(nodes,mod,thread,ancestor,prelim,shift,change,x,y)
  return t
end

function layout_tree_buchheim(t::list)
  tree = typegen(t)
  first_walk(1,tree)
  second_walk(1,-prelim[1],tree)
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
    p = parent(v)
    if find(x->(x==v),tree[p]) > 1
      w = tree[p][1]
      prelim[v] = prelim[w] + distance
      mod[v] = prelim[v] - midpoint
    else
      prelim[v] = midpoint
    end
  end
end

function second_walk(v,m,t)
  prelim = t.prelim
  mod = t.mod
  x = t.x
  y = t.y
  x[v] = prelim[v] + m
  y[v] = level(v,t)
  for w in t.nodes[v]
    second_walk(w,m+mod[v])
  end
end

function find_ancestor(w,v,defaultAncestor,tree)
  if parent(w) == parent(v)
    return tree.ancestor[w]
  else
    return defaultAncestor
  end
end  

function execute_shifts(v,t)
  tree = t.nodes
  shift = t.shift
  change = t.change
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
