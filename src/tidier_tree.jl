
function layout_tree_tr()

end

function setup(tree, depth=0)
  if length(tree.children) == 0
    tree.x = 0
    tree.y = depth
    return tree
  end
  if length(tree.children) == 1
    tree.x = setup(tree.children[0], depth+1).x
    return tree
  end
  left = setup(tree.children[0], depth+1)
  right = setup(tree.children[1], depth+1)
  tree.x = fix_subtrees(left, right)
  return tree
end
