# NetworkLayout.jl
Layout algorithms for graphs and trees in pure Julia.

## Algorithms

### Scalable Force Directed Placement

Spring-Electric Force Directed Placement algorithm as explained in [Efficient and High Quality Force-Directed Graph Drawing](http://yifanhu.net/PUB/graph_draw_small.pdf) by Yifan Hu.

#### Usage

```julia
layout_fdp(adjacency_matrix,dimension,intial_position;tolerance,C,K)
```
##### arguments
  * `adjacency_matrix` - sparse/full adjacency matrix that represents the graph
  * `dimension` - dimension in which the layouting code has to be generated
  * `initial_position` - co-ordinates of the layout to start with. By default, a random layout is used
  * `tolerance` - permitted distance between current and calculated co-ordinate. Lower the tolerance, more the number of iterations (kwarg)
  * `C, K` - used to scale the layout (kwarg)

##### returns
  `network` - co-ordinates of nodes in the layout

##### iterator

A user can move between iterations using a `Layout` object.


#### Example

```julia
using LightGraphs
using NetworkLayout
g = WheelGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout_fdp(a,2,tol=0.1,C=1,K=1) # generate 2D layout
```

### Buchheim Tree Drawing

Buchheim Tree Drawing as explained in [Improving Walker's Algorithm to Run in Linear Time](http://dirk.jivas.de/papers/buchheim02improving.pdf) by Christoph Buchheim, Michael Junger and Sebastian Leipert.

#### Usage

```julia
layout_tree_buchheim(adjacency_list,nodesize)
```

##### arguments
 * `adjacency_list` - adjacency list that represents the tree
 * `nodesize` - sizes of nodes (used to position the nodes)

##### returns
 * `x , y` - x and y co-ordinates of the layout

#### Example

```julia
using NetworkLayout
adj_list = Vector{Int}[   # adjacency list
        [2,3,4],
        [5,6],
        [7],
        [],
        [],
        [],
        []
      ]
 nodesize = [1,2.3,1.2,2,3,1.4,0.8]
 x, y = layout_tree_buchheim(adj_list,nodesize) # generating the layout for the tree
 ```
