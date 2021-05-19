# NetworkLayout.jl
Layout algorithms for graphs and trees in pure Julia.

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagraphs.org/NetworkLayout.jl/stable) -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliagraphs.org/NetworkLayout.jl/dev/) -->
[![Build Status](https://github.com/JuliaGraphs/NetworkLayout.jl/workflows/CI/badge.svg)](https://github.com/JuliaGraphs/NetworkLayout.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGraphs/NetworkLayout.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGraphs/NetworkLayout.jl)

## Algorithms

### Scalable Force Directed Placement

Spring-Electric Force Directed Placement algorithm as explained in [Efficient and High Quality Force-Directed Graph Drawing](http://yifanhu.net/PUB/graph_draw_small.pdf) by Yifan Hu.

Module Name : `SFDP`

#### Usage

```julia
layout(adjacency_matrix,dimension;startpostitions,tol,C,K,iterations)
```
##### arguments
  * `adjacency_matrix` - sparse/full adjacency matrix that represents the graph
  * `dimension` - dimension in which the layouting code has to be generated. `dimension` can be an integer specifying
                  the dimension or a `Point` type, eg. `Point3f0` which denotes 3D.
  * `startpositions` - co-ordinates of the layout to start with. By default, a random layout is used (kwarg)
  * `tol` - permitted distance between current and calculated co-ordinate. Lower the tolerance, more the number of iterations (kwarg)
  * `C, K` - used to scale the layout (kwarg)
  * `iterations` - Number of iterations we apply the forces (kwarg)

##### returns
  `positions` - co-ordinates of nodes in the layout

##### iterator

A user can move between iterations using a `Layout` object.


#### Example

```julia
using LightGraphs
using NetworkLayout:SFDP
g = WheelGraph(10)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a,Point2f0,tol=0.1,C=1,K=1,iterations=10) # generate 2D layout
```
Using Iterator :

```julia
g = WheelGraph(10)
a = adjacency_matrix(g)
tol = 0.1
C = 0.2
K = 1
iterations = 100
network = Layout(a,locs,tol,C,K,iterations)
state = start(network)
while !done(network,state)
  network, state = next(network,state)
end
return network.positions
```
![sfdp](https://cloud.githubusercontent.com/assets/8404278/17638280/a9671850-6106-11e6-912f-be94477f5ecd.png)

The image shows a `LightGraphs.WheelGraph(10)` object layout generated by SFDP Algorithm.

### Buchheim Tree Drawing

Buchheim Tree Drawing as explained in [Improving Walker's Algorithm to Run in Linear Time](http://dirk.jivas.de/papers/buchheim02improving.pdf) by Christoph Buchheim, Michael Junger and Sebastian Leipert.

Module Name : `Buchheim`

#### Usage

```julia
layout(adjacency_list; nodesize)
```

##### arguments
 * `adjacency_list` - adjacency list that represents the tree
 * `nodesize` - sizes of nodes (used to position the nodes) (kwarg)

##### returns
 * `positions` - co-ordinates of the layout

#### Example

```julia
using NetworkLayout:Buchheim
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
 locs = layout(adj_list,nodesize=nodesize) # generating the layout for the tree
 ```
 ![tree](https://cloud.githubusercontent.com/assets/8404278/17638844/afd280a4-610a-11e6-8fea-5c99808bd740.png)

The image shows a `LightGraphs.BinaryTree(4)` object layout by Buchheim Algorithm.

### Spring/Repulsion Model

Spring/Repulsion model of Fruchterman and Reingold (1991). Original code taken from [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl)

Module Name : `Spring`

#### Usage

```julia
layout(adjacency_matrix,dimension;startpositions,C,iterations,initialtemp)
```
##### arguments
 * `adjacency_matrix` - sparse/full adjacency matrix that represents the graph
 * `dimension` - dimension in which the layouting code has to be generated. `dimension` can be an integer specifying
                  the dimension or a `Point` type, eg. `Point3f0` which denotes 3D.
 * `startpositions` - co-ordinates of the layout to start with. By default, a random layout is used (kwarg)
 * `iterations` - Number of iterations we apply the forces (kwarg)
 * `C` - Constant to fiddle with density of resulting layout (kwarg)
 * `initialtemp` - Initial "temperature", controls movement per iteration (kwarg)

##### returns
 `positions` - co-ordinates of nodes in the layout

##### iterator

A user can move between iterations using a `Layout` object.


#### Example

```julia
using LightGraphs
using NetworkLayout:Spring
g = WheelGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a,Point2f0,C=2.0,iterations=100,K=2.0) # generate 2D layout
```
Using Iterator :

```julia
g = WheelGraph(30)
a = adjacency_matrix(g)
iterations = 200
C = 2.0
initialtemp = 2.0
network = Layout(a,locs,C,iterations,initialtemp)
state = start(network)
while !done(network,state)
 network, state = next(network,state)
end
return network.positions
```
![spring](https://cloud.githubusercontent.com/assets/8404278/17638354/1c20cc56-6107-11e6-82ed-8873431d8d33.png)

The image shows a `LightGraphs.WheelGraph(10)` object layout generated by Spring Algorithm.

### Stress Majorization

Based on the algorithm explained in "Graph Drawing by Stress Majorization" by Emden R Gansner, Yehuda Koren and Stephen North. Original code taken from [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl)

Module Name : `Stress`

#### Usage

```julia
layout(δ,dimension;startpositions,weights,iterations,abstols,reltols,abstolx)
```
##### arguments
 * `δ` - Matrix of pairwise distances (Adjacency Matrix can be used)
 * `dimension` - dimension in which the layouting code has to be generated. `dimension` can be an integer specifying
                  the dimension or a `Point` type, eg. `Point3f0` which denotes 3D.
 * `weights` - Matrix of weights (kwarg)
 * `startpositions` - co-ordinates of the layout to start with. By default, a random layout is used (kwarg)
 * `iterations` - Number of iterations we apply the forces (kwarg)
 * `abstols` - Absolute tolerance for convergence of stress (kwarg)
 * `reltols` - Relative tolerance for convergence of stress (kwarg)
 * `abstolx` - Absolute tolerance for convergence of layout (kwarg)

##### returns
 `positions` - co-ordinates of nodes in the layout

##### iterator

A user can move between iterations using a `Layout` object.


#### Example

```julia
using LightGraphs
using NetworkLayout:Stress
g = CompleteGraph(10)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a,2) # generate 2D layout
```
Using Iterator :

```julia
g = CompleteGraph(10)
δ = adjacency_matrix(g)
startpositions=rand(Point{3, Float64}, size(δ,1))
iter = Layout(δ, Point{3,Float64}; startpositions=startpositions)
state = start(iter)
while !done(iter, state)
    iter, state = next(iter, state)
end
iter.positions
```

![stress](https://cloud.githubusercontent.com/assets/8404278/17638554/5e65e26c-6108-11e6-9522-30e6fa044d26.png)

The image shows a `LightGraphs.CompleteGraph(10)` object layout using Stress Algorithm.

### Spectral Layout Algorithm

Uses the technique of Spectral Graph Drawing, which is an under-appreciated method of graph layouts; easier, simpler, and faster than the more common spring-based methods. Original code taken from [PlotRecipes.jl](https://github.com/JuliaPlots/PlotRecipes.jl)

Module Name : `Spectral`

#### Usage

```julia
layout(adjacency_matrix; node_weights, kw...)
```
##### arguments
 * `adjacency_matrix` - Adjacency Matrix in dense/sparse format
 * `node_weights` - weights for different nodes (kwarg)

##### returns
 `positions` - co-ordinates of nodes in the layout

#### Example

```julia
using LightGraphs
using NetworkLayout:Spectral
g = CompleteGraph(10)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a) # generate 3D layout
```
![spectral](https://cloud.githubusercontent.com/assets/8404278/17638718/a0b451ca-6109-11e6-9a66-fd22332b8541.png)

The image shows a `LightGraphs.CompleteGraph(10)` object layout by Spectral Algorithm.

### Circular Layout Algorithm

Position nodes on a circle. Original code taken from [GraphPlot.jl](https://github.com/afternone/GraphPlot.jl)

Module Name : `Circular`

#### Usage

```julia
layout(adjacency_matrix)
```
##### arguments
 * `adjacency_matrix` - Adjacency Matrix in dense/sparse format

##### returns
 `positions` - co-ordinates of nodes in the layout

#### Example

```julia
using LightGraphs
using NetworkLayout:Circular
g = CompleteGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a) # generate 2D layout
```

![circular](https://cloud.githubusercontent.com/assets/8404278/17638609/d8eb4428-6108-11e6-934b-f326f07cf044.png)

The image shows a `LightGraphs.CompleteGraph(10)` object layout using Circular Algorithm.

### Shell Layout Algorithm

Position nodes in concentric circles. Original code taken from [GraphPlot.jl](https://github.com/afternone/GraphPlot.jl)

Module Name : `Shell`

#### Usage

```julia
layout(adjacency_matrix;nlist)
```
##### arguments
 * `adjacency_matrix` - Adjacency Matrix in dense/sparse format
 * `nlist` - Shell-wise separation of nodes (kwarg)

##### returns
 `positions` - co-ordinates of nodes in the layout

#### Example

```julia
using LightGraphs
using NetworkLayout:Shell
g = CompleteGraph(30)
n = Array(Vector{Int},2)
n[1] = [1:15]
n[2] = [16:30]
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a,nlist=n) # generate 2D layout
```
![shell](https://cloud.githubusercontent.com/assets/8404278/17638171/efac921e-6105-11e6-9e48-33471bf3b27e.png)

This figure shows a `LightGraphs.CompleteGraph(30)` object in 2 shells.

## Benchmarks

The iterative algorithms have been benchmarked using 3 different graphs: `LightGraphs.WheelGraph(10)`, `LightGraphs.WheelGraph(100)` and `jagmesh1`. The number of iterations is fixed on 100. The following graph is obtained which shows SFDP to be the fastest in a general scenario, but Stress Algorithm is faster when the number of edges per graph is comparatively less, as in `jagmesh1`.

![bench](https://cloud.githubusercontent.com/assets/8404278/17642254/fd6f1718-615b-11e6-9a30-8c1a362aead7.png)



*NOTE* : All screenshots are generated using [NetworkViz.jl](https://github.com/abhijithanilkumar/NetworkViz.jl), [ThreeJS.jl](https://github.com/rohitvarkey/ThreeJS.jl) and [Escher.jl](https://github.com/shashi/Escher.jlhttps://github.com/rohitvarkey/ThreeJS.jl). The plot used is generated using [Gadfly.jl](https://github.com/dcjones/Gadfly.jl)
