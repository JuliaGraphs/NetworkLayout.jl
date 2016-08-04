# NetworkLayout.jl
Layout algorithms for graphs and trees in pure Julia.

[![Coverage Status](https://coveralls.io/repos/github/JuliaGraphs/NetworkLayout.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGraphs/NetworkLayout.jl?branch=master)

Linux, OSX : [![Build Status](https://travis-ci.org/JuliaGraphs/NetworkLayout.jl.svg?branch=master)](https://travis-ci.org/JuliaGraphs/NetworkLayout.jl)

Windows : [![Build status](https://ci.appveyor.com/api/projects/status/328ph0ct3t8fc91u/branch/master?svg=true)](https://ci.appveyor.com/project/abhijithanilkumar/networklayout-jl-b6gcd/branch/master)

## Algorithms

### Scalable Force Directed Placement

Spring-Electric Force Directed Placement algorithm as explained in [Efficient and High Quality Force-Directed Graph Drawing](http://yifanhu.net/PUB/graph_draw_small.pdf) by Yifan Hu.

Module Name : `SFDP`

#### Usage

```julia
layout(adjacency_matrix,dimension,intial_position;tolerance,C,K,iterations)
```
##### arguments
  * `adjacency_matrix` - sparse/full adjacency matrix that represents the graph
  * `dimension` - dimension in which the layouting code has to be generated
  * `initial_position` - co-ordinates of the layout to start with. By default, a random layout is used
  * `tolerance` - permitted distance between current and calculated co-ordinate. Lower the tolerance, more the number of iterations (kwarg)
  * `C, K` - used to scale the layout (kwarg)
  * `iterations` - Number of iterations we apply the forces (kwarg)

##### returns
  `network` - co-ordinates of nodes in the layout

##### iterator

A user can move between iterations using a `Layout` object.


#### Example

```julia
using LightGraphs
using NetworkLayout:SFDP
g = WheelGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a,2,tol=0.1,C=1,K=1,iterations=10) # generate 2D layout
```
Using Iterator :

```julia
g = WheelGraph(30)
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

### Buchheim Tree Drawing

Buchheim Tree Drawing as explained in [Improving Walker's Algorithm to Run in Linear Time](http://dirk.jivas.de/papers/buchheim02improving.pdf) by Christoph Buchheim, Michael Junger and Sebastian Leipert.

Module Name : `Buchheim`

#### Usage

```julia
layout(adjacency_list,nodesize)
```

##### arguments
 * `adjacency_list` - adjacency list that represents the tree
 * `nodesize` - sizes of nodes (used to position the nodes)

##### returns
 * `locs` - co-ordinates of the layout

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
 locs = layout(adj_list,nodesize) # generating the layout for the tree
 ```

### Spring/Repulsion Model

Spring/Repulsion model of Fruchterman and Reingold (1991). Original code taken from [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl)

Module Name : `Spring`

#### Usage

```julia
layout(adjacency_matrix,dimension,intial_position;C,iterations,INITTEMP)
```
##### arguments
 * `adjacency_matrix` - sparse/full adjacency matrix that represents the graph
 * `dimension` - dimension in which the layouting code has to be generated
 * `initial_position` - co-ordinates of the layout to start with. By default, a random layout is used
 * `iterations` - Number of iterations we apply the forces (kwarg)
 * `C` - Constant to fiddle with density of resulting layout (kwarg)
 * `INITTEMP` - Initial "temperature", controls movement per iteration (kwarg)

##### returns
 `network` - co-ordinates of nodes in the layout

##### iterator

A user can move between iterations using a `Layout` object.


#### Example

```julia
using LightGraphs
using NetworkLayout:Spring
g = WheelGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a,2,C=2.0,iterations=100,K=2.0) # generate 2D layout
```
Using Iterator :

```julia
g = WheelGraph(30)
a = adjacency_matrix(g)
iterations = 200
C = 2.0
INITTEMP = 2.0
network = Layout(a,locs,C,iterations,INITTEMP)
state = start(network)
while !done(network,state)
 network, state = next(network,state)
end
return network.positions
```

### Stress Majorization

Based on the algorithm explained in "Graph Drawing by Stress Majorization" by Emden R Gansner, Yehuda Koren and Stephen North. Original code taken from [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl)

Module Name : `Stress`

#### Usage

```julia
layout(δ,dimension;weights,intial_position,iterations,abstols,reltols,abstolx)
```
##### arguments
 * `δ` - Matrix of pairwise distances (Adjacency Matrix can be used)
 * `dimension` - dimension in which the layouting code has to be generated
 * `weights` - Matrix of weights (kwarg)
 * `initial_position` - co-ordinates of the layout to start with. By default, a random layout is used (kwarg)
 * `iterations` - Number of iterations we apply the forces (kwarg)
 * `abstols` - Absolute tolerance for convergence of stress (kwarg)
 * `reltols` - Relative tolerance for convergence of stress (kwarg)
 * `abstolx` - Absolute tolerance for convergence of layout (kwarg)

##### returns
 `network` - co-ordinates of nodes in the layout

##### iterator

A user can move between iterations using a `Layout` object.


#### Example

```julia
using LightGraphs
using NetworkLayout:Stress
g = WheelGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a,2) # generate 2D layout
```
Using Iterator :

```julia
g = WheelGraph(30)
δ = adjacency_matrix(g)
startpositions=rand(Point{p, Float64}, size(δ,1))
iter = Layout(δ, Point{N,T}; startpositions=startpositions)
state = start(iter)
while !done(iter, state)
    iter, state = next(iter, state)
end
iter.positions
```

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
 `network` - co-ordinates of nodes in the layout

#### Example

```julia
using LightGraphs
using NetworkLayout:Spectral
g = WheelGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a) # generate 3D layout
```

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
 `network` - co-ordinates of nodes in the layout

#### Example

```julia
using LightGraphs
using NetworkLayout:Circular
g = WheelGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a) # generate 2D layout
```

### Shell Layout Algorithm

Position nodes in concentric circles. Original code taken from [GraphPlot.jl](https://github.com/afternone/GraphPlot.jl)

Module Name : `Shell`

#### Usage

```julia
layout(adjacency_matrix)
```
##### arguments
 * `adjacency_matrix` - Adjacency Matrix in dense/sparse format

##### returns
 `network` - co-ordinates of nodes in the layout

#### Example

```julia
using LightGraphs
using NetworkLayout:Shell
g = WheelGraph(30)
a = adjacency_matrix(g) # generates a sparse adjacency matrix
network = layout(a) # generate 2D layout
```
