# NetworkLayout.jl
Layout algorithms for graphs and trees in pure Julia.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagraphs.org/NetworkLayout.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliagraphs.org/NetworkLayout.jl/dev/)
[![Build Status](https://github.com/JuliaGraphs/NetworkLayout.jl/workflows/CI/badge.svg)](https://github.com/JuliaGraphs/NetworkLayout.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGraphs/NetworkLayout.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGraphs/NetworkLayout.jl)

## Installation
``` julia
pkg> add NetworkLayout
```
## Algorithms
The available algorithms and their parameters can be found in the
[docs](https://juliagraphs.org/NetworkLayout.jl/stable).

All of the algorithms represent mappings `adjacency matrix ↦ vector of
positions` where the positions are represented by the `Point` datatype from
[`GeometryBasics.jl](https://github.com/JuliaGeometry/GeometryBasics.jl)

``` julia
using NetworkLayout
using Graphs

adj_matrix = adjacency_matrix(wheel_graph(10))

pos = spring(adj_matrix; iterations=20)
pos = algorithm(adj_matrix)
```
There is also a "delayed" functor version of each algorithm:
```julia
layout = Spring(; iterations=20)
pos = layout(adj_matrix)
```
Instead of passing a adjacency matrix on can also pass `Graphs.jl` graphs directly.
