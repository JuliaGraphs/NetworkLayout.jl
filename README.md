# NetworkLayout.jl
Layout algorithms for graphs and trees in pure Julia.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagraphs.org/NetworkLayout.jl/stable)
[![Build Status](https://github.com/JuliaGraphs/NetworkLayout.jl/workflows/CI/badge.svg)](https://github.com/JuliaGraphs/NetworkLayout.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGraphs/NetworkLayout.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGraphs/NetworkLayout.jl)

## Installation
``` julia
pkg> add NetworkLayout.jl
```
## Algorithms
The available algorithms and their parameters can be found in the
[docs](https://juliagraphs.org/NetworkLayout.jl/stable).


All of the algorithms represent mappings `adjacency matrix ↦ vector of
positions` where the positions are represented by the `Point` datatype from
[`GeometryBasics.jl](https://github.com/JuliaGeometry/GeometryBasics.jl)

``` julia
using NetworkLayout
using LightGraphs

adj_matrix = adjacency_matrix(wheel_graph(10))

algorithm = NetworkLayout.Spring(; iterations=20)
pos = algorithm(adj_matrix)
```
