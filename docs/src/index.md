```@meta
CurrentModule = NetworkLayout
```

# NetworkLayout
This is the Documentation for [NetworkLayout](https://github.com/JuliaGraphs/NetworkLayout.jl).

All example images on this page are created using [`Makie.jl`](https://github.com/JuliaPlots/Makie.jl) and the `graphplot` recipe from [`GraphMakie.jl`](https://github.com/JuliaPlots/GraphMakie.jl).

```@example layouts
using CairoMakie
CairoMakie.activate!(type="png") # hide
set_theme!(resolution=(800, 400)) #hide
CairoMakie.inline!(true) # hide
using NetworkLayout
using GraphMakie, Graphs
nothing #hide
```

# Basic Usage & Algorithms
All of the algorithms follow the [Layout Interface](@ref). Each layout algorithm
is represented by a type `LayoutAlgorithm <: AbstractLayout`. The parameters of each
layout can be set with keyword arguments. The `LayoutAlgorithm` object itself is
callable and transforms the adjacency matrix and returns a list of `Point{N,T}` from [`GeometryBasics.jl`](https://github.com/JuliaGeometry/GeometryBasics.jl).

```
alg = LayoutAlgorithm(; p1="foo", p2=:bar)
positions = alg(adj_matrix)
```
Each of the layouts comes with a lowercase function version:
```
positions = layoutalgorithm(adj_matrix; p1="foo", b2=:bar)
```

Instead of using the adjacency matrix you can use `AbstractGraph` types from [`Graphs.jl`](https://github.com/JuliaGraphs/Graphs.jl) directly.
```
g = complete_graph(10)
positions = layoutalgorithm(g)
```

## Scalable Force Directed Placement
```@docs
SFDP
```
### Example
```@example layouts
g = wheel_graph(10)
layout = SFDP(Ptype=Float32, tol=0.01, C=0.2, K=1)
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f
```

### Iterator Example
```@example layouts
iterator = LayoutIterator(layout, g)
record(f, "sfdp_animation.mp4", iterator; framerate = 10) do pos
    p[:node_pos][] = pos
    autolimits!(ax)
end
nothing #hide
```
![sfdp animation](sfdp_animation.mp4)

## Buchheim Tree Drawing
```@docs
Buchheim
```
### Example
```@example layouts
adj_matrix = [0 1 1 0 0 0 0 0 0 0;
              0 0 0 0 1 1 0 0 0 0;
              0 0 0 1 0 0 1 0 1 0;
              0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 1 0 1;
              0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0 0 0]
g = SimpleDiGraph(adj_matrix)
layout = Buchheim()
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```

## Spring/Repulsion Model
```@docs
Spring
```
### Example
```@example layouts
g = smallgraph(:cubical)
layout = Spring(Ptype=Float32)
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```
### Iterator Example
```@example layouts
iterator = LayoutIterator(layout, g)
record(f, "spring_animation.mp4", iterator; framerate = 10) do pos
    p[:node_pos][] = pos
    autolimits!(ax)
end
nothing #hide
```
![spring animation](spring_animation.mp4)

## Stress Majorization
```@docs
Stress
```
### Example
```@example layouts
g = SimpleGraph(936)
for l in eachline(joinpath(@__DIR__,"..","..","test","jagmesh1.mtx"))
    s = split(l, " ")
    src, dst = parse(Int, s[1]), parse(Int, s[2])
    src != dst && add_edge!(g, src, dst)
end

layout = Stress(Ptype=Float32)
f, ax, p = graphplot(g; layout=layout, node_size=3, edge_width=1)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```

### Iterator Example
```@example layouts
iterator = LayoutIterator(layout, g)
record(f, "stress_animation.mp4", iterator; framerate = 7) do pos
    p[:node_pos][] = pos
    autolimits!(ax)
end
nothing #hide
```
![stress animation](stress_animation.mp4)

##  Shell/Circular Layout
```@docs
Shell
```
```@example layouts
g = smallgraph(:petersen)
layout = Shell(nlist=[6:10,])
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```

## SquareGrid Layout
```@docs
SquareGrid
```
```@example layouts
g = Grid((12,4))
layout = SquareGrid(cols=12)
f, ax, p = graphplot(g, layout=layout, nlabels=repr.(1:nv(g)), nlabels_textsize=10, nlabels_distance=5)
ylims!(-4.5,.5); hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```

## Spectral Layout
```@docs
Spectral
```
```@example layouts
g = watts_strogatz(1000, 5, 0.03; seed=5)
layout = Spectral(dim=2)
f, ax, p = graphplot(g, layout=layout, node_size=0.0, edge_width=1.0)
hidedecorations!(ax); hidespines!(ax); f #hide
f #hide
```
```@example layouts
set_theme!(resolution=(800, 800)) #hide
using Random; Random.seed!(5) # hide
layout = Spectral()
f, ax, p = graphplot(g, layout=layout, node_size=0.0, edge_width=1.0)
f #hide
```
