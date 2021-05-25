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
using GraphMakie, LightGraphs
import Random; Random.seed!(2) # hide
nothing #hide
```

# Basic Usage & Algorithms
All of the algorithms follow the [Layout Interface](@ref). Each layout algorithm
is represented by a type `Algorithm <: AbstractLayout`. The parameters of each
algorithm can be set with keyword arguments. The `Algorithm` object itself is
callable and transforms the adjacency matrix and returns a list of `Point{N,T}` from [`GeometryBasics.jl`](https://github.com/JuliaGeometry/GeometryBasics.jl).

```
alg = Algorithm(; p1="foo", p2=:bar)
positions = alg(adj_matrix)
```

## Scalable Force Directed Placement
```@docs
SFDP
```
### Example
```@example layouts
using NetworkLayout: SFDP
g = wheel_graph(10)
layout = SFDP(tol=0.01, C=0.2, K=1)
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f
```

### Iterator Example
```@example layouts
iterator = LayoutIterator(layout, adjacency_matrix(g))
record(f, "sfdp_animation.mp4", iterator; framerate = 10) do pos
    p[:node_positions][] = pos
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
using NetworkLayout: Buchheim

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
using NetworkLayout: Spring
g = smallgraph(:cubical)
layout = Spring()
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```
### Iterator Example
```@example layouts
iterator = LayoutIterator(layout, adjacency_matrix(g))
record(f, "spring_animation.mp4", iterator; framerate = 10) do pos
    p[:node_positions][] = pos
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
using NetworkLayout: Stress
g = complete_graph(10)
layout = Stress()
Random.seed!(1)
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```

### Iterator Example
```@example layouts
iterator = LayoutIterator(layout, adjacency_matrix(g))
record(f, "stress_animation.mp4", iterator; framerate = 100) do pos
    p[:node_positions][] = pos
    autolimits!(ax)
end
nothing #hide
```
![stress animation](stress_animation.mp4)

##  Circular Layout
```@docs
Circular
```
```@example layouts
using NetworkLayout: Circular
g = smallgraph(:karate)
layout = Circular()
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```

##  Shell Layout
```@docs
Shell
```
```@example layouts
using NetworkLayout: Shell
g = smallgraph(:petersen)
layout = Shell(nlist=[6:10,])
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```

## Spectral
```@docs
Spectral
```
Spectral needs 3d which isn't ready yet on the GraphMakie side.
```
using JSServe
Page(exportable=true, offline=true)
using WGLMakie #hide
WGLMakie.activate!() # hide
set_theme!(resolution=(800, 400)) # hide
scatter([1,2,3], [1,2,3])
```
```
g = smallgraph(:petersen)
layout = Shell(nlist=[6:10,])
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f #hide
```
