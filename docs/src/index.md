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

## `pin` Positions in Interative Layouts
Sometimes it is desired to fix the positions of a few nodes while arranging the rest "naturally" around them.
The iterative layouts [`Stress`](@ref), [`Spring`](@ref) and [`SFDP`](@ref) allow to pin
nodes to certain positions, i.e. those node will stay fixed during the iteration.
```@example layouts
g = SimpleGraph(vcat(hcat(zeros(4,4), ones(4,4)), hcat(ones(4,4), zeros(4,4))))
nothing #hide
```

The keyword argument `pin` takes a Vector or a Dict of key - value pairs. The key has to be
the index of the node. The value can take three forms:
  - `idx => Point2(x,y) ` or `idx => (x,y)` overwrites the initial position of that vertex and pins it there,
  - `idx => true/false` pins or unpins the vertex, position is taken from `initialpos`-keyword argument or random,
  - `idx => (false, true)` allows for fine control over which coordinate to pin.

```@example layouts
initialpos = Dict(1=>Point2f(-1,0.5),
                  3=>Point2f(1,0),
                  4=>Point2f(1,0))
pin = Dict(1=>true,
           2=>(-1,-0.5),
           3=>(true, false),
           4=>(true, false))
nothing #hide
```
Example animation on how those keyword arguments effect different iterative layouts:
```@example layouts
springl = Spring(;initialpos, pin, seed=11) #2
sfdpl   = SFDP(;initialpos, pin, tol=0.0)
stressl = Stress(;initialpos, pin, reltols=0.0, abstolx=0.0, iterations=100)

f = Figure(resolution=(1200,500))
ax1 = f[1,1] = Axis(f; title="Spring")
ax2 = f[1,2] = Axis(f; title="SFDP")
ax3 = f[1,3] = Axis(f; title="SFDP")

for ax in [ax1, ax2, ax3]
    xlims!(ax,-2,2); ylims!(ax,-1.4,1.4); vlines!(ax, 1; color=:red); hidespines!(ax); hidedecorations!(ax)
end

node_color = vcat(:green, :green, :red, :red, [:black for _ in 1:4])
node_size = vcat([40 for _ in 1:4], [20 for _ in 1:4])
nlabels = vcat("1", "2", "3", "4", ["" for _ in 1:4])
nlabels_align = (:center, :center)
nlabels_color = :white
p1 = graphplot!(ax1, g; layout=springl, node_color, node_size, nlabels, nlabels_align, nlabels_color)
p2 = graphplot!(ax2, g; layout=sfdpl, node_color, node_size, nlabels, nlabels_align, nlabels_color)
p3 = graphplot!(ax3, g; layout=stressl, node_color, node_size, nlabels, nlabels_align, nlabels_color)

iterators = [LayoutIterator(l, g) for l in (springl, sfdpl, stressl)]
record(f, "pin_animation.mp4", zip(iterators...); framerate = 10) do (pos1, pos2, pos3)
    p1[:node_pos][] = pos1
    p2[:node_pos][] = pos2
    p3[:node_pos][] = pos3
end
nothing #hide
```
![pin animation](pin_animation.mp4)
