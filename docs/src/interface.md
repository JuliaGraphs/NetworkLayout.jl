```@meta
CurrentModule = NetworkLayout
```

# Layout Interface

At its core, each layout algorithm is a mapping 
```
graph ↦ node_positions
```
where each layout has several parameters. The main goal of the following interface is to keep the separation between parameters and function call. 
Each layout is implemented as subtype of [`AbstractLayout`](@ref).

```@docs
AbstractLayout
```

Therefore, each `Layout <: AbstractLayout` is a functor and can be passed around as a function `graph ↦ node_positions` which encapsulates all the parameters. This is handy for plotting libraries such as [GraphMakie.jl](http://juliaplots.org/GraphMakie.jl/previews/PR9/).

There are some additional guidelines:
- All of the parameters should be keyword arguments, i.e. it should be allways
  possible to call `Layout()` without specifying any parameters.
- Algorithms should allways return `Vector{Point{Dim,Ptype}}`. If the type or
  dimensions can be altered use the keywords `dim` and `Ptype` for it.
- Some parameters may depend on the specific network (i.e. length of start
  positions vector). If possible, there should be a fallback option (i.e.
  truncate the list of start positions if network is to small or append with
  random values).

## Iterative Layouts
Iterative layouts are a specific type of layouts which produce a sequence of positions rather than a single list of positions. Those algorithms are implemented as subtypes of [`IterativeLayout`](@ref):

```@docs
IterativeLayout
```

One can instantiate an iterable object [`LayoutIterator`](@ref)
```@docs
LayoutIterator
```
