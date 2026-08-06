export Springify, springify

using CUDA
mv2cuda(arr::CuArray) = arr
mv2cuda(arr) = CuArray(arr)

"""
    Springify(; kwargs...)(adj_matrix)
    springify(adj_matrix; kwargs...)

Modification of the `Spring` algorithm. It does a little bit more separation.

- Attractive force:  `f_a(d) =  gamma * d₁₂ * d` 
- Repulsive force:   `f_r(d) = -alpha * d₁₂ / d`

where `d` is distance between two vertices, d₁₂ is the x or y axis relative 
distance of the vertices compared to `d`

Takes adjacency matrix representation of a network and returns coordinates of
the nodes.

## Keyword Arguments
- `iterations=100`: number of iterations
- `alpha=0.1f0`: Repulsive force constant
- `gamma=0.01f0`: Attractive force constant
- `initialpos=Point{2,Float32}[]`

  Provide list of initial positions. If length does not match Network size the
  initial positions will be dropped and inited with random values between [-1,1]
  in every coordinate.
"""

@addcall mutable struct Springify <: AbstractLayout{2,Float32}
    iterations::Int
    alpha::Float32
    gamma::Float32
    initialpos::Vector{Point{2,Float32}}
end

function Springify(; alpha=0.1f0, gamma=0.01f0, iterations=200, initialpos=Point{2,Float32}[])
    return Springify(iterations, alpha, gamma, initialpos)
end


function NetworkLayout.layout(algo::Springify, adj_matrix::AbstractMatrix)
    N = assertsquare(adj_matrix)
  
    force = Array{Float32,2}(undef,N,2) |> mv2cuda 
    if N != length(algo.initialpos)
      coords = randn(Float32,N,2) 
    else
      coords = Array{Float32,2}(undef,N,2) 
      for (i,xy) in  enumerate(algo.initialpos)
        coords[i,1] = xy[1]
        coords[i,2] = xy[2]
      end
    end
      # coordc .= reshape(vcat([pts[1] for pts in plt[:node_pos][]]...,[pts[2] for pts in plt[:node_pos][]]...),N,2) 
      coord_gpu = coords |> mv2cuda 
    adj_m = Array(adj_matrix) |> mv2cuda
    alpha = algo.alpha
      gamma = algo.gamma
      nth = 128
    for _ in 1:algo.iterations
          @cuda threads=nth blocks=cld(N, nth) pull_and_push_forces(adj_m,coord_gpu,force,alpha, gamma,N)
          coord_gpu.+=force
      end
    res_coord = Array(coord_gpu)
  
    if !isempty(algo.initialpos)
      for i in 1:size(res_coord,1)
        algo.initialpos[i] = Point{2,Float32}(res_coord[i,1], res_coord[i,2])
      end
    else
      algo.initialpos = [Point{2,Float32}(res_coord[i,1], res_coord[i,2]) for i in 1:size(res_coord,1)]
    end
    return algo.initialpos
  end
  
@inline has_edge(mat,I,i) = @inbounds (mat[i,I] == 1 || mat[I,i] == 1)
pull_and_push_forces(adj_m,coord,force,alpha, gamma,N) = @inbounds begin
    I = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if I > N; return end 
    force[I,1] = 0f0
    force[I,2] = 0f0
    for i in 1:I-1
        isconn = has_edge(adj_m,i,I)
        calc_forces(isconn,coord,force, alpha, gamma, i, I, 1)
        calc_forces(isconn,coord,force, alpha, gamma, i, I, 2)
    end
    for i in I+1:N
        isconn = has_edge(adj_m,i,I)
        calc_forces(isconn,coord,force, alpha, gamma, i, I, 1)
        calc_forces(isconn,coord,force, alpha, gamma, i, I, 2)
    end
    nothing
end
@inline calc_forces(isconn, coord,force, alpha, gamma, i, I, xy) = @inbounds begin
    dd = sqrt((coord[i,1] - coord[I,1])^2 + (coord[i,2] - coord[I,2])^2)
    dx = coord[i,xy] - coord[I,xy]
    e=dx/dd
    force[I,xy] -= alpha / (1f-6 + dd)*e
    force[I,xy] += isconn==1 ? gamma*dd*e : 0f0 # -alpha / (1f-6 + dd*dd)*e # - dx*(1/dd)*1000.1f0
end
