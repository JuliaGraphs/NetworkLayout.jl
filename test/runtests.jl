using NetworkLayout
using LightGraphs
using BaseTestNext
using GeometryTypes

@testset "Testing NetworkLayout" begin

  @testset "Testing SFDP" begin

    @testset "Testing Jagmesh1 graph" begin
      array = round(Int,open(readdlm,"jagmesh1.mtx"))
      row = array[:,1]
      col = array[:,2]
      entry = [1 for i in 1:3600]
      adj_matrix = sparse(row,col,entry)
      positions = layout_fdp(adj_matrix, 2, tol=0.1, K=1)
      @test typeof(positions) == Array{FixedSizeArrays.Point{2,Float64},1}
      positions = layout_fdp(adj_matrix, 3, tol=0.1, K=1)
      @test typeof(positions) == Array{FixedSizeArrays.Point{3,Float64},1}
    end

    @testset "Testing WheelGraph" begin
      g = WheelGraph(10)
      adj_matrix = adjacency_matrix(g)
      positions = layout_fdp(adj_matrix, 2, tol=0.1, K=1)
      @test typeof(positions) == Array{FixedSizeArrays.Point{2,Float64},1}
      positions = layout_fdp(adj_matrix, 3, tol=0.1, K=1)
      @test typeof(positions) == Array{FixedSizeArrays.Point{3,Float64},1}
    end

  end

  @testset "Testing Buchheim Tree Drawing" begin

    @testset "Test a Random tree" begin
      adj_list = Vector{Int}[
        [2,3,4],
        [5,6],
        [7],
        [],
        [],
        [],
        []
      ]
      nodesize = [1,2,1.5,3,0.5,1,1]
      x, y = layout_tree_buchheim(adj_list,nodesize)
      @test eltype(x) == eltype(y) == Float64
    end

    @testset "Test a Binary tree" begin
      g = BinaryTree(10)
      n = Vector{Int32}[]
      a = adjacency_matrix(g)
      for i in 1:size(a,1)
         p = Int32[]
         for e in collect(edges(g))
             if e[1] == i
                 push!(p,e[2])
             end
         end
         push!(n,p)
      end
      x, y = layout_tree_buchheim(n)
      @test eltype(x) == eltype(y) == Float64
    end

  end

end
