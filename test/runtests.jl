using NetworkLayout:SFDP
using NetworkLayout:Spring
using NetworkLayout:Stress
using NetworkLayout:Buchheim
using NetworkLayout:Spectral
using NetworkLayout:Shell
using NetworkLayout:Circular
using LightGraphs
using BaseTestNext
using GeometryTypes

const jagmesh = joinpath(dirname(@__FILE__), "jagmesh1.mtx")

@testset "Testing NetworkLayout" begin

  @testset "Testing SFDP" begin

    @testset "Testing Jagmesh1 graph" begin
      array = round(Int, open(readdlm, jagmesh))
      row = array[:,1]
      col = array[:,2]
      entry = [1 for i in 1:3600]
      adj_matrix = sparse(row,col,entry)
      positions = SFDP.layout(adj_matrix, 2, tol=0.1, K=1)
      @test typeof(positions) == Vector{Point{2, Float64}}
      positions = SFDP.layout(adj_matrix, 3, tol=0.1, K=1)
      @test typeof(positions) == Vector{Point{3, Float64}}
    end

    @testset "Testing WheelGraph" begin
      g = WheelGraph(10)
      adj_matrix = adjacency_matrix(g)
      positions = SFDP.layout(adj_matrix, 2, tol=0.1, K=1)
      @test typeof(positions) == Vector{Point{2, Float64}}
      positions = SFDP.layout(adj_matrix, 3, tol=0.1, K=1)
      @test typeof(positions) == Vector{Point{3, Float64}}
    end

  end

  @testset "Testing Stress Majorization Algorithm" begin

    @testset "Testing Jagmesh1 graph" begin
      array = round(Int,open(readdlm, jagmesh))
      row = array[:,1]
      col = array[:,2]
      entry = [1 for i in 1:3600]
      adj_matrix = sparse(row,col,entry)
      positions = Stress.layout(adj_matrix, 2)
      @test typeof(positions) == Vector{Point{2, Float64}}
      positions = Stress.layout(adj_matrix, 3)
      @test typeof(positions) == Vector{Point{3, Float64}}
    end

    @testset "Testing WheelGraph" begin
      g = WheelGraph(10)
      adj_matrix = adjacency_matrix(g)
      positions = Stress.layout(adj_matrix, 2)
      @test typeof(positions) == Vector{Point{2, Float64}}
      positions = Stress.layout(adj_matrix, 3)
      @test typeof(positions) == Vector{Point{3, Float64}}
    end

  end

  @testset "Testing Spring Algorithm" begin

    @testset "Testing WheelGraph" begin
      g = WheelGraph(10)
      adj_matrix = adjacency_matrix(g)
      positions = Spring.layout(adj_matrix, 2, C=2.0, MAXITER=100, INITTEMP=2.0)
      @test typeof(positions) == Vector{Point{2, Float64}}
      positions = Spring.layout(adj_matrix, 3, C=2.0, MAXITER=100, INITTEMP=2.0)
      @test typeof(positions) == Vector{Point{3, Float64}}
    end

  end

  @testset "Testing Spectral Algorithm" begin

    @testset "Testing WheelGraph" begin
      g = WheelGraph(10)
      adj_matrix = adjacency_matrix(g)
      positions = Spectral.layout(adj_matrix)
      @test typeof(positions) == Vector{Point{3, Float64}}
    end

  end

  @testset "Testing Circular Layout Algorithm" begin

    @testset "Testing WheelGraph" begin
      g = WheelGraph(10)
      adj_matrix = adjacency_matrix(g)
      positions = Circular.layout(adj_matrix)
      @test typeof(positions) == Vector{Point{2, Float64}}
    end

  end

  @testset "Testing Shell Layout Algorithm" begin

    @testset "Testing WheelGraph" begin
      g = WheelGraph(10)
      adj_matrix = adjacency_matrix(g)
      positions = Shell.layout(adj_matrix)
      @test typeof(positions) == Vector{Point{2, Float64}}
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
      locs = Buchheim.layout(adj_list,nodesize)
      @test typeof(locs) == Vector{Point{2, Float64}}
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
      locs = Buchheim.layout(n)
      @test typeof(locs) == Vector{Point{2, Float64}}
    end

  end

end
