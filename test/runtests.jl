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

end
