using NetworkLayout
using LightGraphs
using BaseTestNext

@testset "Testing NetworkLayout" begin

  @testset "Testing SFDP" begin

    @testset "Testing Jagmesh1 graph" begin
      array = round(Int,open(readdlm,"jagmesh1.mtx"))
      row = array[:,1]
      col = array[:,2]
      entry = [1 for i in 1:3600]
      adj_matrix = sparse(row,col,entry)
      x, y = layout_fdp(adj_matrix, tol=0.1, K=1)
      @test size(x) == size(y)
      @test typeof(x) == typeof(y) == Array{Float64,1}
    end

    @testset "Testing WheelGraph" begin
      g = WheelGraph(10)
      adj_matrix = adjacency_matrix(g)
      x, y = layout_fdp(adj_matrix, tol=0.1, K=1)
      @test size(x) == size(y)
      @test typeof(x) == typeof(y) == Array{Float64,1}
    end

  end

end
