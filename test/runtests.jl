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

function jagmesh()
    jagmesh_path = joinpath(dirname(@__FILE__), "jagmesh1.mtx")
    array = round(Int, open(readdlm, jagmesh_path))
    row = array[:,1]
    col = array[:,2]
    entry = [(1:3600)...]
    sparse(row,col,entry)
end
jagmesh_adj = jagmesh()

@testset "Testing NetworkLayout" begin

    @testset "Testing SFDP" begin
        println("SFDP")
        @testset "Testing Jagmesh1 graph" begin
            positions = @time SFDP.layout(jagmesh_adj, 2, tol=0.9, K=1)
            @test typeof(positions) == Vector{Point{2, Float64}}
            positions = @time SFDP.layout(jagmesh_adj, 3, tol=0.9, K=1)
            @test typeof(positions) == Vector{Point{3, Float64}}
        end

        @testset "Testing WheelGraph" begin
            println("Wheelgraph")
            g = WheelGraph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time SFDP.layout(adj_matrix, 2, tol=0.1, K=1)
            @test typeof(positions) == Vector{Point{2, Float64}}
            positions = @time SFDP.layout(adj_matrix, 3, tol=0.1, K=1)
            @test typeof(positions) == Vector{Point{3, Float64}}
        end

    end

    @testset "Testing Stress Majorization Algorithm" begin
    println("Stress")

    @testset "Testing Jagmesh1 graph" begin
        positions = @time Stress.layout(jagmesh_adj, Point2{Float64}, iterations=10)
        @test typeof(positions) == Vector{Point{2, Float64}}
        positions = @time Stress.layout(jagmesh_adj, Point3f0, iterations=10)
        @test typeof(positions) == Vector{Point{3, Float32}}
    end

    @testset "Testing WheelGraph" begin
        println("Stress WheelGraph")
        g = WheelGraph(10)
        adj_matrix = adjacency_matrix(g)
        positions = @time Stress.layout(adj_matrix, Point2f0, iterations=10)
        @test typeof(positions) == Vector{Point2f0}
        positions = @time Stress.layout(adj_matrix, Point{3, Float64}, iterations=10)
        @test typeof(positions) == Vector{Point{3, Float64}}
    end

    end

    @testset "Testing Spring Algorithm" begin
        println("Spring WheelGraph")
        @testset "Testing WheelGraph" begin
            g = WheelGraph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Spring.layout(adj_matrix, 2, C=2.0, MAXITER=100, INITTEMP=2.0)
            @test typeof(positions) == Vector{Point{2, Float64}}
            positions = @time Spring.layout(adj_matrix, 3, C=2.0, MAXITER=100, INITTEMP=2.0)
            @test typeof(positions) == Vector{Point{3, Float64}}
        end

    end

    @testset "Testing Spectral Algorithm" begin
        println("Spectral WheelGraph")
        @testset "Testing WheelGraph" begin
          g = WheelGraph(10)
          adj_matrix = adjacency_matrix(g)
          positions = @time Spectral.layout(adj_matrix)
          @test typeof(positions) == Vector{Point{3, Float64}}
        end
    end

    @testset "Testing Circular Layout Algorithm" begin
        println("Circular WheelGraph")
        @testset "Testing WheelGraph" begin
          g = WheelGraph(10)
          adj_matrix = adjacency_matrix(g)
          positions = @time Circular.layout(adj_matrix)
          @test typeof(positions) == Vector{Point{2, Float64}}
        end

    end

    @testset "Testing Shell Layout Algorithm" begin
        println("Shell WheelGraph")

        @testset "Testing WheelGraph" begin
            g = WheelGraph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Shell.layout(adj_matrix)
            @test typeof(positions) == Vector{Point{2, Float64}}
        end
    end

    @testset "Testing Buchheim Tree Drawing" begin
        println("Buchheim Buchheim")

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
            locs = @time Buchheim.layout(adj_list,nodesize)
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
          locs = @time Buchheim.layout(n)
          @test typeof(locs) == Vector{Point{2, Float64}}
        end
    end

end
