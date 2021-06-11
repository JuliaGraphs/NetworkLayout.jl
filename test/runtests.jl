using NetworkLayout
using LightGraphs
using GeometryBasics
using DelimitedFiles: readdlm
using SparseArrays: sparse
using Test

function jagmesh()
    jagmesh_path = joinpath(dirname(@__FILE__), "jagmesh1.mtx")
    array = round.(Int, open(readdlm, jagmesh_path))
    row = array[:, 1]
    col = array[:, 2]
    entry = [(1:3600)...]
    return sparse(row, col, entry)
end
jagmesh_adj = jagmesh()

@testset "Testing NetworkLayout" begin
    @testset "Testing SFDP" begin
        println("SFDP")
        @testset "SFDP construction" begin
            algo = SFDP()
            @test algo isa SFDP{2,Float64}
            algo = SFDP(; dim=3, Ptype=Int)
            @test algo isa SFDP{3,Int}
            ip = [(1, 2.0, 3.0), (1, 2.0, 3)]
            algo = SFDP(; initialpos=ip)
            @test algo isa SFDP{3,Float64}
            ip = [Point2f0(1, 2)]
            algo = SFDP(; initialpos=ip)
            @test algo isa SFDP{2,Float32}
        end

        @testset "iterator size" begin
            for l in [1, 10, 100]
                adj_matrix = adjacency_matrix(wheel_graph(10))
                it = LayoutIterator(SFDP(; iterations=l, tol=0.0), adj_matrix)
                vec = Any[]
                for p in it
                    push!(vec, p)
                end
                @test length(vec) == l
                @test length(unique!(vec)) == l
            end
        end

        @testset "Testing Jagmesh1 graph" begin
            println("SFDP Jagmesh1")
            positions = @time SFDP(; dim=2, Ptype=Float32, tol=0.9, K=1, iterations=10)(jagmesh_adj)
            @test typeof(positions) == Vector{Point2f0}
            positions = @time SFDP(; dim=3, Ptype=Float32, tol=0.9, K=1, iterations=10)(jagmesh_adj)
            @test typeof(positions) == Vector{Point3f0}
        end

        @testset "Testing wheel_graph" begin
            println("SFDP Wheelgraph")
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time SFDP(; dim=2, Ptype=Float32, tol=0.1, K=1)(adj_matrix)
            @test typeof(positions) == Vector{Point2f0}
            positions = @time SFDP(; dim=3, Ptype=Float32, tol=0.1, K=1)(adj_matrix)
            @test typeof(positions) == Vector{Point3f0}
        end
    end

    @testset "Testing Stress Majorization Algorithm" begin
        println("Stress")
        @testset "Stress construction" begin
            algo = Stress()
            @test algo isa Stress{2,Float64}
            algo = Stress(; dim=3, Ptype=Int)
            @test algo isa Stress{3,Int}
            ip = [(1, 2.0, 3.0), (1, 2.0, 3)]
            algo = Stress(; initialpos=ip)
            @test algo isa Stress{3,Float64}
            ip = [Point2f0(1, 2)]
            algo = Stress(; initialpos=ip)
            @test algo isa Stress{2,Float32}
        end

        @testset "iterator size" begin
            for l in [1, 10, 100]
                adj_matrix = adjacency_matrix(wheel_graph(10))
                it = LayoutIterator(Stress(; iterations=l, abstolx=0.0, reltols=0.0, abstols=0.0), adj_matrix)
                vec = Any[]
                for p in it
                    push!(vec, p)
                end
                @test length(vec) == l
                @test length(unique!(vec)) == l
            end
        end

        @testset "Testing Jagmesh1 graph" begin
            println("Stress Jagmesh1")
            positions = @time Stress(; iterations=10, Ptype=Float32)(jagmesh_adj)
            @test typeof(positions) == Vector{Point2f0}
            positions = @time Stress(; iterations=10, dim=3, Ptype=Float32)(jagmesh_adj)
            @test typeof(positions) == Vector{Point3f0}
        end

        @testset "Testing wheel_graph" begin
            println("Stress wheel_graph")
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Stress(; iterations=10, Ptype=Float32)(adj_matrix)
            @test typeof(positions) == Vector{Point2f0}
            positions = @time Stress(; iterations=10, dim=3, Ptype=Float32)(adj_matrix)
            @test typeof(positions) == Vector{Point3f0}
        end
    end

    @testset "Testing Spring Algorithm" begin
        println("Spring wheel_graph")
        @testset "Spring construction" begin
            algo = Spring()
            @test algo isa Spring{2,Float64}
            algo = Spring(; dim=3, Ptype=Int)
            @test algo isa Spring{3,Int}
            ip = [(1, 2.0, 3.0), (1, 2.0, 3)]
            algo = Spring(; initialpos=ip)
            @test algo isa Spring{3,Float64}
            ip = [Point2f0(1, 2)]
            algo = Spring(; initialpos=ip)
            @test algo isa Spring{2,Float32}
        end

        @testset "iterator size" begin
            for l in [1, 10, 100]
                adj_matrix = adjacency_matrix(wheel_graph(10))
                it = LayoutIterator(Spring(; iterations=l), adj_matrix)
                vec = Any[]
                for p in it
                    push!(vec, p)
                end
                @test length(vec) == l
                @test length(unique!(vec)) == l
            end
        end

        @testset "Testing wheel_graph" begin
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Spring(; C=2.0, iterations=100, initialtemp=2.0, Ptype=Float32)(adj_matrix)
            @test typeof(positions) == Vector{Point2f0}
            positions = @time Spring(; C=2.0, iterations=100, initialtemp=2.0, Ptype=Float32, dim=3)(adj_matrix)
            @test typeof(positions) == Vector{Point3f0}
        end
    end

    @testset "Testing Spectral Algorithm" begin
        println("Spectral wheel_graph")
        @testset "Testing wheel_graph" begin
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Spectral()(adj_matrix)
            @test typeof(positions) == Vector{Point{3,Float64}}
            positions = @time Spectral(; Ptype=Float32)(adj_matrix)
            @test typeof(positions) == Vector{Point{3,Float32}}
        end
    end

    @testset "Testing Shell Layout Algorithm" begin
        println("Shell wheel_graph")

        @testset "Testing wheel_graph" begin
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Shell()(adj_matrix)
            @test typeof(positions) == Vector{Point{2,Float64}}
        end
        @testset "Testing Base Case" begin
            g = Graph(1)
            adj_matrix = adjacency_matrix(g)
            positions = @time Shell()(adj_matrix)
            @test typeof(positions) == Vector{Point{2,Float64}}
        end
    end

    @testset "Testing Buchheim Tree Drawing" begin
        println("Buchheim")
        using NetworkLayout: adj_mat_to_list

        @testset "matrix -> list conversion" begin
            g = wheel_graph(5)
            mat = adjacency_matrix(g)
            list = adj_mat_to_list(mat)
            @test list == [[2, 3, 4, 5], [1, 3, 5], [1, 2, 4], [1, 3, 5], [1, 2, 4]]

            g = SimpleDiGraph(7)
            add_edge!(g, 1, 2)
            add_edge!(g, 1, 3)
            add_edge!(g, 1, 4)
            add_edge!(g, 2, 5)
            add_edge!(g, 2, 6)
            add_edge!(g, 3, 7)
            mat = adjacency_matrix(g)

            list = adj_mat_to_list(mat)
            @test list == [[2, 3, 4], [5, 6], [7], [], [], [], []]
        end

        @testset "Test a Random tree" begin
            println("Buchheim Random")
            adj_list = Vector{Int}[[2, 3, 4], [5, 6], [7], [], [], [], []]
            nodesize = [1, 2, 1.5, 3, 0.5, 1, 1]
            locs = @time Buchheim(; nodesize)(adj_list)
            @test typeof(locs) == Vector{Point{2,Float64}}
            locs = @time Buchheim(; Ptype=Float32)(adj_list)
            @test typeof(locs) == Vector{Point{2,Float32}}
        end

        @testset "Test a Binary tree" begin
            println("Buchheim binary_tree")
            g = binary_tree(10)
            dirg = SimpleDiGraph(collect(edges(g)))
            a = adjacency_matrix(dirg)
            locs = @time Buchheim()(a)
            @test typeof(locs) == Vector{Point{2,Float64}}
        end

        # This testset is to cover a piece of code
        # that was not covered anymore after fixing a bug
        @testset "Test another Binary tree" begin
            tree = [[5, 6], [6], [7, 8], [8], [], [9], [9], [], []]
            locs = @time Buchheim()(tree)
            @test typeof(locs) == Vector{Point{2,Float64}}
        end
    end

    @testset "Testing Square Grid Layout" begin
        println("SquareGrid")
        @testset "Testing col length" begin
            M = adjacency_matrix(SimpleGraph(4))
            positions = SquareGrid(; Ptype=Int)(M)
            @test positions == Point2.([(0, 0), (1, 0), (0, -1), (1, -1)])

            M = adjacency_matrix(SimpleGraph(3))
            positions = SquareGrid(; Ptype=Int)(M)
            @test positions == Point2.([(0, 0), (1, 0), (0, -1)])

            M = adjacency_matrix(SimpleGraph(5))
            positions = SquareGrid(; Ptype=Int, cols=2)(M)
            @test positions == Point2.([(0, 0), (1, 0), (0, -1), (1, -1), (0, -2)])
        end

        @testset "Testing dx,dy" begin
            M = adjacency_matrix(SimpleGraph(4))
            positions = SquareGrid(; Ptype=Int, dx=2, dy=3)(M)
            @test positions == Point2.([(0, 0), (2, 0), (0, 3), (2, 3)])
        end

        @testset "Testing skip" begin
            M = adjacency_matrix(SimpleGraph(4))
            positions = SquareGrid(; Ptype=Int, skip=[(1, 1), (2, 2)])(M)
            @test positions == Point2.([(1, 0), (2, 0), (0, -1), (2, -1)])
        end
    end

    @testset "test LightGraphs glue code" begin
        println("Test LightGraphs glue code")
        g = complete_graph(10)
        pos = Spring()(g)
        @test pos isa Vector{Point{2,Float64}}
    end

    @testset "test assert square" begin
        using NetworkLayout: assertsquare
        M1 = rand(2,4)
        @test_throws ArgumentError assertsquare(M1)
        M2 = rand(4,4)
        @test assertsquare(M2) == 4
        @test_throws ArgumentError layout(Buchheim(), M1)
        @test_throws ArgumentError layout(SFDP(), M1)
        @test_throws ArgumentError layout(Shell(), M1)
        @test_throws ArgumentError layout(Spectral(), M1)
        @test_throws ArgumentError layout(Spring(), M1)
        @test_throws ArgumentError layout(SquareGrid(), M1)
        @test_throws ArgumentError layout(Stress(), M1)
    end
end
