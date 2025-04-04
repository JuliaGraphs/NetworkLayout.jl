using NetworkLayout
using NetworkLayout: AbstractLayout, @addcall
using Graphs
using GeometryBasics
using DelimitedFiles: readdlm
using SparseArrays: sparse
using StaticArrays
using StableRNGs
using Test
using Random
using LinearAlgebra

NetworkLayout.DEFAULT_RNG[] = StableRNG

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
            ip = [Point2f(1, 2)]
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
            @test typeof(positions) == Vector{Point2f}
            positions = @time SFDP(; dim=3, Ptype=Float32, tol=0.9, K=1, iterations=10)(jagmesh_adj)
            @test typeof(positions) == Vector{Point3f}
        end

        @testset "Testing wheel_graph" begin
            println("SFDP Wheelgraph")
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time SFDP(; dim=2, Ptype=Float32, tol=0.1, K=1)(adj_matrix)
            @test typeof(positions) == Vector{Point2f}
            @test positions == sfdp(adj_matrix; dim=2, Ptype=Float32, tol=0.1, K=1)
            positions = @time SFDP(; dim=3, Ptype=Float32, tol=0.1, K=1)(adj_matrix)
            @test typeof(positions) == Vector{Point3f}
            @test positions == sfdp(adj_matrix; dim=3, Ptype=Float32, tol=0.1, K=1)

            NetworkLayout.DEFAULT_RNG[] = StableRNG
            l = SFDP()
            @test l.rng isa StableRNG
            li = LayoutIterator(l, g)
            p1, p2 = iterate(li)[1], iterate(li)[1]
            @test p1 == p2
            NetworkLayout.DEFAULT_RNG[] = MersenneTwister
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
            ip = [Point2f(1, 2)]
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
            @test typeof(positions) == Vector{Point2f}
            positions = @time Stress(; iterations=10, dim=3, Ptype=Float32)(jagmesh_adj)
            @test typeof(positions) == Vector{Point3f}
        end

        @testset "Testing wheel_graph" begin
            println("Stress wheel_graph")
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Stress(; iterations=10, Ptype=Float32)(adj_matrix)
            @test typeof(positions) == Vector{Point2f}
            @test positions == stress(adj_matrix; iterations=10, Ptype=Float32)
            positions = @time Stress(; iterations=10, dim=3, Ptype=Float32)(adj_matrix)
            @test typeof(positions) == Vector{Point3f}
            @test positions == stress(adj_matrix; iterations=10, dim=3, Ptype=Float32)

            NetworkLayout.DEFAULT_RNG[] = StableRNG
            l = Stress()
            @test l.rng isa StableRNG
            li = LayoutIterator(l, g)
            p1, p2 = iterate(li)[1], iterate(li)[1]
            @test p1 == p2
            NetworkLayout.DEFAULT_RNG[] = MersenneTwister
        end

        @testset "test pairwise_distance" begin
            using NetworkLayout: make_symmetric!, pairwise_distance
            δ = [0 1 0 0 1;
                 1 0 0 1 0;
                 0 0 0 1 1;
                 0 1 1 0 0;
                 1 0 1 0 0]
            d = pairwise_distance(δ)

            @test d == make_symmetric!([0 1 2 2 1;
                                        0 0 2 1 2;
                                        0 0 0 1 1;
                                        0 0 0 0 2;
                                        0 0 0 0 0])

        end

        @testset "test unconnected graphs" begin
            _δ = [0 1 1;
                  1 0 1;
                  1 1 0]
            N = 10# 10 fully connected 3-node graphs
            δ = kron(Matrix(I, N, N), _δ);
            g = SimpleGraph(δ)
            pos = Stress()(g)
            @test all(norm.(pos) .< 2)
            # graphplot(g; layout=Stress())

            sgkeys = [:bull, :chvatal, :cubical, :desargues,
                      :diamond, :dodecahedral, :frucht, :heawood,
                      :house, :housex, :icosahedral, :karate, :krackhardtkite,
                      :moebiuskantor, :octahedral, :pappus, :petersen,
                      :sedgewickmaze, :tetrahedral, :truncatedcube,
                      :truncatedtetrahedron, :truncatedtetrahedron_dir, :tutte]
            gs = filter(!is_directed, smallgraph.(sgkeys))
            absdim = mapreduce(nv, +, gs)
            adj = zeros(Int, absdim, absdim);
            let i = 1
                for g in gs
                    r = i:i+nv(g)-1
                    adj[r, r] .= adjacency_matrix(g)
                    i += nv(g)
                end
            end
            g = SimpleGraph(adj)

            pos = Stress()(g)
            @test all(norm.(pos) .< 20)
            # graphplot(g; layout=Stress())
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
            ip = [Point2f(1, 2)]
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
            @test typeof(positions) == Vector{Point2f}
            @test positions == spring(adj_matrix; C=2.0, iterations=100, initialtemp=2.0, Ptype=Float32)
            positions = @time Spring(; C=2.0, iterations=100, initialtemp=2.0, Ptype=Float32, dim=3)(adj_matrix)
            @test typeof(positions) == Vector{Point3f}
            @test positions ==
                  spring(adj_matrix; C=2.0, iterations=100, initialtemp=2.0, Ptype=Float32, dim=3)

            NetworkLayout.DEFAULT_RNG[] = StableRNG
            l = Spring()
            @test l.rng isa StableRNG
            li = LayoutIterator(l, g)
            p1, p2 = iterate(li)[1], iterate(li)[1]
            @test p1 == p2
            NetworkLayout.DEFAULT_RNG[] = MersenneTwister
        end

        @testset "test single node graph" begin
            g = SimpleGraph(1)
            pos = Spring()(g)
            @test length(pos) == 1
            @test !isnan(pos[1])
        end
    end

    @testset "Testing Spectral Algorithm" begin
        println("Spectral wheel_graph")
        @testset "Testing wheel_graph" begin
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Spectral()(adj_matrix)
            @test typeof(positions) == Vector{Point{3,Float64}}
            @test positions == spectral(adj_matrix)
            positions = @time Spectral(; Ptype=Float32)(adj_matrix)
            @test typeof(positions) == Vector{Point{3,Float32}}
            @test positions == spectral(adj_matrix; Ptype=Float32)
        end
    end

    @testset "Testing Shell Layout Algorithm" begin
        println("Shell wheel_graph")

        @testset "Testing wheel_graph" begin
            g = wheel_graph(10)
            adj_matrix = adjacency_matrix(g)
            positions = @time Shell()(adj_matrix)
            @test typeof(positions) == Vector{Point{2,Float64}}
            @test positions == shell(adj_matrix)
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
            @test locs == buchheim(adj_list; nodesize)
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

        @testset "test requirements" begin
            # more than one parent
            g = SimpleDiGraph(3)
            add_edge!(g, 1, 2)
            add_edge!(g, 1, 3)
            add_edge!(g, 2, 3)
            @test_throws ArgumentError Buchheim()(g)

            # node 1 not parent
            g = SimpleDiGraph(3)
            add_edge!(g, 2, 1)
            add_edge!(g, 2, 3)
            @test_throws ArgumentError Buchheim()(g)

            # not all nodes reached
            g = SimpleDiGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            @test_throws ArgumentError Buchheim()(g)
        end
    end

    @testset "Testing Square Grid Layout" begin
        println("SquareGrid")
        @testset "Testing col length" begin
            M = adjacency_matrix(SimpleGraph(4))
            positions = SquareGrid(; Ptype=Int)(M)
            @test positions == Point2.([(0, 0), (1, 0), (0, -1), (1, -1)])
            @test positions == squaregrid(M; Ptype=Int)

            M = adjacency_matrix(SimpleGraph(3))
            positions = SquareGrid(; Ptype=Int)(M)
            @test positions == Point2.([(0, 0), (1, 0), (0, -1)])
            @test positions == squaregrid(M; Ptype=Int)

            M = adjacency_matrix(SimpleGraph(5))
            positions = SquareGrid(; Ptype=Int, cols=2)(M)
            @test positions == Point2.([(0, 0), (1, 0), (0, -1), (1, -1), (0, -2)])
            @test positions == squaregrid(M; Ptype=Int, cols=2)
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

    @testset "test Graphs glue code" begin
        println("Test Graphs glue code")
        g = complete_graph(10)
        pos = Spring()(g)
        @test pos isa Vector{Point{2,Float64}}
    end

    @testset "test assert square" begin
        using NetworkLayout: assertsquare
        M1 = rand(2, 4)
        @test_throws ArgumentError assertsquare(M1)
        M2 = rand(4, 4)
        @test assertsquare(M2) == 4
        @test_throws ArgumentError buchheim(M1)
        @test_throws ArgumentError sfdp(M1)
        @test_throws ArgumentError shell(M1)
        @test_throws ArgumentError spectral(M1)
        @test_throws ArgumentError spring(M1)
        @test_throws ArgumentError squaregrid(M1)
        @test_throws ArgumentError stress(M1)
    end

    @testset "make_symmetric" begin
        using LinearAlgebra: issymmetric
        using NetworkLayout: make_symmetric!

        M = [1 0; 0 1]
        make_symmetric!(M)
        @test issymmetric(M)

        M = [0 1; 0 0]
        make_symmetric!(M)
        @test issymmetric(M)

        M = [0 0; 1 0]
        make_symmetric!(M)
        @test issymmetric(M)

        M = [0 -1; 1 0]
        @test_throws ArgumentError make_symmetric!(M)

        M = [1  0  2;
             6  2  4;
             2  0  3]
        make_symmetric!(M)
        @test M == [1  6  2;
                    6  2  4;
                    2  4  3]
    end

    @testset "infer ptype" begin
        using NetworkLayout: infer_pointtype
        @test infer_pointtype(Point2f(1,2)) == (2, Float32)
        @test infer_pointtype(Point2(1.0,2)) == (2, Float64)
        @test infer_pointtype(Point(1,2,4)) == (3, Int64)
        @test_throws ArgumentError infer_pointtype([(1,2), (2,3.2,1)])
        @test infer_pointtype([(1,2), (2,3.2)]) == (2, Float64)

        dany = Dict(1=>Point2(1,1), 4=>[1.0,2.0], 7=>(1.0, 4.0))
        @test infer_pointtype(dany) == (2, Float64)

        dany[2] = (1,2,3)
        @test_throws ArgumentError infer_pointtype(dany)
    end

    @testset "Sanitize initialpos pin"  begin
        using NetworkLayout: _sanitize_initialpos_pin
        pos = [(0,0),(1,1),(2,2)]
        pin = []
        _pos, _pin = _sanitize_initialpos_pin(2, Float64, pos, pin)
        @test _pos == Dict(pairs(Point2f.(pos)))
        @test _pin == Dict{Int, SVector{2, Bool}}()

        pos = [(0,0),(1,1),(2,2)]
        pin = Dict(1=>(true,false), 3=>true, 2=>(3.0, 3.0))
        _pos, _pin = _sanitize_initialpos_pin(2, Float64, pos, pin)
        @test _pos == Dict(1=>Point2f(0,0), 2=>Point2f(3.0,3.0), 3=>Point2f(2.0,2.0))
        @test _pin == Dict(1=>SVector(true,false), 2=>SVector(true,true), 3=>SVector(true,true))
    end

    @testset "test pin" begin
        for algo in [sfdp, spring, stress]
            g = complete_graph(10)
            ep = algo(g; pin=[(0,0), (0,0)])
            @test ep[1] == [0,0]
            @test ep[2] == [0,0]

            ep = algo(g; initialpos=Dict(4=>(0,0,0), 5=>(1,2,3)), pin=Dict(4=>true, 5=>(true, false, true)))
            @test ep[4] == [0,0,0]
            @test ep[5][[1,3]] == [1,3]
            @test ep[5][2] != [2]
        end
    end

    @testset "Align" begin
        @addcall struct Manual{Dim, Ptype} <: AbstractLayout{Dim, Ptype}
            positions :: Vector{Point{Dim, Ptype}}
        end
        NetworkLayout.layout(algo::Manual, ::AbstractMatrix) = copy(algo.positions)
        
        g = Graph(2); add_edge!(g, 1, 2)
        pos = Align(Manual([Point2f(1, 2), Point2f(2, 3)]), 0.0)(g)
        @test all(r->abs(r[2])<1e-12, pos)
        @test norm(pos[1]-pos[2]) == norm(Point2f(1, 2)-Point2f(2, 3))

        g = Graph(3); add_edge!(g, 1, 2); add_edge!(g, 2, 3); add_edge!(g, 3, 1)
        pos = Align(Manual([Point2f(0, 4), Point2f(-1, -2), Point2f(1, -2)]), 0.0)(g)
        @test pos ≈ [Point2f(4, 0), Point2f(-2, 1), Point2f(-2, -1)]

        pos = Align(Manual([Point2f(0, 4), Point2f(-1, -2), Point2f(1, -2)]), π/2)(g)
        @test pos ≈ [Point2f(0, 4), Point2f(-1, -2), Point2f(1, -2)]
    end
end
