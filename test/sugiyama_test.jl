@testset "Testing Sugiyama Layout" begin
    using NetworkLayout: SugiGraph, _add_vertex!, _add_edge!, find_edge, _nv, slack
    using NetworkLayout: remove_cycles!, rank!, init_cutvalues!, init_low_lim!
    using NetworkLayout: leave_edge, enter_edge
    using NetworkLayout: insert_dummy_vertices!, bilayer_cross_count, total_crossings, SugiOrder
    using NetworkLayout: order_layer, barycenter, ordering
    using NetworkLayout: reset_alignment!, mark_type1_conflicts!, create_vertical_alignments!,
                          reverse_graph!, place_blocks

    println("Sugiyama")

    # ---- helpers: build a SugiGraph from 0-based (tail, head) edges --------
    function sugi_from_edges(edges::Vector{Tuple{Int,Int}}, n::Int)
        g = SugiGraph()
        for _ in 1:n
            _add_vertex!(g)
        end
        for (t, h) in edges
            _add_edge!(g, t + 1, h + 1)
        end
        return g
    end

    function set_tree_edges!(g, tree_edges::Vector{Tuple{Int,Int}})
        for (t, h) in tree_edges
            g.verts[t + 1].is_tree_vertex = true
            g.verts[h + 1].is_tree_vertex = true
            g.edges[find_edge(g, t + 1, h + 1)].is_tree_edge = true
        end
    end

    @testset "Sugiyama construction" begin
        algo = Sugiyama()
        @test algo isa Sugiyama{Float64}
        @test algo.ranking_type == :networksimplex
        @test algo.crossing_minimization == :barycenter
        @test algo.direction == :down

        algo = Sugiyama(; Ptype=Float32)
        @test algo isa Sugiyama{Float32}

        @test_throws ArgumentError Sugiyama(; direction=:sideways)
        @test_throws ArgumentError Sugiyama(; ranking_type=:foo)
        @test_throws ArgumentError Sugiyama(; crossing_minimization=:foo)
        @test_throws ArgumentError Sugiyama(; minimum_length=0)
    end

    # ---- Phase 1: ranking / network simplex, ported from rust-sugiyama's --
    # ---- p1_layering test fixtures (src/algorithm/p1_layering/tests.rs)  --
    @testset "ranking: network simplex (ported rust fixtures)" begin
        EXAMPLE_GRAPH = [(0, 1), (1, 2), (2, 3), (3, 7), (4, 6), (5, 6), (6, 7), (0, 4), (0, 5)]
        FEASIBLE_TREE_NEG = [(0, 1), (1, 2), (2, 3), (3, 7), (4, 6), (5, 6), (6, 7)]
        FEASIBLE_TREE_POS = [(0, 1), (0, 4), (1, 2), (2, 3), (3, 7), (4, 6), (5, 6)]
        LOW_LIM_GRAPH = [(0, 1), (1, 2), (1, 3), (0, 4), (4, 5), (5, 6), (4, 7), (4, 8)]

        @testset "cut values, tree with a negative cut value" begin
            g = sugi_from_edges(EXAMPLE_GRAPH, 8)
            set_tree_edges!(g, FEASIBLE_TREE_NEG)
            init_cutvalues!(g)
            expected = [3, 3, 3, 3, 0, 0, -1]
            for (i, (t, h)) in enumerate(FEASIBLE_TREE_NEG)
                @test g.edges[find_edge(g, t + 1, h + 1)].cut_value == expected[i]
            end
        end

        @testset "cut values, all-positive tree" begin
            g = sugi_from_edges(EXAMPLE_GRAPH, 8)
            set_tree_edges!(g, FEASIBLE_TREE_POS)
            init_cutvalues!(g)
            expected = [2, 1, 2, 2, 2, 1, 0]
            for (i, (t, h)) in enumerate(FEASIBLE_TREE_POS)
                @test g.edges[find_edge(g, t + 1, h + 1)].cut_value == expected[i]
            end
        end

        @testset "low/lim: parent structure and subtree-containment invariant" begin
            # the exact low/lim *numbers* are traversal-order dependent (an
            # implementation detail); what matters is that they correctly
            # encode subtree containment, which is what network simplex relies on.
            g = sugi_from_edges(LOW_LIM_GRAPH, 9)
            set_tree_edges!(g, LOW_LIM_GRAPH)
            init_low_lim!(g)

            expected_parent = Dict(0 => nothing, 1 => 0, 2 => 1, 3 => 1, 4 => 0, 5 => 4,
                                    6 => 5, 7 => 4, 8 => 4)
            for (id, parent) in expected_parent
                @test g.verts[id + 1].parent == (parent === nothing ? nothing : parent + 1)
            end

            function subtree_of(root, tree_edges, n)
                children = [Int[] for _ in 1:n]
                for (t, h) in tree_edges
                    push!(children[t + 1], h + 1)
                end
                seen, stack = Set([root]), [root]
                while !isempty(stack)
                    v = pop!(stack)
                    for c in children[v]
                        c in seen || (push!(seen, c); push!(stack, c))
                    end
                end
                return seen
            end

            for x in 0:8, y in 0:8
                sub = subtree_of(x + 1, LOW_LIM_GRAPH, 9)
                vx, vy = g.verts[x + 1], g.verts[y + 1]
                @test ((y + 1) in sub) == (vx.low <= vy.lim <= vx.lim)
            end
        end

        @testset "leave_edge / enter_edge" begin
            g = sugi_from_edges(EXAMPLE_GRAPH, 8)
            set_tree_edges!(g, FEASIBLE_TREE_NEG)
            init_cutvalues!(g)
            init_low_lim!(g)

            e = leave_edge(g)
            @test e !== nothing
            @test (g.edges[e].tail, g.edges[e].head) == (7, 8) # rust node 6->7

            swap = enter_edge(g, e, 1)
            @test g.edges[swap].tail == 1 # rust node 0
            @test g.edges[swap].head in (5, 6) # rust node 4 or 5

            g2 = sugi_from_edges(EXAMPLE_GRAPH, 8)
            set_tree_edges!(g2, FEASIBLE_TREE_POS)
            init_cutvalues!(g2)
            init_low_lim!(g2)
            @test leave_edge(g2) === nothing
        end

        @testset "full network simplex is optimal & consistent" begin
            function is_correct(g, minimum_length)
                tree_cv = [g.edges[eid].cut_value for eid in eachindex(g.edges)
                           if g.edges[eid] !== nothing && g.edges[eid].is_tree_edge]
                all(cv -> cv !== nothing && cv >= 0, tree_cv) || return false
                tree_slacks = [slack(g, eid, minimum_length) for eid in eachindex(g.edges)
                               if g.edges[eid] !== nothing && g.edges[eid].is_tree_edge]
                all(==(0), tree_slacks) || return false
                return minimum(v.rank for v in g.verts) == 1
            end

            g = sugi_from_edges(EXAMPLE_GRAPH, 8)
            rank!(g, 1, :networksimplex)
            @test is_correct(g, 1)

            # a bigger, denser DAG
            import Random
            rng = Random.MersenneTwister(42)
            n = 80
            edges = Tuple{Int,Int}[]
            for i in 0:(n - 2), _ in 1:2
                push!(edges, (i, rand(rng, (i + 1):min(i + 5, n - 1))))
            end
            g2 = sugi_from_edges(unique(edges), n)
            rank!(g2, 1, :networksimplex)
            @test is_correct(g2, 1)
        end

        @testset "verify_looks_good graph matches rust's expected width/height" begin
            # rust's `verify_looks_good` test asserts width==4.0, height==6.0
            # (max nodes on a rank, and number of ranks) for this graph.
            edges = [(0, 1), (1, 2), (2, 3), (2, 4), (3, 5), (3, 6), (3, 7), (3, 8),
                     (4, 5), (4, 6), (4, 7), (4, 8), (5, 9), (6, 9), (7, 9), (8, 9)]
            g = sugi_from_edges(edges, 10)
            remove_cycles!(g)
            rank!(g, 1, :networksimplex)
            @test length(unique(v.rank for v in g.verts)) == 6
            insert_dummy_vertices!(g, 1, 0.0)
            layers = ordering(g, :barycenter, true)
            @test maximum(length.(layers)) == 4
        end
    end

    # ---- Phase 2: crossing reduction, ported from p2_reduce_crossings/tests.rs
    @testset "ordering / crossing minimization (ported rust fixtures)" begin
        @testset "bilayer_cross_count" begin
            g = sugi_from_edges([(0, 4), (1, 3), (2, 3)], 5) # 0-based ids -> n0,n1,n2 / s0,s1
            for v in 1:3
                g.verts[v].rank = 0
            end
            for v in 4:5
                g.verts[v].rank = 1
            end
            order = SugiOrder([[1, 2, 3], [4, 5]], 5)
            @test bilayer_cross_count(g, order, 1) == 2

            g2 = sugi_from_edges([(0, 7), (1, 6), (2, 5), (3, 4)], 8)
            for v in 1:4
                g2.verts[v].rank = 0
            end
            for v in 5:8
                g2.verts[v].rank = 1
            end
            order2 = SugiOrder([[1, 2, 3, 4], [5, 6, 7, 8]], 8)
            @test bilayer_cross_count(g2, order2, 1) == 6
        end

        @testset "total crossings" begin
            edges = [(0, 6), (1, 7), (1, 8), (2, 6), (2, 9), (2, 10), (3, 6), (3, 9),
                     (4, 9), (5, 8), (5, 10)]
            g = sugi_from_edges(edges, 11)
            for v in 1:6
                g.verts[v].rank = 0
            end
            for v in 7:11
                g.verts[v].rank = 1
            end
            order = SugiOrder([collect(1:6), collect(7:11)], 11)
            @test total_crossings(g, order) == 12
        end

        @testset "insert_dummy_vertices" begin
            edges = [(0, 1), (1, 2), (2, 3), (3, 7), (4, 6), (5, 6), (6, 7), (0, 4), (0, 5)]
            g = sugi_from_edges(edges, 8)
            for (v, r) in [(0, 0), (1, 1), (2, 2), (3, 3), (4, 1), (5, 1), (6, 2), (7, 4)]
                g.verts[v + 1].rank = r
            end
            insert_dummy_vertices!(g, 1, 0.0)
            @test _nv(g) == 9
            @test count(v -> v.is_dummy, g.verts) == 1
        end

        @testset "barycenter reordering" begin
            g = sugi_from_edges([(0, 8), (1, 8), (2, 8), (2, 10), (3, 9), (4, 10), (5, 10),
                                  (6, 11), (7, 12)], 13)
            for v in 1:8
                g.verts[v].rank = 0
            end
            for v in 9:13
                g.verts[v].rank = 1
            end
            inner = [[1, 3, 5, 4, 7, 8, 2, 6], collect(9:13)]
            order = SugiOrder(inner, 13)
            result = order_layer(g, false, order, barycenter)
            @test result.layers[1] == collect(1:8)
        end
    end

    # ---- Phase 3: Brandes & Köpf alignment, ported from
    # ---- p3_calculate_coordinates/tests.rs's `create_test_layout` fixture
    @testset "coordinate assignment (ported rust fixtures)" begin
        edges30 = [(0, 2), (0, 6), (0, 18), (1, 16), (1, 17), (3, 8), (16, 8), (4, 8),
                   (17, 19), (18, 20), (5, 8), (5, 9), (6, 8), (6, 21), (7, 10), (7, 11),
                   (7, 12), (19, 23), (20, 24), (21, 12), (9, 22), (9, 25), (10, 13),
                   (10, 14), (11, 14), (22, 13), (23, 15), (24, 15), (12, 15), (25, 15)]
        layers0 = [[0, 1], [2, 3, 16, 4, 17, 18, 5, 6], [7, 8, 19, 20, 21, 9],
                   [10, 11, 22, 23, 24, 12, 25], [13, 14, 15]]

        function fixture()
            g = sugi_from_edges(edges30, 26)
            layers = [[v + 1 for v in row] for row in layers0]
            for (rank, row) in enumerate(layers), (pos, v) in enumerate(row)
                vert = g.verts[v]
                vert.rank, vert.pos = rank, pos
                vert.root = vert.align = vert.sink = v
                vert.is_dummy = (v - 1) >= 16
                vert.width = vert.height = (v - 1) >= 16 ? 1.0 : 10.0
            end
            return g, layers
        end
        root_of(g, id) = g.verts[id + 1].root - 1

        @testset "type-1 conflicts" begin
            g, l = fixture()
            mark_type1_conflicts!(g, l)
            for (t, h) in [(6, 8), (7, 12), (5, 8), (9, 22)]
                @test g.edges[find_edge(g, t + 1, h + 1)].has_type1_conflict
            end
        end

        @testset "type-1 conflicts: last vertex of a row must be checked (rust-port #27)" begin
            # Two ranks: upper = [A,B,C,D,Dm] (Dm a dummy continuing a chain into
            # the lower rank), lower = [E,Dm2,F] (Dm2 continues that chain, F is
            # an ordinary vertex and also the last position in the lower rank).
            # A->F should be flagged: drawing it straight would cross the Dm-Dm2
            # inner segment.
            g = SugiGraph()
            for _ in 1:8
                _add_vertex!(g)
            end
            A, B, C, D, Dm, E, Dm2, F = 1:8
            _add_edge!(g, B, E)   # ordinary edge, well inside any window
            _add_edge!(g, Dm, Dm2) # inner segment (dummy -> dummy)
            _add_edge!(g, A, F)   # crosses the inner segment: must be flagged

            layers = [[A, B, C, D, Dm], [E, Dm2, F]]
            g.verts[Dm].is_dummy = true
            g.verts[Dm2].is_dummy = true
            reset_alignment!(g, layers)
            mark_type1_conflicts!(g, layers)

            @test !g.edges[find_edge(g, B, E)].has_type1_conflict
            @test g.edges[find_edge(g, A, F)].has_type1_conflict
        end

        @testset "down-right alignment (exact root & align)" begin
            g, l = fixture()
            mark_type1_conflicts!(g, l)
            reset_alignment!(g, l)
            create_vertical_alignments!(g, l)

            exp_root = Dict(0=>0,1=>1,2=>0,3=>3,4=>4,5=>5,6=>6,7=>7,8=>4,9=>9,10=>7,
                             11=>11,12=>6,13=>7,14=>11,15=>18,16=>1,17=>17,18=>18,19=>17,
                             20=>18,21=>6,22=>22,23=>17,24=>18,25=>9)
            for (id, r) in exp_root
                @test root_of(g, id) == r
            end

            exp_align = Dict(0=>2,1=>16,2=>0,3=>3,4=>8,5=>5,6=>21,7=>10,8=>4,9=>25,10=>13,
                              11=>14,12=>6,13=>7,14=>11,15=>18,16=>1,17=>19,18=>20,19=>23,
                              20=>24,21=>12,22=>22,23=>17,24=>15,25=>9)
            for (id, a) in exp_align
                @test g.verts[id + 1].align - 1 == a
            end
        end

        @testset "down-left / up-right / up-left alignment (block membership)" begin
            g, l = fixture()
            mark_type1_conflicts!(g, l)
            foreach(reverse!, l)
            reset_alignment!(g, l)
            create_vertical_alignments!(g, l)
            for (root, members) in Dict(0=>[0,6], 4=>[4,8], 17=>[17,19,23], 18=>[18,20,24],
                                         5=>[5,9,25], 7=>[7,11,14], 21=>[21,12,15], 10=>[10,13])
                for m in members
                    @test root_of(g, m) == root
                end
            end

            g2, l2 = fixture()
            mark_type1_conflicts!(g2, l2)
            reverse_graph!(g2)
            reverse!(l2)
            reset_alignment!(g2, l2)
            create_vertical_alignments!(g2, l2)
            for (root, members) in Dict(13=>[13,10], 14=>[14,11,7], 15=>[15,23,19,17],
                                         24=>[24,20,18,0], 12=>[12,21], 25=>[25,9,5], 8=>[8,3])
                for m in members
                    @test root_of(g2, m) == root
                end
            end

            g3, l3 = fixture()
            mark_type1_conflicts!(g3, l3)
            reverse_graph!(g3)
            reverse!(l3)
            foreach(reverse!, l3)
            reset_alignment!(g3, l3)
            create_vertical_alignments!(g3, l3)
            for (root, members) in Dict(15=>[15,25,9], 13=>[13,22], 12=>[12,21,6],
                                         24=>[24,20,18], 23=>[23,19,17,1], 11=>[11,7], 8=>[8,4])
                for m in members
                    @test root_of(g3, m) == root
                end
            end
        end

        @testset "place_blocks sinks" begin
            g, l = fixture()
            mark_type1_conflicts!(g, l)
            create_vertical_alignments!(g, l)
            x = place_blocks(g, l)
            @test length(x) == 26
            for v in [0,1,2,3,4,5,6,8,9,12,15,16,17,18,19,20,21,23,24,25]
                @test g.verts[v + 1].sink - 1 == 0
            end
            for v in [7,10,11,22,13,14]
                @test g.verts[v + 1].sink - 1 == 7
            end
        end
    end

    # ---- End-to-end API-level tests -----------------------------------
    @testset "end-to-end layout()" begin
        @testset "empty and trivial graphs" begin
            @test Sugiyama()(zeros(Int, 0, 0)) == Point{2,Float64}[]
            @test length(Sugiyama()(zeros(Int, 1, 1))) == 1
        end

        @testset "self loops are ignored" begin
            adj = [1 1; 0 0]
            pos = Sugiyama()(adj)
            @test length(pos) == 2
            @test all(isfinite, pos[1]) && all(isfinite, pos[2])
        end

        @testset "cyclic graphs are handled (implicit edge reversal)" begin
            adj = zeros(Int, 5, 5)
            for (i, j) in [(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 3)]
                adj[i, j] = 1
            end
            pos = Sugiyama()(adj)
            @test length(pos) == 5
            @test all(p -> all(isfinite, p), pos)
        end

        @testset "disconnected components are tiled without overlap" begin
            adj = [0 1 0 0 0;
                   0 0 0 0 0;
                   0 0 0 1 0;
                   0 0 0 0 1;
                   0 0 0 0 0]
            pos = Sugiyama()(adj)
            @test length(pos) == 5
            @test length(unique(pos)) == 5
        end

        @testset "wheel_graph via Graphs.jl" begin
            g = wheel_graph(10)
            dirg = SimpleDiGraph(collect(edges(g)))
            adj = adjacency_matrix(dirg)
            pos = @time Sugiyama()(adj)
            @test typeof(pos) == Vector{Point{2,Float64}}
            @test pos == sugiyama(adj)
            @test pos == Sugiyama()(dirg)
        end

        @testset "Ptype and direction keywords" begin
            adj = adjacency_matrix(SimpleDiGraph(path_digraph(5)))
            pos = Sugiyama(; Ptype=Float32)(adj)
            @test typeof(pos) == Vector{Point{2,Float32}}

            for d in (:down, :up, :left, :right)
                pos = Sugiyama(; direction=d)(adj)
                @test length(pos) == 5
                @test all(p -> all(isfinite, p), pos)
            end
            # :down and :up should be vertical mirror images (same |y| magnitudes)
            pd = Sugiyama(; direction=:down)(adj)
            pu = Sugiyama(; direction=:up)(adj)
            @test getindex.(pd, 2) == -getindex.(pu, 2)
        end

        @testset "ranking_type / crossing_minimization / transpose combinations" begin
            adj = adjacency_matrix(SimpleDiGraph(path_digraph(6)))
            for rt in (:networksimplex, :longestpath, :up, :down),
                cm in (:barycenter, :median), tr in (true, false)

                pos = Sugiyama(; ranking_type=rt, crossing_minimization=cm, transpose=tr)(adj)
                @test length(pos) == 6
            end
        end

        @testset "nodesize / nodespacing / dummysize keywords" begin
            adj = adjacency_matrix(SimpleDiGraph(path_digraph(4)))
            pos = Sugiyama(; nodesize=[1.0, 2.0, 0.5], nodespacing=2.0, dummysize=0.3)(adj)
            @test length(pos) == 4
        end
    end

    @testset "assert square" begin
        M1 = rand(2, 4)
        @test_throws ArgumentError sugiyama(M1)
    end
end
