using ArcheoGra
using Test
using SparseArrays
using LinearAlgebra

@testset "LCPA.jl - Integration Tests" begin

    # =========================================================================
    # 1. COST FUNCTIONS & UTILITIES
    # =========================================================================
    @testset "1. Cost Functions & Math" begin
        @test Tobler_Hiking(0.0) > 0.0
        @test KS_Hiking(0.0) > 0.0
        @test adjusted_Hiking(0.0) > 0.0
        @test f_LloberaSluckin(0.0) == 1.0
        
        # Distance beteen Equator and 1° North should be approximately 111 km
        @test 110.0 < distance((0.0, 0.0), (1.0, 0.0)) < 112.0
        
        # Test for Indicator Function
        @test characteristic(5, 0, 10, 10) == 1
    end

    # =========================================================================
    # 2. CLASSIC GRAPH PIPELINE
    # Tests: ind2coords, get_xy, Haversine_Config, weighted_adjacency_matrix_combi,
    #         dijkstra_shortest_path (MinHeap), get_path
    # =========================================================================
    @testset "2. Classic Grid & Pathfinding Pipeline" begin
        # Create a synthetic 16x16 map
        cost_map = rand(16, 16)
        top_data = TopographicalData(cost_map, (16, 16), (0.0, 0.0), (1.0, 1.0), (1, 1))
        
        # Pipeline 1: Graph-Creation
        # Triggers thousands of internal calls to ind2coords, distance, and the WeightConfig
        adj_combi = weighted_adjacency_matrix_combi(top_data, Haversine_Config, Cost_func=Tobler_Hiking)
        
        @test size(adj_combi) == (256, 256)
        @test SparseArrays.nnz(adj_combi) > 0
        
        # Pipeline 2: Dijkstra
        # Triggers the new MinHeap and the Early-Stopping logic
        costs, pred = dijkstra_shortest_path(adj_combi, 1, targets=[256])
        
        @test length(costs) == 256
        @test costs[1] == 0.0          # Start node has cost 0
        @test costs[256] != -1.0       # Target was reached
        
        # Pipeline 3: Path Extraction
        path = get_path(256, 1, pred)
        @test length(path) >= 2        # At least start and target
        @test path[1] == 256           # Path starts (backwards) at the target
        @test path[end] == 1           # Path ends at the start
    end

    # =========================================================================
    # Tests: compression_2D, all Mℓ/Gℓ Wavelet matrices, N_2 Wavelet,
    #         build_G, adaptive_leaf_cells, leaves_to_hat_nodes,
    #         hat_adjacency_matrix, get_weighted_hat_matrix
    #         hat_adjacency_matrix, get_weighted_hat_matrix
    # =========================================================================
    @testset "3. Wavelets & Adaptive Hat-Graph Pipeline" begin
        # Create a 32x32 map (Wavelets often require 2^L dimensions)
        L_geom = 5 # 2^5 = 32
        cost_map = rand(31, 31)
        top_data = TopographicalData(cost_map, (31, 31), (0.0, 0.0), (1.0, 1.0), (1, 1))
        
        # Step A: Wavelet Compression & Matrizen
        A, comp_rate = get_adaptive_matrices(cost_map, N_2, r=0.1)
        @test comp_rate >= 0.0
        @test haskey(A, L_geom)
        
        # Step B: Quadtree Construction
        G = build_G(A)
        leaves = adaptive_leaf_cells(G, L_geom)
        @test length(leaves) > 0
        @test length(leaves) <= 32*32 # Cannot have more leaves than pixels
        
        # Step C: Hat Nodes
        hat_nodes = leaves_to_hat_nodes(leaves, L_geom)
        @test length(hat_nodes) == length(leaves)
        
        # Step D: Adjacency and Weighting
        Ad_hat01 = hat_adjacency_matrix(hat_nodes, L_geom)
        Ad_weighted = get_weighted_hat_matrix(Ad_hat01, hat_nodes, top_data, L_geom)

        @test size(Ad_weighted) == (length(hat_nodes), length(hat_nodes))
        @test nnz(Ad_weighted) > 0 # Graph must have edges
    end

    # =========================================================================
    # 4. IN-OUT FUNCTIONS
    # =========================================================================
    @testset "4. InOut Functions" begin
        dummy_vec = [1.5, 2.5, 3.5]
        test_file = "test_output_dummy.txt"
        
        write_Output(dummy_vec, test_file)
        @test isfile(test_file)
        
        # Aufräumen
        rm(test_file)
    end
end

@testset "ConNetw.jl - Integration Tests" begin

    # =========================================================================
    # 1. CONNECTION PROBABILITIES & MATH
    # =========================================================================
    @testset "1. Math & Rayleigh SISO" begin
        # Base property: distance 0 -> probability 1
        @test RayleighSISO(0.0, 50.0) == 1.0
        
        # Base property: exactly at normalization distance
        @test RayleighSISO(50.0, 50.0, beta=10.0, eta=2.0) ≈ exp(-10.0)
        
        # Base property: half distance
        @test RayleighSISO(25.0, 50.0, beta=10.0, eta=2.0) ≈ exp(-10.0 * 0.25)
    end

    # =========================================================================
    # 2. CORE DATA STRUCTURES & ALLOCATIONS
    # Tests: prepare_network_data, distance_allocation, adjacency_data
    # =========================================================================
    @testset "2. Data Structures & Graph Allocation" begin
        # 2.1 Setup Mock Data
        mat_list = MaterialDataList(Dict(
            "Flint_LBK" => MaterialEntry(1, ("LBK", 100.0)),
            "Obsidian_LBK" => MaterialEntry(2, ("LBK", 300.0))
        ))

        site1 = ArchDataEntry(id=1, name="SiteA", coords=(6.0, 50.0), time_period="LBK", 
                              material_presence=["Flint_LBK"], num_material_presence=[1], 
                              material_multiplicity=[10], grid_id=1)
        
        site2 = ArchDataEntry(id=2, name="SiteB", coords=(7.0, 51.0), time_period="LBK", 
                              material_presence=["Flint_LBK", "Obsidian_LBK"], num_material_presence=[1, 2], 
                              material_multiplicity=[5, 2], grid_id=2)
        
        dataset = Arch_DataSet([site1, site2])
        
        # 2.2 Test Pre-Indexing (PreparedTimePeriod)
        prep_data = prepare_network_data(dataset, mat_list)
        @test haskey(prep_data, "LBK")
        
        lbk_prep = prep_data["LBK"]
        @test length(lbk_prep.sites) == 2
        @test haskey(lbk_prep.material_map, 1) # Material ID 1 (Flint) is present
        
        # 2.3 Test Geographic Distance Allocation
        dist_mat = distance_allocation(dataset.Data)
        @test size(dist_mat) == (2, 2)
        @test dist_mat[1, 1] == 0.0
        @test dist_mat[1, 2] == dist_mat[2, 1]
        @test dist_mat[1, 2] > 100.0 # Distance between roughly (6,50) and (7,51) is > 100km
        
        # 2.4 Test Material Adjacency Extraction
        adj_dict = adjacency_data(mat_list, lbk_prep)
        @test haskey(adj_dict, "Flint_LBK")
        
        flint_adj, flint_max_dist = adj_dict["Flint_LBK"]
        @test flint_max_dist == 100.0
        @test flint_adj[1, 2] == 15.0 # site1(10) + site2(5)
        @test flint_adj[2, 1] == 15.0 # Symmetry check
    end

    # =========================================================================
    # 3. NETWORK FEATURES & GRAPH STATISTICS
    # Tests: extract_density, extract_degree_distribution, clustering_coefficient
    # =========================================================================
    @testset "3. Graph Feature Extraction" begin
        # Create a synthetic fully connected 3-node graph (Triangle)
        Ad = sparse([
            0.0 1.0 1.0;
            1.0 0.0 1.0;
            1.0 1.0 0.0
        ])
        
        # Density (Edges / Possible Edges)
        @test extract_density(Ad) == 1.0 
        
        # Degree Distributions
        deg, two_hop = extract_degree_distribution(Ad)
        @test deg == [2, 2, 2]
        @test two_hop == [2.0, 2.0, 2.0]
        
        # Clustering Coefficient (Triangles / Possible Triangles)
        cc, two_hop_cc = extract_clustering_coefficient(Ad)
        @test cc == [1.0, 1.0, 1.0]
        @test two_hop_cc == [1.0, 1.0, 1.0]
    end

    # =========================================================================
    # 4. CLUSTERING QUALITY & UTILITIES
    # Tests: coverage, modularity, column_normalise!, mcl_inflate_prune!
    # =========================================================================
    @testset "4. Clustering Quality & Matrix Tools" begin
        # 4.1 Test Markov Matrix Utilities
        M = sparse([1.0 2.0; 3.0 4.0])
        
        M_norm = column_normalise!(copy(M))
        @test sum(M_norm[:, 1]) ≈ 1.0
        @test sum(M_norm[:, 2]) ≈ 1.0
        @test M_norm[1, 1] ≈ 0.25 # 1.0 / (1.0 + 3.0)
        
        M_inf = mcl_inflate_prune!(copy(M_norm), 2, prune_ratio=0.0)
        @test sum(M_inf[:, 1]) ≈ 1.0
        @test M_inf[2, 1] > M_norm[2, 1] # Squaring inflates larger probabilities
        
        # 4.2 Test Clustering Quality Metrics
        Ad = sparse([
            0.0 1.0 0.0;
            1.0 0.0 0.0;
            0.0 0.0 0.0
        ])
        clusters = [[1, 2], [3]]
        
        # Coverage: internal edges / total edges
        cov = coverage(Ad, clusters)
        @test cov == 1.0 # 100% of edges are internal to the clusters

        # Modularity: fraction of internal edges minus expected fraction
        mod_score = modularity(Ad, clusters)
        @test mod_score == 0.0 # Expected internal edges exactly match actual for this tiny layout
    end

    # =========================================================================
    # 5. COMPONENT UNIFICATION
    # Tests: components, connect_components!
    # =========================================================================
    @testset "5. Component Connectivity" begin
        # Create a disconnected graph: Node 1-2 are connected, Node 3 is isolated
        Ad = sparse([
            0.0 1.0 0.0;
            1.0 0.0 0.0;
            0.0 0.0 0.0
        ])
        
        dist_mat = [
            0.0 10.0 50.0;
            10.0 0.0 100.0;
            50.0 100.0 0.0
        ]
        
        # Step A: Identify disjoint components
        labels, comp_list = components(Ad)
        @test length(comp_list) == 2
        @test labels[1] == labels[2]
        @test labels[1] != labels[3]
        
        # Step B: Connect components via shortest distance
        connect_components!(Ad, labels, comp_list, dist_mat)
        
        # Node 3 should now be connected to Node 1 (distance 50 is shorter than 100)
        @test Ad[1, 3] == 1.0
        @test Ad[3, 1] == 1.0
        
        # Re-check components to ensure full connectedness
        _, new_comp_list = components(Ad)
        @test length(new_comp_list) == 1
    end

end