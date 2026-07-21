module ArcheoGra
    
    using DocStringExtensions, Plots, Random, GeoDatasets
    import Dates, Random.seed!, Plots.plot

    export  ### base 
            distance,
            set_seed,
    	    
            ### Examples
            example_networks_given_data_fullLBK_magda_new,
            example_networks_from_own_data,
            example_spectral_plots_given_data,
            example_spectral_plots_from_own_data,
            example_cluster_results_given_data,
            example_cluster_results_from_own_data,
            example_site_to_site_lcpa_uniform_given_file,
            example_site_to_site_lcpa_uniform_own_file,
            example_site_to_site_lcpa_adaptive_given_file,
            example_site_to_site_lcpa_adaptive_own_file,
            example_benchmark_adaptive_lcpa,
            
            ###########################################################################
            ### ConNetw
            ## types
            AbstractDataSet,
            AbstractDataEntry,
            ArchDataEntry,
            SourceDataEntry,
            Source_DataSet,
            Arch_DataSet,
            MaterialEntry,
            MaterialDataList,
            file_props,
            PreparedTimePeriod,

            ## constants
            earth_radius,
            maximum_distance,
            material_multiplier,

            ## common functions
            distance,
            setup_global_vars,
            prepare_network_data,
            get_sites_for_time,
            write_results,
            _output_path,

            ## read_data
            num_entries, 
            add_data!, 
            get_coordinates, 
            check_presence, 
            get_material_id, 
            add_material!,
            read_data,
            read_data_TOML,
            read_source_data,
            plot_sites,
            plot_material,

            ## Haversine_networks
            distance_allocation,
            Haversine_is_connected,
            Haversine_is_connected_sources,
            Haversine_is_connected_source2source,
            Haversine_connect_nodes,
            Haversine_connect_nodes_sources,
            Haversine_plot_connections,
            Haversine_plot_connections_Sources,

            ## connect_nodes_funcs
            adjacency_data,
            RayleighSISO,
            connect_nodes,
            adjacency_sources,

            ## testfunctions
            # ArcheoGra_graph_comparison,
            # Output_function_1104,
            # test_LCPAdistance,
            # test_LCPAdistance_nonsymm,
            # test_crosstime,
            
            ## connect_components_adj
            components, 
            connect_components!,
            custom_argmin,

            ## Plot_Prob
            # plot_prob,
            # total_prob,

            ## compare_graphs
            matrix_spectral_distance,
            effective_resistance,
            feature_extraction,
            canberra_distance,
            feature_aggregation,
            SRGG,
            get_matrices,
            produce_results,

            ### Clusterings
            ## Utility and Basic Operations
            extract_clusters,
            mcl_inflate_prune!,
            column_normalise!,
            build_coarse_adjacency,

            ## Spectral Clustering and K-Means
            identify_spectral_gap,
            spectral_kmeans,  
            FluidC_clustering,
            FluidC_clustering_modified,
            MCL,
            R_MCL,
            restarted_MCL,
            restarted_R_MCL,
            normalized_cut,

            ## Cluster Matching
            jaccard_distance_matrix,
            hungarian_match,
            match_clusters,
            
            ## Quality Measures
            modularity,
            silhouette_index,
            conductance,
            coverage,
            performance,

            ## Orchestrator and Plotting
            ClusterMethod,
            archaeological_clusterings,
            plot_clusters,

            ## Graph Data and GNN
            GraphData,
            GNNLayer,
            apply_layer,
            strength_aggregation_adj,
            strength_of_connection_gnn,

            ## DRA and Quality Control
            DRA_QC_adj,
            DRA_QC_CE_adj,
            DRA_QC_CE_detailed_adj,

            ##  Internal Functions
            _ncut_bisect,
            _dra_tentative_adj,
            _compute_gamma,
            _qc_extract,
            _bad_vertices_removal,
            _ext_weight,
            _int_weight,
            _criterion_16,
            _cholesky_quality_test,
            _build_AG_XG,
            _subgroup_extraction,           
            _dra_qc_ce_internal_adj,
            build_coarse_graph,
            iterative_deep_clustering,

            ## Spectral Plots
            smoothed_eigenvalue_dist,
            Lorentz_kernel,
            spectral_plots_comparison,
            archaeological_spectral_plots,

            ## archeoLCPA
            # connect_nodes_LCPA_symm,
            # archeoLCPA,
            # all_pairs_LCPA_GMT,
            # connect_nodes_LCPA_non_symm,
            # connect_nodes_LCPA,

            ## cross_timeperiods
            # create_connections,
            # one_to_rest_LCPA,

            ## full_material_comparison
            compare_all_is_connected, 
            compare_all_connect_nodes,
            compare_all_is_connected_sources,
            compare_all_connect_nodes_sources,

            ## FeatureComparison
            extract_degree_distribution,
            extract_clustering_coefficient,
            extract_egonet_features,
            extract_path_features,
            extract_density,
            # _heappush!,
            # _heappop!,
            feature_aggregation,
            network_features,
            mahalanobis_pinv,
            mahalanobis_distance,
            standardized_distance_squared, 
            empirical_quantile,
            empirical_cdf_interpolated, 
            create_graph_statistic_Haversine,
            is_average,
            graph_sensitivity,
            plot_network_robustness,
            
            ######################################################################################
            ##LCPA
            # Cost Functions (Constants)
            LloberaSluckin,
            Tobler_Hiking_adjusted,

            # Adjacency Matrix Configurations
            WeightConfig,
            Kondo_Average,
            Haversine_Config,
            
            # Types
            AbstractLCPAData,
            TopographicalData,
            ConnectionFunction,

            # AdjMatrix_Funcs
            weighted_adjacency_matrix,
            weighted_adjacency_matrix_haversine,
            weighted_adjacency_matrix_combi,

            # Cost_Funcs
            Tobler_Hiking,
            KS_Hiking,
            adjusted_Hiking,
            f_LloberaSluckin,
            distance,

            # LCPA_Main
            neighbours,
            delConnections,
            compress,
            get_path,
            all_pairs_LCPA,
            all_pairs_LCPA_GMT,
            all_pairs_LCPA_adaptive,
            all_pairs_LCPA_adaptive_GMT,
            extract_quadratic_map,
            time_comparison,
            # example_lcpa_main,
            dijkstra_shortest_path,

            # Indices_Numbering_etc
            characteristic,
            coords2lexi,
            coords2leaf_index,
            ind2coords,
            lexiMat2Vec,
            get_xy,
            lexi2ind,

            # InOut_Funcs
            tiff_to_matrix,
            write_Output,

            # WaterMap
            WaterLandBool,
            adaptive_WaterLandBool,

            # Wavelets
            coeffi_mat_Mℓ_0,
            coeffi_mat_Mℓ_1,
            N_m_and_Wavelet,
            Mℓ_and_Gℓ,
            N_1,
            N_2,
            N_2_boundary,
            create_Tℓ_inv,
            create_Tℓ,
            compression_1D,
            compression_2D,
            largestNTerm,
            relErrorTerm,
            heatmap_custom,
            example_compression,

            # Adaptive Refinement
            build_G,
            adaptive_leaf_cells,
            G_to_int,
            get_A_j,
            get_adaptive_matrices,
            leaf_edge_masks,
            cut_leaves,
            plot_quadtree_leaves,
            leaf_cells_to_level_matrix,
            group_leaves_by_level,
            leaf_adjacency_matrix_directed,
            get_leaves_adjacency_matrix,
            get_weighted_adaptive_matrix,
            example_weighted_adaptive_matrix,
            run_adaptive_hat_pipeline,
            benchmark_adaptive_lcpa,
            example_adaptive,
            leaves_to_hat_nodes,
            hat_adjacency_matrix,
            get_weighted_hat_matrix
    
    """
        set_seed(seed::Int)

    Set the global random seed for the reproducibility of stochastic network generation.

    # Arguments
    - `seed::Int`: The integer seed value.

    # Output
    - `Nothing`
    """
    function set_seed(seed::Int)
        Random.seed!(seed)
    end

    ######################################################################
    ######################################################################
    """
        earth_radius::Float64

    The radius of the Earth in kilometers, primarily calibrated for Central Europe.

    This constant serves as the baseline for spatial distance calculations, ensuring correct adjustments for the Earth's curvature when computing direct connections (e.g., Haversine distance).
    """
    const earth_radius = 6371 ## Earth radius in central europe (https://de.wikipedia.org/wiki/Erdradius)

    ######################################################################
    ######################################################################

    include("common_function.jl")
    include("Examples.jl")

    include("LCPA/LCPA.jl")

    include("ConNetw/ConNetw.jl")
    
    using .ConNetw
    using .LCPA
end