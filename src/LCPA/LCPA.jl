module LCPA

    using DocStringExtensions  

    using XLSX,  DelimitedFiles,  
        SparseArrays, LinearAlgebra, Plots,
        ProgressMeter, 
        GeoDatasets, ColorSchemes, Images, BenchmarkTools, Printf
    import CSV, GMT, Dates

    import ArcheoGra.distance, ArcheoGra.components, ArcheoGra.earth_radius

    export  # Cost Functions (Constants)
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

    abstract type AbstractLCPAData end

    """
        TopographicalData
    
    A mutable structure representing a topographical dataset for least-cost path analysis.

    $(FIELDS)
    """
    mutable struct TopographicalData <: AbstractLCPAData
        "Cell Raster"
        CostMap :: AbstractArray
        "Size of ElevationMap (should be the size of finest res)"
        size :: Tuple
        "Left bottom"
        lBot :: Tuple
        "Right top"
        rTop :: Tuple
        "Resolution"
        kl :: Tuple
    end

    """
        ConnectionFunction

    A wrapper type for travel cost functions, indicating whether the function is directionally symmetric.

    $(FIELDS)
    """
    struct ConnectionFunction
        "Function"
        f :: Function
        "Boolean Parameter for Symmetry or not"
        symmetric :: Bool
    end


    """
        distance(u::Tuple, v::Tuple)

    Calculates the great-circle distance between two geographic coordinates using the Haversine formula.

    # Arguments
    - `u::Tuple`: A tuple `(latitude, longitude)` for the first point in degrees.
    - `v::Tuple`: A tuple `(latitude, longitude)` for the second point in degrees.

    # Output
    - `Float64`: The calculated geographic distance in kilometers.
    """
    function distance(u::Tuple, v::Tuple)
        d2g(x)=deg2rad(x)
        return 2*earth_radius*asin(sqrt(0.5*(1-cos(d2g(v[2])-d2g(u[2]))+
                cos(d2g(u[2]))*cos(d2g(v[2]))*(1-cos(d2g(v[1])-d2g(u[1]))))))
    end

    """
        N_m_and_Wavelet

    A structure defining the properties and filter masks of a specific wavelet configuration.

    $(FIELDS)
    """
    struct N_m_and_Wavelet
        "Translation indices of phi"
        Delta_l::Function
        "Translation indices of psi"
        Nabla_l::Function
        "Mask a_k"
        Mask
        "Mask b_k"
        Wavelet_Mask
        "Inverse for mask a_k"
        inv_Mask
        "Inverse for mask b_k"
        inv_Wavelet_Mask
    end

    """
    N_1

    Base-level wavelet configuration using simple dyadic scaling (e.g., Haar-like characteristics).

    # Details
    `N_m_and_Wavelet(x -> 2^x, x -> 2^x, [1, 1], [1, -1], [1, 1], [1, -1])`
    """
    const N_1 = N_m_and_Wavelet(
        x -> 2^x,
        x -> 2^x,
        [1, 1],
        [1, -1],
        [1, 1],
        [1, -1]
    )

    """
    N_2

    Second-level wavelet configuration with offset scaling for smoother representations.

    # Details
    `N_m_and_Wavelet(x -> 2^x - 1, x -> 2^x, [1/2, 1, 1/2], [1, 0, 0], [0, 1, 0], [-1/2, 1, -1/2])`
    """
    const N_2 = N_m_and_Wavelet(
        x -> 2^x-1,
        x -> 2^x,
        [1/2, 1, 1/2],
        [1, 0, 0],
        [0, 1, 0],
        [-1/2, 1, -1/2]
    )

    """
    N_2_boundary

    Wavelet configuration with boundary corrections for extended finite domains.

    # Details
    Includes extended masks for handling edges without severe boundary artifacts.
    """
    const N_2_boundary = N_m_and_Wavelet(
        x -> 2^x+1,
        x -> 2^x,
        [1/2, 1, 1/2],
        [-1/4, -1/2, 3/2, -1/2, -1/4],
        [-1/4, 1/2, 3/2, 1/2, -1/4],
        [-1/2, 1, -1/2]
    )
    
    #const earth_radius = 6371 ## Earth radius in central europe (https://de.wikipedia.org/wiki/Erdradius)

    include("Adaptive_Refinement.jl")
    include("AdjMatrix_Funcs.jl")
    include("Cost_Funcs.jl")
    include("Indices_Numbering_etc.jl")
    include("InOut_Funcs.jl")
    include("LCPA_Main.jl")
    include("WaterMap.jl")
    include("Wavelets.jl")
    include("Adaptive_Refinement_Hat.jl")
end