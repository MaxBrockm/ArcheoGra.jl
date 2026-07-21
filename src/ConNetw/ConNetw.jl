module ConNetw

    using DocStringExtensions  

    using XLSX, SparseArrays, DelimitedFiles, ProgressMeter, Arpack, LinearAlgebra, Plots, Statistics, PrettyTables, LaTeXStrings, Random, Base.Threads
    import GMT, Graphs.barabasi_albert, Graphs.watts_strogatz, Graphs.adjacency_matrix, Dates, StatsBase.sample, TOML, Colors.distinguishable_colors

    import ..ArcheoGra: distance, earth_radius
    using ..LCPA

    export  ## types
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
            plot_network_robustness
    ######################################################################
    ######################################################################

    abstract type AbstractDataSet end
    abstract type AbstractDataEntry end

    """
        ArchDataEntry

    Type that represents a data entry from an archaeological dataset. Stores essential information for each entry.

    $(FIELDS)
    """
    @kwdef mutable struct ArchDataEntry <: AbstractDataEntry
        "ID of stored entry. Created independently of possible IDs in existing datasets."
        id::Int
        "Name of the archaeological site."
        name::String = ""
        "Coordinates of the archaeological site (longitude, latitude)."
        coords::Tuple{Float64, Float64} 
        "Time period of the archaeological site."
        time_period::String
        "List of names of present materials."
        material_presence::Vector{String}
        "IDs of materials found (ID corresponds to the Material Data Set)."
        num_material_presence::Vector{Int}
        "Amount of findings of the respective material."
        material_multiplicity::Vector{Int}
        "Possible prior classification or group assigned to the entry."
        entry_classification::Any = nothing
        "Spatial grid identifier."
        grid_id::Int
    end

    """
        SourceDataEntry

    Type that represents a data entry from the raw-material source dataset. 

    $(FIELDS)
    """
    @kwdef mutable struct SourceDataEntry <: AbstractDataEntry
        "ID of the raw material source."
        id::Int 
        "Name of the material found at this source."
        material_name::String
        "Coordinates of the raw material source (longitude, latitude)."
        coords::Tuple{Float64, Float64} 
        "Spatial grid identifier."
        grid_id::Int
    end

    """
        Source_DataSet

    Container type for the source dataset.

    $(FIELDS)
    """
    struct Source_DataSet <: AbstractDataSet
        "Vector of source data entries of type `SourceDataEntry`."
        Data::Vector{SourceDataEntry}
        "Multiplier for weighting sources when calculating connections."
        source_mult::Number
    end

    """
        Arch_DataSet

    Container type for the archaeological dataset.

    $(FIELDS)
    """
    struct Arch_DataSet <: AbstractDataSet
        "Vector of archaeological data entries of type `ArchDataEntry`."
        Data::Vector{ArchDataEntry}
    end


    """
        MaterialEntry

    Type representing a specific material data entry.

    $(FIELDS)
    """
    struct MaterialEntry <: AbstractDataEntry
        "Internal ID for the material in the local list of materials."
        material_id::Int
        "Property of the material entry consisting of the time period and distribution distance."
        mat_prop::Tuple{String, Number}
    end

    """
        MaterialDataList

    Container type for material storage. 

    $(FIELDS)
    """
    mutable struct MaterialDataList 
        "Dictionary assigning a unique string identifier to a material entry."
        Material_list::Dict{String, MaterialEntry}
    end

    """
        file_props
        
    Type for storing properties of a dataset, including columns mapping for parsing data entries.

    $(FIELDS)
    """
    @kwdef mutable struct file_props
        "Path of the file. Can be a relative or absolute file path."
        filePath::String
        "File type described as a string (e.g., 'xlsx', 'csv', 'txt')."
        file_type::String
        "Excel only: Sheet name of the dataset."
        sheetname::String = ""
        "For delimited files: The delimiter character."
        delim::Char = ' '
        "Range of columns where material presence is stored."
        material_cols::UnitRange
        "String representation of the time period for the materials in the data."
        material_age::String
        "Transport distance of the material in the data."
        material_distance::Vector{Float64}
        "Column index containing the grid identifier."
        grid_id::Int
        "Column index containing the longitude."
        long_ind::Int
        "Column index containing the latitude."
        lat_ind::Int
        "Column index containing the name of the archaeological site (if available)."
        name_ind::Any = nothing
        "Column index containing prior group classification (if available)."
        class_ind::Any = nothing
    end

    """
        PreparedTimePeriod

    Pre-indexed view of archaeological data for a single time period.
    Built once by `prepare_network_data` and reused by all network-creation functions,
    replacing repeated linear scans over the full dataset with O(1) lookups.

    $(FIELDS)
    """
    struct PreparedTimePeriod
        "Name of the time period."
        time_period::String
        "Sites belonging to this time period."
        sites::Vector{ArchDataEntry}
        "Indices of `sites` into the original `Arch_DataSet.Data` vector."
        site_indices::Vector{Int}
        "Per-site dictionary mapping `material_id` to `multiplicity`."
        site_material_dicts::Vector{Dict{Int, Int}}
        "Per-material mapping `material_id` to a list of tuples `(local_site_index, multiplicity)`."
        material_map::Dict{Int, Vector{Tuple{Int, Int}}}
    end

    """
        prepare_network_data(data::Arch_DataSet, material_data::MaterialDataList)

    Processes the flat `Arch_DataSet` into a dictionary keyed by time period.
    Every network-creation function can obtain filtered sites, material-to-site mappings,
    and per-site material lookups in O(1) instead of rescanning the dataset.

    # Arguments
    - `data::Arch_DataSet`: The full archaeological dataset.
    - `material_data::MaterialDataList`: The material catalogue.

    # Output
    - `Dict{String, PreparedTimePeriod}`: A dictionary containing one indexed entry per unique time period.

    # Examples
    ```julia
    prepared_data = prepare_network_data(my_arch_data, my_material_data)
    ```
    """
    function prepare_network_data(data::Arch_DataSet, material_data::MaterialDataList)
        tp_index = Dict{String, Vector{Int}}()
        for (i, entry) in enumerate(data.Data)
            push!(get!(Vector{Int}, tp_index, entry.time_period), i)
        end

        result = Dict{String, PreparedTimePeriod}()

        for (tp, indices) in tp_index
            sites = data.Data[indices]

            site_mat_dicts = Vector{Dict{Int, Int}}(undef, length(sites))
            for (li, entry) in enumerate(sites)
                d = Dict{Int, Int}()
                for k in eachindex(entry.num_material_presence)
                    d[entry.num_material_presence[k]] = entry.material_multiplicity[k]
                end
                site_mat_dicts[li] = d
            end

            mat_map = Dict{Int, Vector{Tuple{Int, Int}}}()
            for (_, mat) in material_data.Material_list
                if mat.mat_prop[1] == tp
                    mid = mat.material_id
                    bucket = Tuple{Int, Int}[]
                    for (li, d) in enumerate(site_mat_dicts)
                        mult = get(d, mid, 0)
                        if mult != 0
                            push!(bucket, (li, mult))
                        end
                    end
                    mat_map[mid] = bucket
                end
            end

            result[tp] = PreparedTimePeriod(tp, sites, indices, site_mat_dicts, mat_map)
        end
        return result
    end

    """
        get_sites_for_time(data::Arch_DataSet, time_period::String)

    Filters sites by a specific time period via a linear scan. Prefer using `PreparedTimePeriod` for better performance.

    # Arguments
    - `data::Arch_DataSet`: The archaeological dataset.
    - `time_period::String`: The targeted time period.

    # Output
    - `Vector{ArchDataEntry}`: Filtered list of entries.
    """
    function get_sites_for_time(data::Arch_DataSet, time_period::String)
        return filter(e -> e.time_period == time_period, data.Data)
    end

    """
        distance(u::ArchDataEntry, v::ArchDataEntry)
        
    Calculates the Haversine distance between two archaeological sites.

    # Arguments
    - `u::ArchDataEntry`: The first archaeological site.
    - `v::ArchDataEntry`: The second archaeological site.

    # Output
    - `Float64`: The distance in kilometers.
    """
    @inline function distance(u::ArchDataEntry, v::ArchDataEntry)
        d2g(x) = deg2rad(x)
        return 2.0 * earth_radius * asin(sqrt(0.5 * (1.0 - cos(d2g(v.coords[2]) - d2g(u.coords[2])) +
            cos(d2g(u.coords[2])) * cos(d2g(v.coords[2])) *
            (1.0 - cos(d2g(v.coords[1]) - d2g(u.coords[1]))))))
    end

    """
        setup_global_vars(setup_file::String)

    Overrides default parameters for the maximum walking distance and the material weight multiplier.
    Values are read from the `global_vars` section of the provided TOML file.

    # Arguments
    - `setup_file::String`: Path to the TOML configuration file.

    # Output
    - `Nothing`
    """
    function setup_global_vars(setup_file::String)
        cfg = TOML.parsefile(setup_file)["global_vars"]
        global maximum_distance = Float64(cfg["maximum_distance"])
        global material_multiplier = Float64(cfg["material_multiplier"])

        @info "Global variables correctly updated" maximum_distance material_multiplier
        return nothing
    end


    """
        write_results(path::String, adjacency::AbstractArray, adjacency_compare::AbstractArray, distance_mat::AbstractArray, ids::AbstractVector)

    Writes multiple adjacency matrices, difference matrices, and their distance-weighted equivalents to CSV files.

    # Arguments
    - `path::String`: Target directory for the exported CSV files.
    - `adjacency::AbstractArray`: The computed base adjacency matrix.
    - `adjacency_compare::AbstractArray`: The reference adjacency matrix for comparison.
    - `distance_mat::AbstractArray`: The matrix defining distances between points.
    - `ids::AbstractVector`: Identifiers corresponding to matrix rows and columns.

    # Output
    - `Nothing`
    """
    function write_results(path::String, adjacency::AbstractArray, adjacency_compare::AbstractArray, distance_mat::AbstractArray, ids::AbstractVector)
        if !isdir(path)
            @info "Creating new output directory" path
            mkpath(path)
        end
        
        id_row = reshape(ids, 1, :)
        
        writedlm(joinpath(path, "adjacency.csv"), vcat(hcat("", id_row), hcat(ids, adjacency)), ',')
        writedlm(joinpath(path, "adjacency_compare.csv"), vcat(hcat("", id_row), hcat(ids, adjacency_compare)), ',')
        writedlm(joinpath(path, "adjacency_difference.csv"), vcat(hcat("", id_row), hcat(ids, adjacency .- adjacency_compare)), ',')
        writedlm(joinpath(path, "adjacency_weighted.csv"), vcat(hcat("", id_row), hcat(ids, adjacency .* distance_mat)), ',')
        writedlm(joinpath(path, "adjacency_compare_weighted.csv"), vcat(hcat("", id_row), hcat(ids, adjacency_compare .* distance_mat)), ',')
        writedlm(joinpath(path, "adjacency_difference_weighted.csv"), vcat(hcat("", id_row), hcat(ids, (adjacency .- adjacency_compare) .* distance_mat)), ',')
        
        return nothing
    end

    """
        _output_path(parts::String...; output_dir::String = pwd())

    Constructs an output path dynamically and creates the required directory tree if it does not exist.

    # Arguments
    - `parts::String...`: The path segments to join.

    # Keywords
    - `output_dir::String`: Base directory. Defaults to the current working directory.

    # Output
    - `String`: The full resolved output path.
    """
    function _output_path(parts::String...; output_dir::String = pwd())
        path = joinpath(output_dir, "Output", string(Dates.today()), parts...)
        dir_path = dirname(path)
        if !isempty(dir_path) && !isdir(dir_path)
            mkpath(dir_path)
        end
        return path
    end
######################################################################
######################################################################
    #const earth_radius = 6371.0 ## Earth radius in central europe (https://de.wikipedia.org/wiki/Erdradius)

    """
        maximum_distance::Float64

    The default maximum walking distance between archaeological sites.

    This global variable defines the upper spatial threshold for drawing edges between vertices in the network. Its value can be dynamically updated by reading in the corresponding configuration file via `setup_global_vars`.
    """   
    global maximum_distance = 50.0

    """
        material_multiplier::Int

    The weight of raw material connectivity during network generation.

    This global variable determines the impact of shared raw materials within the parameter η of the Rayleigh function. It acts as a multiplier that increases the probability of connection between vertices if they share the presence of raw materials.
    """
    global material_multiplier = 1
######################################################################
######################################################################

    include("compare_graphs.jl")
    include("Haversine_networks.jl")
    include("Plot_Prob.jl")
    include("connected_components_adj.jl")
    include("read_data.jl")
    include("testfunctions.jl")
    include("archaeo_LCPA.jl")
    include("connect_nodes_funcs.jl")
    include("cross_timeperiods.jl")
    include("full_material_comparison.jl")
    include("spectral_plots.jl")
    include("clusterings.jl")
    include("FeatureComparison.jl")
end