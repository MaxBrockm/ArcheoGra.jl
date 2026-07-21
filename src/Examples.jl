# function ArchaeoNetworks_Haversine_All_Sources()

#     # create_graph_statistic_Haversine(".\\src\\DataFolder\\Setup_Magdalenien.toml", num_graphs=500, double_distance = true)
#     # create_graph_statistic_Haversine(".\\src\\DataFolder\\Setup_FullLBK.toml", num_graphs=500)
#     # create_graph_statistic_Haversine(".\\src\\DataFolder\\Setup_oldLBK.toml", num_graphs=500)
#     # create_graph_statistic_Haversine(".\\src\\DataFolder\\Setup_Magdalenien_new.toml", num_graphs=500, double_distance = true)

#     # Ad_magda, _, _ = graph_sensitivity(".\\src\\DataFolder\\Setup_Magdalenien.toml", num_graphs=500, double_distance = true)
#     # _, Ad_oLBK, _ = graph_sensitivity(".\\src\\DataFolder\\Setup_oldLBK.toml", num_graphs=500)
#     # _, _, Ad_fullLBK = graph_sensitivity(".\\src\\DataFolder\\Setup_FullLBK.toml", num_graphs=500)
#     # Ad_magda_new, _, _ = graph_sensitivity(".\\src\\DataFolder\\Setup_Magdalenien_new.toml", num_graphs=500, double_distance = true)

#     try
#         mkpath(".\\src\\OutputFolder\\$(string(Dates.today()))")
#     catch 
#         @info "Possibly overwriting files"
#     end
#     # Random.seed!(350) ## Accepted
#     # ArData, MaData = read_data_TOML(".\\src\\DataFolder\\Setup_Magdalenien.toml")
#     # SoData = read_source_data(".\\src\\DataFolder\\Setup_Magdalenien.toml", MaData)
#     # graphen = Haversine_connect_nodes_sources(ArData, MaData, SoData, Plot_compare=true, plot_results=true, double_distance = true)
#     # cp("src\\DataFolder\\Setup_Magdalenien.toml", "src\\OutputFolder\\$(string(Dates.today()))\\Setup_Magdalenien.toml",force=true)
#     # open("src\\OutputFolder\\$(string(Dates.today()))\\Setup_Magdalenien.toml", "a") do io
#     #     write(io, "\nmaximum_distance = $maximum_distance\nmaterial_multiplier=$material_multiplier" )
#     # end

#     # for i in 1:2:length(graphen)
#     #     graph = graphen[i]
#     #     comparison_graph = graphen[i+1]
#     #     time_period = graph[1][1].time_period
#     #     Ad_weighted = graph[4]
        
#     #     @show sum(Ad_weighted .> 0) ÷ 2
#     #     comparison_Ad_weighted = comparison_graph[4]
#     #     max_dist = maximum(distance_allocation(graph[1]))
#     #     @show max_dist
#     #     println("Features of graph from time period $time_period with materials:")
#     #     network_features(Ad_weighted, print_results=true, max_dist)
#     #     println("\n\nFeatures of graph from time period $time_period withOUT materials:")
#     #     network_features(comparison_Ad_weighted, print_results=true, max_dist)
#     #     avg, _, _ , _  = is_average(Ad_weighted, max_dist, "Haversine", time_period, incl_sources=true)

#     # end
#     # Random.seed!(81) ## Accepted
#     # ArData, MaData = read_data_TOML(".\\src\\DataFolder\\Setup_oldLBK.toml")
#     # SoData = read_source_data(".\\src\\DataFolder\\Setup_oldLBK.toml", MaData)
#     # graphen = Haversine_connect_nodes_sources(ArData, MaData, SoData, Plot_compare=true, plot_results=true, double_distance = false)
#     # cp("src\\DataFolder\\Setup_oldLBK.toml", "src\\OutputFolder\\$(string(Dates.today()))\\Setup_oldLBK.toml",force=true)
#     # open("src\\OutputFolder\\$(string(Dates.today()))\\Setup_oldLBK.toml", "a") do io
#     #     write(io, "\nmaximum_distance = $maximum_distance\nmaterial_multiplier=$material_multiplier" )
#     # end
#     # for i in 1:2:length(graphen)
#     #     graph = graphen[i]
#     #     comparison_graph = graphen[i+1]
#     #     time_period = graph[1][1].time_period
#     #     Ad_weighted = graph[4]
#     #     @show length(components(Ad_weighted)[2])
#     #     comparison_Ad_weighted = comparison_graph[4]
#     #     max_dist = maximum(distance_allocation(graph[1])) ## Flawed but okay

#     #     println("Features of graph from time period $time_period with materials:")
#     #     network_features(Ad_weighted, print_results=true, max_dist)
#     #     println("\n\nFeatures of graph from time period $time_period withOUT materials:")
#     #     network_features(comparison_Ad_weighted, print_results=true, max_dist)
#     #     avg, _, _ , _  = is_average(Ad_weighted, max_dist, "Haversine", time_period, incl_sources=true)
#     # end
    
#     Random.seed!(58) ## Accepted
#     ArData, MaData = read_data_TOML(".\\src\\DataFolder\\Setup_FullLBK.toml")
#     SoData = read_source_data(".\\src\\DataFolder\\Setup_FullLBK.toml", MaData)
#     graphen = Haversine_connect_nodes_sources(ArData, MaData, SoData, Plot_compare=true, plot_results=true, double_distance = false)
#     cp("src\\DataFolder\\Setup_FullLBK.toml", "src\\OutputFolder\\$(string(Dates.today()))\\Setup_FullLBK.toml",force=true)
#     open("src\\OutputFolder\\$(string(Dates.today()))\\Setup_FullLBK.toml", "a") do io
#         write(io, "\nmaximum_distance = $maximum_distance\nmaterial_multiplier=$material_multiplier" )
#     end
#     for i in 1:2:length(graphen)
#         graph = graphen[i]
#         comparison_graph = graphen[i+1]
#         time_period = graph[1][1].time_period
#         Ad_weighted = graph[4]
#         comparison_Ad_weighted = comparison_graph[4]
#         max_dist = maximum(distance_allocation(graph[1]))
#         println("Features of graph from time period $time_period with materials:")
#         network_features(Ad_weighted, print_results=true, max_dist)
#         println("\n\nFeatures of graph from time period $time_period withOUT materials:")
#         network_features(comparison_Ad_weighted, print_results=true, max_dist)
#         avg, _, _ , _  = is_average(Ad_weighted, max_dist, "Haversine", time_period, incl_sources=true)
#     end

#     Random.seed!(13) ## Accepted
#     ArData, MaData = read_data_TOML(".\\src\\DataFolder\\Setup_Magdalenien_new.toml")
#     SoData = read_source_data(".\\src\\DataFolder\\Setup_Magdalenien_new.toml", MaData)
#     graphen = Haversine_connect_nodes_sources(ArData, MaData, SoData, Plot_compare=true, plot_results=true, double_distance = true)
#     cp("src\\DataFolder\\Setup_Magdalenien_new.toml", "src\\OutputFolder\\$(string(Dates.today()))\\Setup_Magdalenien_new.toml",force=true)
#     open("src\\OutputFolder\\$(string(Dates.today()))\\Setup_Magdalenien_new.toml", "a") do io
#         write(io, "\nmaximum_distance = $maximum_distance\nmaterial_multiplier=$material_multiplier" )
#     end

#     for i in 1:2:length(graphen)
#         graph = graphen[i]
#         comparison_graph = graphen[i+1]
#         time_period = graph[1][1].time_period
#         Ad_weighted = graph[4]
#         comparison_Ad_weighted = comparison_graph[4]
#         max_dist = maximum(distance_allocation(graph[1]))
#         println("Features of graph from time period $time_period with materials:")
#         network_features(Ad_weighted, print_results=true, max_dist)
#         println("\n\nFeatures of graph from time period $time_period withOUT materials:")
#         network_features(comparison_Ad_weighted, print_results=true, max_dist)
#         avg, _, _ , _  = is_average(Ad_weighted, max_dist, "Haversine", time_period, incl_sources=true)
#     end

#     return
# end


# ===============================
# Published Example Wrappers
# ===============================

"""
    _maybe_analyze_haversine_graphs(graphen; incl_sources::Bool)

Internal helper for optional post-processing of generated Haversine graphs.
If enabled by the caller, this prints feature summaries and runs `is_average`.
"""
function _maybe_analyze_haversine_graphs(graphen; incl_sources::Bool)
    if length(graphen) < 2
        return nothing
    end

    for i in 1:2:(length(graphen)-1)
        graph = graphen[i]
        comparison_graph = graphen[i+1]

        if length(graph) < 4 || length(comparison_graph) < 4 || isempty(graph[1])
            continue
        end

        time_period = graph[1][1].time_period
        Ad_weighted = graph[4]
        comparison_Ad_weighted = comparison_graph[4]
        max_dist = maximum(distance_allocation(graph[1]))

        println("Features of graph from time period $time_period with materials:")
        network_features(Ad_weighted, max_dist; print_results=true)

        println("\n\nFeatures of graph from time period $time_period withOUT materials:")
        network_features(comparison_Ad_weighted, max_dist; print_results=true)

        is_average(Ad_weighted, max_dist, "Haversine", time_period; incl_sources=incl_sources)
    end

    return nothing
end

"""
    example_networks_given_data_fullLBK_magda_new(; output_dir::String=pwd(), overwrite::Bool=false, plot_compare::Bool=true, plot_results::Bool=false, analyze_features::Bool=false)

Creates Haversine-based archaeological networks for the package-provided datasets
`Setup_FullLBK.toml` and `Setup_Magdalenien_new.toml`.

# Output
- `Dict{String, Any}` with keys `"fullLBK"` and `"Magdalenien_new"`.

# Optional Analysis
- If `analyze_features=true`, this function also prints `network_features` summaries and runs `is_average`.
"""
function example_networks_given_data_fullLBK_magda_new(; output_dir::String=pwd(), overwrite::Bool=false, plot_compare::Bool=true, plot_results::Bool=true, analyze_features::Bool=false)
    out_dir = joinpath(output_dir, string(Dates.today()), "ExampleNetworks")
    if isdir(out_dir)
        if !overwrite
            error("Output directory '$out_dir' already exists. Use overwrite=true to continue.")
        end
    else
        mkpath(out_dir)
    end

    full_cfg = joinpath(@__DIR__, "DataFolder", "Setup_FullLBK.toml")
    magda_cfg = joinpath(@__DIR__, "DataFolder", "Setup_Magdalenien_new.toml")

    full_data, full_materials = read_data_TOML(full_cfg)
    full_sources = read_source_data(full_cfg, full_materials)
    full_graphs = Haversine_connect_nodes_sources(
        full_data,
        full_materials,
        full_sources,
        Plot_compare=plot_compare,
        plot_results=plot_results,
        double_distance=false,
        output_dir=out_dir,
        overwrite=overwrite,
    )

    magda_data, magda_materials = read_data_TOML(magda_cfg)
    magda_sources = read_source_data(magda_cfg, magda_materials)
    magda_graphs = Haversine_connect_nodes_sources(
        magda_data,
        magda_materials,
        magda_sources,
        Plot_compare=plot_compare,
        plot_results=plot_results,
        double_distance=true,
        output_dir=out_dir,
        overwrite=true,
    )

    if analyze_features
        _maybe_analyze_haversine_graphs(full_graphs, incl_sources=true)
        _maybe_analyze_haversine_graphs(magda_graphs, incl_sources=true)
    end

    cp(full_cfg, joinpath(out_dir, basename(full_cfg)), force=true)
    cp(magda_cfg, joinpath(out_dir, basename(magda_cfg)), force=true)

    return Dict(
        "fullLBK" => full_graphs,
        "Magdalenien_new" => magda_graphs,
    )
end


"""
    example_networks_from_own_data(config_file::String; incl_sources::Bool=true, double_distance::Bool=false, output_dir::String=pwd(), overwrite::Bool=false, plot_compare::Bool=true, plot_results::Bool=true, analyze_features::Bool=false)

Template example for custom data.

# Setup
1. Create a TOML config analogous to the files in `src/DataFolder`.
2. Ensure file paths in that TOML point to your archaeological/source tables.
3. Call this function with your TOML path.

# Output
- Graph container returned by `Haversine_connect_nodes` or `Haversine_connect_nodes_sources`.

# Optional Analysis
- If `analyze_features=true`, this function also prints `network_features` summaries and runs `is_average`.
"""
function example_networks_from_own_data(config_file::String; incl_sources::Bool=true, double_distance::Bool=false, output_dir::String=pwd(), overwrite::Bool=false, plot_compare::Bool=true, plot_results::Bool=true, analyze_features::Bool=false)
    cfg = isfile(config_file) ? config_file : joinpath(@__DIR__, "DataFolder", config_file)
    if !isfile(cfg)
        error("Could not find TOML config '$config_file'.")
    end

    out_dir = joinpath(output_dir, string(Dates.today()), "ExampleNetworksOwn")
    if isdir(out_dir)
        if !overwrite
            error("Output directory '$out_dir' already exists. Use overwrite=true to continue.")
        end
    else
        mkpath(out_dir)
    end

    ArData, MaData = read_data_TOML(cfg)
    graphs = if incl_sources
        SoData = read_source_data(cfg, MaData)
        Haversine_connect_nodes_sources(
            ArData,
            MaData,
            SoData,
            Plot_compare=plot_compare,
            plot_results=plot_results,
            double_distance=double_distance,
            output_dir=out_dir,
            overwrite=true,
        )
    else
        Haversine_connect_nodes(
            ArData,
            MaData,
            Plot_compare=plot_compare,
            plot_results=plot_results,
            output_dir=out_dir,
            overwrite=true,
        )
    end

    if analyze_features
        _maybe_analyze_haversine_graphs(graphs, incl_sources=incl_sources)
    end

    cp(cfg, joinpath(out_dir, basename(cfg)), force=true)
    return graphs
end


"""
    example_spectral_plots_given_data(; num_graphs::Int=100, seed::Union{Nothing,Int}=nothing, output_dir::String=pwd(), overwrite::Bool=false)

Creates archaeological spectral plots for package datasets `Setup_FullLBK.toml`
and `Setup_Magdalenien_new.toml`.
"""
function example_spectral_plots_given_data(; num_graphs::Int=100, seed::Union{Nothing,Int}=nothing, output_dir::String=pwd(), overwrite::Bool=false)
    out_dir = joinpath(output_dir, string(Dates.today()), "ExampleSpectralPlots")
    if isdir(out_dir)
        if !overwrite
            error("Output directory '$out_dir' already exists. Use overwrite=true to continue.")
        end
    else
        mkpath(out_dir)
    end

    full_cfg = joinpath(@__DIR__,  "DataFolder", "Setup_FullLBK.toml")
    magda_cfg = joinpath(@__DIR__, "DataFolder", "Setup_Magdalenien_new.toml")

    archaeological_spectral_plots(full_cfg, num_graphs=num_graphs, seed=seed, incl_sources=true, double_distance=false, output_dir=out_dir, overwrite=true)
    archaeological_spectral_plots(magda_cfg, num_graphs=num_graphs, seed=seed, incl_sources=true, double_distance=true, output_dir=out_dir, overwrite=true)

    return out_dir
end


"""
    example_spectral_plots_from_own_data(config_file::String; num_graphs::Int=100, seed::Union{Nothing,Int}=nothing, incl_sources::Bool=true, double_distance::Bool=false, output_dir::String=pwd(), overwrite::Bool=false)

Runs `archaeological_spectral_plots` for your own TOML configuration.

# Setup
- Provide a valid config TOML with your site/material (and optionally source) data.
"""
function example_spectral_plots_from_own_data(config_file::String; num_graphs::Int=100, seed::Union{Nothing,Int}=nothing, incl_sources::Bool=true, double_distance::Bool=false, output_dir::String=pwd(), overwrite::Bool=false)
    cfg = isfile(config_file) ? config_file : joinpath(@__DIR__, "DataFolder", config_file)
    if !isfile(cfg)
        error("Could not find TOML config '$config_file'.")
    end

    out_dir = joinpath(output_dir, string(Dates.today()), "ExampleSpectralPlotsOwn")
    if isdir(out_dir)
        if !overwrite
            error("Output directory '$out_dir' already exists. Use overwrite=true to continue.")
        end
    else
        mkpath(out_dir)
    end

    archaeological_spectral_plots(
        cfg,
        num_graphs=num_graphs,
        seed=seed,
        incl_sources=incl_sources,
        double_distance=double_distance,
        output_dir=out_dir,
        overwrite=true,
    )

    return out_dir
end


"""
    example_cluster_results_given_data(; output_dir::String=pwd(), overwrite::Bool=false, num_clusters::Int=0, max_clusters_param::Int=10, plot_map::Bool=false, seed::Any=nothing)

Creates clustering outputs via `archaeological_clusterings` for
`Setup_FullLBK.toml` and `Setup_Magdalenien_new.toml`.

# Output
- `Dict{String, String}`: Output paths for both dataset runs.
"""
function example_cluster_results_given_data(; output_dir::String=pwd(), overwrite::Bool=false, num_clusters::Int=0, max_clusters_param::Int=10, plot_map::Bool=false, seed::Any=nothing)
    out_dir = joinpath(output_dir, string(Dates.today()), "ExampleClusterResults")
    if isdir(out_dir)
        if !overwrite
            error("Output directory '$out_dir' already exists. Use overwrite=true to continue.")
        end
    else
        mkpath(out_dir)
    end

    full_cfg = "Setup_FullLBK.toml"
    magda_cfg = "Setup_Magdalenien_new.toml"

    full_base_out = joinpath(out_dir, "fullLBK")
    magda_base_out = joinpath(out_dir, "Magdalenien_new")
    @show full_base_out
    archaeological_clusterings(
        full_cfg;
        num_clusters=num_clusters,
        max_clusters_param=max_clusters_param,
        plot_map=plot_map,
        seed=seed,
        incl_sources=true,
        double_distance=false,
        output_dir=full_base_out,
        overwrite=true,
    )

    archaeological_clusterings(
        magda_cfg;
        num_clusters=num_clusters,
        max_clusters_param=max_clusters_param,
        plot_map=plot_map,
        seed=seed,
        incl_sources=true,
        double_distance=true,
        output_dir=magda_base_out,
        overwrite=true,
    )

    return Dict(
        "fullLBK" => joinpath(full_base_out, string(Dates.today()), "Clusterings"),
        "Magdalenien_new" => joinpath(magda_base_out, string(Dates.today()), "Clusterings"),
    )
end


"""
    example_cluster_results_from_own_data(config_file::String; incl_sources::Bool=true, double_distance::Bool=false, output_dir::String=pwd(), overwrite::Bool=false, num_clusters::Int=0, max_clusters_param::Int=10, plot_map::Bool=false, seed::Any=nothing)

Runs `archaeological_clusterings` on your own TOML data.

# Setup
1. Prepare your TOML file and linked data files.
2. Run this function with the TOML path.
3. Find generated files inside the returned output directory.
"""
function example_cluster_results_from_own_data(config_file::String; incl_sources::Bool=true, double_distance::Bool=false, output_dir::String=pwd(), overwrite::Bool=false, num_clusters::Int=0, max_clusters_param::Int=10, plot_map::Bool=false, seed::Any=nothing)
    cfg = isfile(config_file) ? config_file : joinpath(@__DIR__, "DataFolder", config_file)
    if !isfile(cfg)
        error("Could not find TOML config '$config_file'.")
    end

    out_dir = joinpath(output_dir, string(Dates.today()), "ExampleClusterResultsOwn")
    if isdir(out_dir)
        if !overwrite
            error("Output directory '$out_dir' already exists. Use overwrite=true to continue.")
        end
    else
        mkpath(out_dir)
    end

    archaeological_clusterings(
        cfg;
        num_clusters=num_clusters,
        max_clusters_param=max_clusters_param,
        plot_map=plot_map,
        seed=seed,
        incl_sources=incl_sources,
        double_distance=double_distance,
        output_dir=out_dir,
        overwrite=true,
    )

    return joinpath(out_dir, string(Dates.today()), "Clusterings")
end


"""
    _resolve_wavelet_choice(wavelet_choice::Symbol)

Maps the explicit user choice `:N_1` or `:N_2` to the corresponding wavelet object.
"""
function _resolve_wavelet_choice(wavelet_choice::Symbol)
    if wavelet_choice == :N_1
        return N_1
    elseif wavelet_choice == :N_2
        return N_2
    else
        error("Unsupported wavelet_choice=$wavelet_choice. Choose :N_1 or :N_2.")
    end
end


"""
    _site_to_site_lcpa_pipeline(tiff_file::String; start_site, target_site, margin_ratio::Float64, L_ref::Int, wavelet_type=N_1, r_tol::Float64=-1.0, N_tol::Int=-1, adaptive::Bool=false, plot_result::Bool=true, output_dir::String=pwd(), overwrite::Bool=false)

Internal helper for site-to-site LCPA examples.
"""
function _site_to_site_lcpa_pipeline(
        tiff_file::String; 
        start_site::Tuple{Number, Number}, 
        target_site::Tuple{Number, Number}, 
        margin_ratio::Float64=0.5, 
        L_ref::Int=10, 
        wavelet_type=N_1, 
        r_tol::Float64=-1.0,
        N_tol::Int=-1,
        adaptive::Bool=false, 
        plot_result::Bool=true, 
        output_dir::String=pwd(), 
        overwrite::Bool=false
        )
    TopoData_quad, Ad_weighted, hat_nodes, leaves, L_geom, comp_rate = if adaptive
        if (r_tol > 0) == (N_tol > 0)
            error("For adaptive mode, choose exactly one threshold: set either r_tol > 0 (residual) or N_tol > 0 (Best-N-term).")
        end

        run_adaptive_hat_pipeline(
            tiff_file,
            wavelet_type,
            r=r_tol,
            N=N_tol,
            start_site=start_site,
            target_site=target_site,
            margin_ratio=margin_ratio,
        )
    else
        LCPA.run_uniform_reference_pipeline(
            tiff_file,
            L_ref,
            wavelet_type,
            start_site=start_site,
            target_site=target_site,
            margin_ratio=margin_ratio,
            plot_cutout=false,
        )
    end
    mask = GeoDatasets.LandSeaMask(;resolution='f', grid=1.25)
    if !any(pt -> GeoDatasets.is_land(mask, pt[2], pt[1]), [start_site, target_site])
        error("Either the start site or the target site is not on land.")
    end
    start_idx = LCPA.find_nearest_hat_node(start_site, hat_nodes, TopoData_quad)
    target_idx = LCPA.find_nearest_hat_node(target_site, hat_nodes, TopoData_quad)
    
    keep_nodes = adaptive_WaterLandBool(leaves, TopoData_quad, L_geom, [start_idx, target_idx])
    Ad_clean = delConnections(Ad_weighted, keep_nodes)

    costs, pred = dijkstra_shortest_path(Ad_clean, start_idx)
    path_indices = get_path(target_idx, start_idx, pred)
    if isempty(path_indices)
        return (path_coords=Tuple{Number,Number}[], total_cost=Inf, comp_rate=comp_rate, output_file=nothing)
    end

    reverse!(path_indices)
    path_coords = Tuple{Number,Number}[]
    for idx in path_indices
        node = hat_nodes[idx]
        r = clamp(round(Int, node.row), 1, size(TopoData_quad.CostMap, 1))
        c = clamp(round(Int, node.col), 1, size(TopoData_quad.CostMap, 2))
        push!(path_coords, ind2coords((r, c), TopoData_quad))
    end

    output_file = nothing
    if plot_result
        out_dir = joinpath(output_dir, string(Dates.today()), adaptive ? "AdaptiveLCPA" : "UniformLCPA")
        if isdir(out_dir)
            if !overwrite
                error("Output directory '$out_dir' already exists. Use overwrite=true to continue.")
            end
        else
            mkpath(out_dir)
        end
        filename = adaptive ? "adaptive_site_to_site" : "uniform_site_to_site"
        LCPA.plot_lcp_gmt(tiff_file, path_coords, start_site, target_site, filename, output_dir=out_dir, overwrite=true)
        output_file = joinpath(out_dir, "$(filename).png")
    end

    return (path_coords=path_coords, total_cost=costs[target_idx], comp_rate=comp_rate, output_file=output_file)
end


"""
    example_site_to_site_lcpa_uniform_given_file(; tiff_file::String=joinpath("src", "DataFolder", "-13_35_46_70.tiff"), start_site::Tuple{Number, Number}=(-5.4, 38.6), target_site::Tuple{Number, Number}=(10.0, 45.0), L_ref::Int=10, margin_ratio::Float64=0.5, output_dir::String=pwd(), overwrite::Bool=false)

Uniform-grid site-to-site LCPA example using package-provided TIFF.
"""
function example_site_to_site_lcpa_uniform_given_file(; tiff_file::String=joinpath(dirname(@__DIR__), "src", "DataFolder", "-13_35_46_70.tiff"), start_site::Tuple{Number, Number}=(-5.4, 38.6), target_site::Tuple{Number, Number}=(10.0, 45.0), L_ref::Int=10, margin_ratio::Float64=0.5, output_dir::String=pwd(), overwrite::Bool=false)
    return _site_to_site_lcpa_pipeline(
        tiff_file,
        start_site=start_site,
        target_site=target_site,
        margin_ratio=margin_ratio,
        L_ref=L_ref,
        adaptive=false,
        plot_result=true,
        output_dir=output_dir,
        overwrite=overwrite,
    )
end


"""
    example_site_to_site_lcpa_uniform_own_file(tiff_file::String; start_site, target_site, L_ref::Int=10, margin_ratio::Float64=0.5, output_dir::String=pwd(), overwrite::Bool=false)

Uniform-grid site-to-site LCPA example for your own TIFF file.

# Setup
- Provide a georeferenced TIFF with the same coordinate convention as ArcheoGra LCPA routines.
"""
function example_site_to_site_lcpa_uniform_own_file(tiff_file::String; start_site::Tuple{Number,Number}, target_site::Tuple{Number,Number}, L_ref::Int=10, margin_ratio::Float64=0.5, output_dir::String=pwd(), overwrite::Bool=false)
    return _site_to_site_lcpa_pipeline(
        tiff_file,
        start_site=start_site,
        target_site=target_site,
        margin_ratio=margin_ratio,
        L_ref=L_ref,
        adaptive=false,
        plot_result=true,
        output_dir=output_dir,
        overwrite=overwrite,
    )
end


"""
    example_site_to_site_lcpa_adaptive_given_file(; tiff_file::String=joinpath("src", "DataFolder", "-13_35_46_70.tiff"), start_site::Tuple{Number,Number}=(-5.4, 38.6), target_site::Tuple{Number,Number}=(10.0, 45.0), wavelet_choice::Symbol=:N_1, r_tol::Float64=-1.0, N_tol::Int=5000, margin_ratio::Float64=0.5, output_dir::String=pwd(), overwrite::Bool=false)

Adaptive site-to-site LCPA example using package-provided TIFF.

# Important
- Set `wavelet_choice=:N_1` or `wavelet_choice=:N_2` explicitly.
- Thresholding: choose exactly one mode.
: Residual mode with `r_tol > 0` and `N_tol = -1`.
: Best-N-term mode with `N_tol > 0` and `r_tol = -1.0`.
"""
function example_site_to_site_lcpa_adaptive_given_file(; tiff_file::String=joinpath("src", "DataFolder", "-13_35_46_70.tiff"), start_site::Tuple{Number,Number}=(-5.4, 38.6), target_site::Tuple{Number,Number}=(10.0, 45.0), wavelet_choice::Symbol=:N_1, r_tol::Float64=-1.0, N_tol::Int=5000, margin_ratio::Float64=0.5, output_dir::String=pwd(), overwrite::Bool=false)
    wavelet_type = _resolve_wavelet_choice(wavelet_choice)
    return _site_to_site_lcpa_pipeline(
        tiff_file,
        start_site=start_site,
        target_site=target_site,
        margin_ratio=margin_ratio,
        wavelet_type=wavelet_type,
        r_tol=r_tol,
        N_tol=N_tol,
        adaptive=true,
        plot_result=true,
        output_dir=output_dir,
        overwrite=overwrite,
    )
end


"""
    example_site_to_site_lcpa_adaptive_own_file(tiff_file::String; start_site, target_site, wavelet_choice::Symbol=:N_1, r_tol::Float64=-1.0, N_tol::Int=5000, margin_ratio::Float64=0.5, output_dir::String=pwd(), overwrite::Bool=false)

Adaptive site-to-site LCPA example for your own TIFF file.

# Setup
- Provide a georeferenced TIFF with the same coordinate convention as ArcheoGra LCPA routines.
- Set `wavelet_choice=:N_1` or `wavelet_choice=:N_2` explicitly.
- Thresholding: choose exactly one mode.
: Residual mode with `r_tol > 0` and `N_tol = -1`.
: Best-N-term mode with `N_tol > 0` and `r_tol = -1.0`.
"""
function example_site_to_site_lcpa_adaptive_own_file(tiff_file::String; start_site::Tuple{Number,Number}, target_site::Tuple{Number,Number}, wavelet_choice::Symbol=:N_1, r_tol::Float64=-1.0, N_tol::Int=5000, margin_ratio::Float64=0.5, output_dir::String=pwd(), overwrite::Bool=false)
    wavelet_type = _resolve_wavelet_choice(wavelet_choice)
    return _site_to_site_lcpa_pipeline(
        tiff_file,
        start_site=start_site,
        target_site=target_site,
        margin_ratio=margin_ratio,
        wavelet_type=wavelet_type,
        r_tol=r_tol,
        N_tol=N_tol,
        adaptive=true,
        plot_result=true,
        output_dir=output_dir,
        overwrite=overwrite,
    )
end


"""
    example_benchmark_adaptive_lcpa(filepath::String; start_site::Tuple{Number,Number}=(-5.4, 38.6), target_site::Tuple{Number,Number}=(10.0, 45.0), margin_ratio::Float64=0.5, L_ref::Int=10, plot_cutout::Bool=true, cutout_plot_backend::Symbol=:auto, wavelet_choices::Vector{Symbol}=[:N_1, :N_2], r_values::Vector{Float64}=[0.01, 0.025, 0.05, 0.1, 0.2], N_values::Vector{Int}=[5000, 20000, 50000, 75000, 150000], output_dir::String=pwd(), overwrite::Bool=false)

Convenience wrapper for `benchmark_adaptive_lcpa` with improved usage documentation.

# Setup
1. Choose a TIFF file (`filepath`) and start/target coordinates.
2. Choose wavelets with `wavelet_choices=[:N_1]`, `[:N_2]`, or `[:N_1, :N_2]`.
3. Set either `r_values`, `N_values`, or both for the benchmark sweep.
4. Run this function.

# Results
- Returns the benchmark tuple vector from `benchmark_adaptive_lcpa`.
- Also writes a CSV summary to:
  `joinpath(output_dir, string(Dates.today()), "AdaptiveBenchmark", "benchmark_adaptive_lcpa.csv")`
"""
function example_benchmark_adaptive_lcpa(filepath::String; start_site::Tuple{Number,Number}=(-5.4, 38.6), target_site::Tuple{Number,Number}=(10.0, 45.0), margin_ratio::Float64=0.5, L_ref::Int=10, plot_cutout::Bool=true, cutout_plot_backend::Symbol=:auto, wavelet_choices::Vector{Symbol}=[:N_1, :N_2], r_values::Vector{Float64}=[0.01, 0.025, 0.05, 0.1, 0.2], N_values::Vector{Int}=[5000, 20000, 50000, 75000, 150000], output_dir::String=pwd(), overwrite::Bool=false)
    out_dir = joinpath(output_dir, string(Dates.today()), "AdaptiveBenchmark")
    if isdir(out_dir)
        if !overwrite
            error("Output directory '$out_dir' already exists. Use overwrite=true to continue.")
        end
    else
        mkpath(out_dir)
    end

    results = benchmark_adaptive_lcpa(
        filepath,
        start_site,
        target_site,
        margin_ratio=margin_ratio,
        L_ref=L_ref,
        plot_cutout=plot_cutout,
        cutout_plot_backend=cutout_plot_backend,
        wavelet_choices=wavelet_choices,
        r_values=r_values,
        N_values=N_values,
    )

    csv_file = joinpath(out_dir, "benchmark_adaptive_lcpa.csv")
    open(csv_file, "w") do io
        write(io, "Method,Threshold_Type,Threshold_Value,Comp_Rate,Vertices,Edges,Cost,Distance,Rel_Error\n")
        for row in results
            write(io, string(row[1], ",", row[2], ",", row[3], ",", row[4], ",", row[5], ",", row[6], ",", row[7], ",", row[8], ",", row[9], "\n"))
        end
    end

    return (results=results, csv_file=csv_file)
end