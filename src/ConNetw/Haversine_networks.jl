"""
    distance_allocation(data::AbstractArray{ArchDataEntry})

Calculates the Haversine distance (as the bird flies) between every archaeological site. 
The distance is stored in a symmetric matrix.

# Arguments
- `data::AbstractArray{ArchDataEntry}`: Array containing entries of archaeological sites.

# Output
- `Matrix{Float64}`: Matrix containing distances between all given sites.
"""
function distance_allocation(data::AbstractArray{ArchDataEntry})
    num_e = length(data)
    distance_matrix = zeros(Float64, num_e, num_e)
    for start_node in 1:num_e
        for end_node in (start_node+1):num_e 
            d = distance(data[start_node], data[end_node])
            distance_matrix[start_node, end_node] = d
            distance_matrix[end_node, start_node] = d
        end
    end
    return distance_matrix
end

"""
    distance_allocation(data::AbstractArray{ArchDataEntry}, source_data::Source_DataSet)

Calculates the Haversine distance between every raw material source and archaeological site.

# Arguments
- `data::AbstractArray{ArchDataEntry}`: Array containing entries of archaeological sites.
- `source_data::Source_DataSet`: Dataset containing raw material sources.

# Output
- `Matrix{Float64}`: Distance matrix (rows: sources, columns: sites).
"""
function distance_allocation(data::AbstractArray{ArchDataEntry}, source_data::Source_DataSet)
    num_e = length(data)
    num_s = length(source_data.Data)
    distance_matrix = zeros(Float64, num_s, num_e)
    for (start_node, source) in enumerate(source_data.Data)
        for (end_node,site) in enumerate(data) 
            d = distance([source.coords...], [site.coords...])
            distance_matrix[start_node, end_node] = d
        end
    end
    return distance_matrix
end


"""
    distance_allocation(source_data::Source_DataSet)

Calculates the Haversine distance between all given raw material sources.

# Arguments
- `source_data::Source_DataSet`: Dataset containing raw material sources.

# Output
- `Matrix{Float64}`: Distance matrix between all sources.
"""
function distance_allocation(source_data::Source_DataSet)
    num_s = length(source_data.Data)
    distance_matrix = zeros(Float64, num_s, num_s)
    for start_node in 1:num_s
        for end_node in (start_node+1):num_s
            d = distance([source_data.Data[start_node].coords...], [source_data.Data[end_node].coords...])
            if d == 0
                d = 0.1
            end
            distance_matrix[start_node, end_node] = d
            distance_matrix[end_node, start_node] = d
        end
    end
    return distance_matrix
end


"""
    Haversine_is_connected(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, material_weight_matrices::Dict; relevant_sites::Vector{ArchDataEntry})

Determines if two vertices should be connected based on spatial distance and shared material presence.
The probability function relies on a Rayleigh model, where distance acts as a baseline, and shared materials increase connectivity likelihood. Group classifications can further increase connection probability.

# Arguments
- `start_node::Int`: Index of the first vertex.
- `end_node::Int`: Index of the second vertex.
- `distance_matrix::AbstractArray{Float64}`: Matrix containing distances between all vertices.
- `material_weight_matrices::Dict`: Dictionary of material matrices and their corresponding distribution ranges.

# Keywords
- `relevant_sites::Vector{ArchDataEntry}`: List of all relevant archaeological sites to evaluate group classifications.

# Output
- `Tuple{Bool, Bool}`: Returns two booleans: `(material_connection, default_connection)`.
"""
function Haversine_is_connected(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, material_weight_matrices::Dict; relevant_sites::Vector{ArchDataEntry})
    r = rand()

    eff_max_distance = maximum_distance 
    
    if !isnothing(relevant_sites) && (relevant_sites[start_node].entry_classification == relevant_sites[end_node].entry_classification)
        eff_max_distance *= 2
    end

    Prob = RayleighSISO(distance_matrix[start_node,end_node], eff_max_distance; eta=2)
    Prob1 = copy(Prob)

    for (name, mat) in material_weight_matrices
        (mat[1][start_node,end_node] == 0) && (continue)
        Prob += (1-Prob)*RayleighSISO(distance_matrix[start_node,end_node], mat[2], eta=2+material_multiplier*mat[1][start_node, end_node])
    end
    return  r < Prob, r < Prob1
end


"""
    Haversine_is_connected_sources(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, distance_matrix_sources::AbstractArray{Float64}, material_weight_matrices::Dict, sourcedata::Source_DataSet; relevant_sites::Vector{ArchDataEntry})

Determines connectivity between two archaeological sites by factoring in intermediate paths via raw material sources.

# Arguments
- `start_node::Int`: Index of the first site.
- `end_node::Int`: Index of the second site.
- `distance_matrix::AbstractArray{Float64}`: Distance matrix between sites.
- `distance_matrix_sources::AbstractArray{Float64}`: Distance matrix between sources and sites.
- `material_weight_matrices::Dict`: Dictionary containing material weights.
- `sourcedata::Source_DataSet`: Dataset containing material source locations.

# Keywords
- `relevant_sites::Vector{ArchDataEntry}`: List of relevant archaeological sites.

# Output
- `Tuple{Bool, Bool}`: Returns `(material_connection, default_connection)`.
"""
function Haversine_is_connected_sources(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, distance_matrix_sources::AbstractArray{Float64}, material_weight_matrices::Dict, sourcedata::Source_DataSet; relevant_sites::Vector{ArchDataEntry})
    r = rand()

    eff_max_distance = maximum_distance 
    
    if !isnothing(relevant_sites) && (relevant_sites[start_node].entry_classification == relevant_sites[end_node].entry_classification)
        eff_max_distance *= 2
    end

    Prob = RayleighSISO(distance_matrix[start_node,end_node], eff_max_distance; eta=2)
    Prob1 = copy(Prob)

    for (i,s) in enumerate(sourcedata.Data)
        haskey(material_weight_matrices, s.material_name) || continue
        mat = material_weight_matrices[s.material_name]
        (mat[1][start_node,end_node] == 0) && (continue)
        Prob += (1-Prob)*RayleighSISO(distance_matrix_sources[i,start_node]+distance_matrix_sources[i,end_node], mat[2], eta=2+material_multiplier*mat[1][start_node, end_node])
    end
    return  r < Prob, r < Prob1
end


"""
    Haversine_is_connected_source2source(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, shared_sites::AbstractArray, source_mult::Number; eff_distance::Number = 500, beta::Number = 25, double_distance::Bool = false)

Evaluates the direct connection probability between two material sources.

# Arguments
- `start_node::Int`: First source index.
- `end_node::Int`: Second source index.
- `distance_matrix::AbstractArray{Float64}`: Distance matrix between sources.
- `shared_sites::AbstractArray`: Matrix representing shared connected sites between sources.
- `source_mult::Number`: Multiplier modifying connection probability based on shared properties.

# Keywords
- `eff_distance::Number`: Base effective distance for the Rayleigh function. Defaults to 500.
- `beta::Number`: Beta parameter for the Rayleigh function. Defaults to 25.
- `double_distance::Bool`: If true, bypasses the probabilistic function and connects purely based on shared site overlap. Defaults to false.

# Output
- `Bool`: Connection status.
"""
function Haversine_is_connected_source2source(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, shared_sites::AbstractArray, source_mult::Number; eff_distance::Number = 500, beta::Number=25, double_distance=false)
    r = rand()
    if !double_distance
        Prob = RayleighSISO(distance_matrix[start_node,end_node], eff_distance; eta=2 + source_mult*shared_sites[start_node, end_node], beta=beta)
        return  r < Prob
    else
        connected = shared_sites[start_node, end_node] > 0 ? true : false
        return connected
    end
end


"""
    Haversine_connect_nodes(data::Arch_DataSet, MaterialData::MaterialDataList; Plot_compare::Bool = true, plot_results::Bool = true, prepared::Union{Dict{String,PreparedTimePeriod},Nothing} = nothing, output_dir::String = pwd(), overwrite::Bool = false)

Generates networks for all time periods from archaeological sites and their materials.
If disconnected components emerge, shortest-path auxiliary edges are injected to ensure a single connected component.

# Arguments
- `data::Arch_DataSet`: Complete archaeological dataset.
- `MaterialData::MaterialDataList`: List of materials and properties.

# Keywords
- `Plot_compare::Bool`: If true, differences between baseline and material-driven networks are highlighted. Defaults to true.
- `plot_results::Bool`: Determines if the resulting geographic maps are rendered. Defaults to true.
- `prepared::Union{Dict{String,PreparedTimePeriod},Nothing}`: Pre-indexed data structures to optimize generation. Defaults to nothing.
- `output_dir::String`: Base directory to save CSV and map exports. Defaults to current working directory.
- `overwrite::Bool`: Enables file overwrite privileges. Throws an error if set to false and the directory already exists.

# Output
- `Vector{Any}`: A collection of graph states and matrices mapped per time period.
"""
function Haversine_connect_nodes(data::Arch_DataSet, MaterialData::MaterialDataList; Plot_compare::Bool = true, plot_results::Bool = true, prepared::Union{Dict{String,PreparedTimePeriod},Nothing} = nothing, output_dir::String = pwd(), overwrite::Bool = false)
    prep_data = isnothing(prepared) ? prepare_network_data(data, MaterialData) : prepared
    time_periods = collect(keys(prep_data))
    graphen = []
    
    for time_period in time_periods
        prep = prep_data[time_period]
        @info "Generating network for time period: $time_period"
        
        base_path = joinpath(output_dir, string(Dates.today()), "Haversine", time_period)
        if isdir(base_path) && !overwrite
            error("Output directory $base_path already exists. Set `overwrite=true` to proceed.")
        else
            mkpath(base_path)
        end
        
        material_connectivity_arr = adjacency_data(MaterialData, prep)
        relevant_sites = prep.sites

        distance_matrix = distance_allocation(relevant_sites)
        anz_rel = length(relevant_sites)
        adjacency_matrix_daten = spzeros(Float64, anz_rel, anz_rel)
        adjacency_matrix_daten_compare = spzeros(Float64, anz_rel, anz_rel)
        
        for i in 1:anz_rel
            for j in (i+1):anz_rel
                res_mat, res_def = Haversine_is_connected(i, j, distance_matrix, material_connectivity_arr; relevant_sites=relevant_sites)
                adjacency_matrix_daten[i, j] = res_mat ? 1.0 : 0.0
                adjacency_matrix_daten_compare[i, j] = res_def ? 1.0 : 0.0
            end
        end
        
        adjacency_matrix_daten += adjacency_matrix_daten'
        adjacency_matrix_daten_compare += adjacency_matrix_daten_compare'
        
        non_conn_adj = copy(adjacency_matrix_daten)
        non_conn_adj_comp = copy(adjacency_matrix_daten_compare)

        writedlm(joinpath(base_path, "non_conn_$time_period.csv"), adjacency_matrix_daten, ',')
        writedlm(joinpath(base_path, "no_mat_non_conn_$time_period.csv"), adjacency_matrix_daten_compare, ',')
        writedlm(joinpath(base_path, "material_diff_$time_period.csv"), adjacency_matrix_daten .- adjacency_matrix_daten_compare, ',')
        writedlm(joinpath(base_path, "material_diff_distance_$time_period.csv"), sparse((adjacency_matrix_daten .- adjacency_matrix_daten_compare) .* distance_matrix), ',')
        
        connect_components!(adjacency_matrix_daten, components(adjacency_matrix_daten)..., distance_matrix)
        connect_components!(adjacency_matrix_daten_compare, components(adjacency_matrix_daten_compare)..., distance_matrix)
        
        writedlm(joinpath(base_path, "$time_period.csv"), adjacency_matrix_daten, ',')
        writedlm(joinpath(base_path, "no_mat_$time_period.csv"), adjacency_matrix_daten_compare, ',')

        weighted_adjacency_matrix = adjacency_matrix_daten .* distance_matrix
        writedlm(joinpath(base_path, "weight_$time_period.csv"), weighted_adjacency_matrix, ',')
        
        weighted_adjacency_matrix_comp = adjacency_matrix_daten_compare .* distance_matrix
        writedlm(joinpath(base_path, "weight_no_mat_$time_period.csv"), weighted_adjacency_matrix_comp, ',')

        push!(graphen, [relevant_sites, adjacency_matrix_daten, non_conn_adj, weighted_adjacency_matrix])
        push!(graphen, [relevant_sites, adjacency_matrix_daten_compare, non_conn_adj_comp, weighted_adjacency_matrix_comp, distance_matrix])
    end
    
    if plot_results
        Haversine_plot_connections(graphen, compare=Plot_compare, output_dir=output_dir, overwrite=overwrite)
    end
    return graphen
end

"""
    Haversine_connect_nodes_sources(data::Arch_DataSet, MaterialData::MaterialDataList, SourceDataDict::Dict; Plot_compare::Bool = true, plot_results::Bool = true, double_distance::Bool = true, prepared::Union{Dict{String,PreparedTimePeriod},Nothing} = nothing, output_dir::String = pwd(), overwrite::Bool = false)

Generates spatial networks inclusive of raw material sources mapping to archaeological sites.

# Arguments
- `data::Arch_DataSet`: Complete archaeological dataset.
- `MaterialData::MaterialDataList`: List of materials and properties.
- `SourceDataDict::Dict`: Dictionary of source records split by time periods.

# Keywords
- `Plot_compare::Bool`: Indicates if difference maps are generated. Defaults to true.
- `plot_results::Bool`: Evaluates geographic map outputs. Defaults to true.
- `double_distance::Bool`: Skips probabilistic filtering for internal source connections. Defaults to true.
- `prepared::Union{Dict{String,PreparedTimePeriod},Nothing}`: Pre-indexed memory views to bypass slow filtering. Defaults to nothing.
- `output_dir::String`: Root path for CSV arrays and maps. Defaults to current working directory.
- `overwrite::Bool`: Grants delete and replace permissions. Defaults to false.

# Output
- `Vector{Any}`: A collection of graph states mapping mixed-mode matrices per epoch.
"""
function Haversine_connect_nodes_sources(data::Arch_DataSet, MaterialData::MaterialDataList, SourceDataDict::Dict; Plot_compare::Bool = true, plot_results::Bool = true, double_distance::Bool = true, prepared::Union{Dict{String,PreparedTimePeriod},Nothing} = nothing, output_dir::String = pwd(), overwrite::Bool = false)
    prep_data = isnothing(prepared) ? prepare_network_data(data, MaterialData) : prepared
    time_periods = collect(keys(prep_data))
    graphen = []

    for time_period in time_periods
        prep = prep_data[time_period]
        @info "Generating source-inclusive network for time period: $time_period"
        
        SourceData = SourceDataDict[time_period][1] 
        source_base = joinpath(output_dir, string(Dates.today()), "Haversine", "Sources", time_period)
        
        if isdir(source_base) && !overwrite
            error("Output directory $source_base already exists. Set `overwrite=true` to proceed.")
        else
            mkpath(source_base)
        end
       
        source_adjacency, source_maxdistance = adjacency_sources(prep, SourceData, MaterialData)
        material_connectivity_arr = adjacency_data(MaterialData, prep)
        relevant_sites = prep.sites
        relevant_ids = [site.grid_id for site in relevant_sites] 
        source_ids = [source.grid_id for source in SourceData.Data]
        ids = [relevant_ids; source_ids]
     
        distance_site_site = distance_allocation(relevant_sites)
        distance_source_site = distance_allocation(relevant_sites, SourceData)
        distance_source_source = distance_allocation(SourceData)
        
        anz_rel = length(relevant_sites)
        anz_sources = length(SourceData.Data)

        s_adj_I, s_adj_J, _ = findnz(source_adjacency)
        adjacency_source2site = sparse(s_adj_I, s_adj_J, ones(length(s_adj_J)), size(source_adjacency)...)

        adjacency_site2site = spzeros(anz_rel, anz_rel)
        adjacency_site2site_compare = spzeros(anz_rel, anz_rel)
        for i in 1:anz_rel
            for j in (i+1):anz_rel
                if double_distance
                    res_mat, res_def = Haversine_is_connected_sources(i, j, distance_site_site, distance_source_site, material_connectivity_arr, SourceData; relevant_sites=relevant_sites)
                else
                    res_mat, res_def = Haversine_is_connected(i, j, distance_site_site, material_connectivity_arr; relevant_sites=relevant_sites)
                end
                adjacency_site2site[i, j] = res_mat ? 1.0 : 0.0
                adjacency_site2site_compare[i, j] = res_def ? 1.0 : 0.0
            end
        end
        adjacency_site2site += adjacency_site2site' 
        adjacency_site2site_compare += adjacency_site2site_compare'

        only_spatial = copy(adjacency_site2site)
        shared_sites_mat = source_adjacency * source_adjacency'
        adj_source2source = spzeros(Float64, anz_sources, anz_sources)
        
        for i in 1:anz_sources
            for j in i+1:anz_sources
                adj_source2source[i, j] = Haversine_is_connected_source2source(i, j, distance_source_source, shared_sites_mat, SourceDataDict[time_period][2], eff_distance=500.0, beta=25.0, double_distance=double_distance) ? 1.0 : 0.0
            end
        end
        adj_source2source += adj_source2source'
        
        if !double_distance
            for i in 1:anz_sources
                for j in 1:anz_rel
                    if adjacency_source2site[i,j] != 0
                        adjacency_source2site[i,j] = rand() < RayleighSISO(distance_source_site[i,j], source_maxdistance[i], eta=4, beta = 20) ? 1 : 0
                    end
                end
            end
        end

        Adj_total = [adjacency_site2site adjacency_source2site'; adjacency_source2site adj_source2source]
        Adj_compare = [adjacency_site2site_compare adjacency_source2site'; adjacency_source2site adj_source2source]
        distances_total = [distance_site_site distance_source_site'; distance_source_site distance_source_source]
        
        dropzeros!(Adj_total)
        dropzeros!(Adj_compare)

        only_spatial_total = [only_spatial zeros(Float64, size(adjacency_source2site')); zeros(Float64, size(adjacency_source2site)) zeros(Float64, size(adj_source2source))]
        
        write_results(joinpath(source_base, "only_spatial"), only_spatial_total, only_spatial_total, distances_total, ids)
        write_results(joinpath(source_base, "non_connected_components"), Adj_total, Adj_compare, distances_total, ids)
        
        Adj_total_non_conn = copy(Adj_total)
        Adj_compare_non_conn = copy(Adj_compare)
   
        connect_components!(Adj_total, components(Adj_total)..., distances_total)
        connect_components!(Adj_compare, components(Adj_compare)..., distances_total)

        write_results(joinpath(source_base, "connected_components"), Adj_total, Adj_compare, distances_total, ids)

        push!(graphen, [relevant_sites, Adj_total, Adj_total_non_conn, Adj_total .* distances_total, SourceData.Data, distances_total])
        push!(graphen, [relevant_sites, Adj_compare, Adj_compare_non_conn, Adj_compare .* distances_total, SourceData.Data, distances_total])
    end
    
    if plot_results
        Haversine_plot_connections_Sources(graphen, compare=Plot_compare, output_dir=output_dir, overwrite=overwrite)
    end
    
    return graphen
end


"""
    Haversine_plot_connections(graphen::AbstractArray; compare::Bool = true, output_dir::String = pwd())

Plots geographic lines representing connections between archaeological sites.
Black edges represent the default model, red edges represent new material-driven connectivity, and cyan indicates component unification.

# Arguments
- `graphen::AbstractArray`: Nested matrix array corresponding to graph evaluation steps.

# Keywords
- `compare::Bool`: If true, differences between models are highlighted via multi-color representation. Defaults to true.
- `output_dir::String`: Target export folder. Defaults to active user working location.

# Output
- `Nothing`: Executes physical system disk commits recording graphical maps.
"""
function Haversine_plot_connections(graphen::AbstractArray; compare::Bool = true, output_dir::String = pwd())
    
    for i in length(graphen)-1:-2:1
        graph = graphen[i] 
        graph2 = graphen[i+1] 
        time_period = graph[1][1].time_period

        out_dir = joinpath(output_dir, string(Dates.today()), "Haversine", "Plots", time_period)
        if isdir(out_dir) && !overwrite
            error("Output directory $out_dir already exists. Set `overwrite=true` to proceed.")
        else
            mkpath(out_dir)
        end

        @debug "Rendering connection plot for: $time_period"

        long = [entry.coords[1] for entry in graph[1]]
        lat = [entry.coords[2] for entry in graph[1]]

        left_border, right_border = extrema(long)
        bottom_border, top_border = extrema(lat)
        Δlong = 2.0; Δlat = 1.0;
        
        GMT.coast(
            region="$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+2.0*Δlong)/$(top_border+Δlat)+r", 
            proj=(name=:laea, center=[10, 50]), frame=:ag, res=:low, area=500, shore=:thin, borders=1,
            water="deepskyblue3@50", land="palegreen4@50",
        )
        
        if compare
            for knoten in axes(graph[2], 1)
                for adjacent_knoten in graph[2][:, knoten].nzind
                    if (adjacent_knoten in graph2[3][:, knoten].nzind) && (adjacent_knoten in graph[3][:, knoten].nzind)
                        GMT.lines!([graph[1][knoten].coords[1] graph[1][knoten].coords[2]; graph[1][adjacent_knoten].coords[1] graph[1][adjacent_knoten].coords[2]], pen=(0.5, :black))
                    elseif adjacent_knoten in graph[3][:, knoten].nzind
                        GMT.lines!([graph[1][knoten].coords[1] graph[1][knoten].coords[2]; graph[1][adjacent_knoten].coords[1] graph[1][adjacent_knoten].coords[2]], pen=(0.5, :red))
                    else
                        GMT.lines!([graph[1][knoten].coords[1] graph[1][knoten].coords[2]; graph[1][adjacent_knoten].coords[1] graph[1][adjacent_knoten].coords[2]], pen=(1.0, :cyan))
                    end
                end
            end
        else 
            for knoten in axes(graph[2], 1) 
                for adjacent_knoten in graph[2][:, knoten].nzind
                    GMT.lines!([graph[1][knoten].coords[1] graph[1][knoten].coords[2]; graph[1][adjacent_knoten].coords[1] graph[1][adjacent_knoten].coords[2]], pen=(0.5, :black))
                end
            end
        end 
        
        GMT.scatter!(long, lat, marker=:circle, markeredgecolor=0, size=0.1, markerfacecolor=:blue)
        
        legend_entries = (
            vspace = -0.25,
            header = (text="Archaeological graphs; Only site 2 site connections", font=(12, "Times-Bold")),
            hline1 = (pen=1, offset=0.5),
            symbol1 = (marker=:circle, size=0.4, dx_left=0.25, fill=:blue, dx_right=0.75, text="Site from $time_period"),
            hline2 = (pen=1, offset=0.5),
            symbol3 = (marker="-", size=0.6, dx_left=0.25, pen=(0.5, :black), dx_right=0.75, text="Default movement"),
            symbol4 = (marker="-", size=0.6, dx_left=0.25, pen=(0.5, :red), dx_right=0.75, text="Material connectivity"),
            symbol5 = (marker="-", size=0.6, dx_left=0.25, pen=(0.5, :cyan), dx_right=0.75, text="Connect components")
        )

        GMT.legend!(
            legend_entries,
            pos = (offset=(0.25, -0.25), anchor=:TL, width=10.0, justify=:BL, spacing=1.2),
            clearance = (0.25, 0.25), box = (pen=0.5, fill=:azure1),
            figsize = 10, show = false, format = :png, 
            figname = joinpath(out_dir, time_period)
        )
    end
    @info "Network plots generated successfully."
    return nothing
end


"""
    Haversine_plot_connections_Sources(graphen::AbstractArray; compare::Bool = true, output_dir::String = pwd(), overwrite::Bool = false)

Plots connections integrating archaeological sites and physical raw material sources.

# Arguments
- `graphen::AbstractArray`: Graph structure array containing sites, sources, and adjacency components.

# Keywords
- `compare::Bool`: Highlighting network variations based on materials versus base thresholds. Defaults to true.
- `output_dir::String`: Standard export path setup. Defaults to current working directory.
- `overwrite::Bool`: Grants delete and replace capabilities over the file structure. Defaults to false.

# Output
- `Nothing`: Persists visualization elements directly to disk (PDF execution).
"""
function Haversine_plot_connections_Sources(graphen::AbstractArray; compare::Bool = true, output_dir::String = pwd(), overwrite::Bool = false)


    for i in length(graphen)-1:-2:1
        graph = graphen[i] 
        graph2 = graphen[i+1] 
        time_period = graph[1][1].time_period

        out_dir = joinpath(output_dir, string(Dates.today()), "Haversine", "Sources", "Plots", time_period)
        if isdir(out_dir) && !overwrite
            error("Output directory $out_dir already exists. Set `overwrite=true` to proceed.")
        else
            mkpath(out_dir)
        end

        @debug "Rendering source connection plot for: $time_period"

        long1, lat1 = map(e -> (e.coords[1], e.coords[2]), graph[1]) |> x -> (getindex.(x, 1), getindex.(x, 2))
        long2, lat2 = map(s -> (s.coords[1], s.coords[2]), graph[5]) |> x -> (getindex.(x, 1), getindex.(x, 2))

        (left_border, right_border, bottom_border, top_border) = (1.89, 29.0, 44.9, 53.7)
        Δlong = 2.0; Δlat = 1.0; 

        GMT.coast(
            region="$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+2.0*Δlong)/$(top_border+Δlat)+r", 
            proj=(name=:laea, center=[10, 50]), frame=:ag, res=:low, area=500, shore=:thin, borders=1,
            water="deepskyblue3@50", land="palegreen4@50",
        )
        
        anz_sites = length(graph[1])
        anz_sources = length(graph[5])

        if compare
            for site in 1:anz_sites
                for adjacent_site in graph[2][site+1:anz_sites, site].nzind .+ site
                    if (adjacent_site in graph2[3][1:anz_sites, site].nzind) && (adjacent_site in graph[3][1:anz_sites, site].nzind)
                        GMT.lines!([graph[1][site].coords[1] graph[1][site].coords[2]; graph[1][adjacent_site].coords[1] graph[1][adjacent_site].coords[2]], pen=(0.5, :black))
                    elseif adjacent_site in graph[3][1:anz_sites, site].nzind
                        GMT.lines!([graph[1][site].coords[1] graph[1][site].coords[2]; graph[1][adjacent_site].coords[1] graph[1][adjacent_site].coords[2]], pen=(0.5, :red))
                    else
                        GMT.lines!([graph[1][site].coords[1] graph[1][site].coords[2]; graph[1][adjacent_site].coords[1] graph[1][adjacent_site].coords[2]], pen=(0.5, :cyan))
                    end
                end
            end
            
            for source in 1:anz_sources
                for adjacent_site in graph[2][1:anz_sites, anz_sites + source].nzind
                    if adjacent_site in graph[3][1:anz_sites, anz_sites + source].nzind
                        GMT.lines!([graph[5][source].coords[1] graph[5][source].coords[2]; graph[1][adjacent_site].coords[1] graph[1][adjacent_site].coords[2]], pen=(0.5, :pink))
                    else
                        GMT.lines!([graph[5][source].coords[1] graph[5][source].coords[2]; graph[1][adjacent_site].coords[1] graph[1][adjacent_site].coords[2]], pen=(0.5, :cyan))
                    end
                end
            end
            
            for source in 1:anz_sources
                for adjacent_source in graph[2][anz_sites+source+1:end, anz_sites+source].nzind .+ source
                    if adjacent_source in graph[3][anz_sites+1:end, anz_sites+source].nzind
                        GMT.lines!([graph[5][source].coords[1] graph[5][source].coords[2]; graph[5][adjacent_source].coords[1] graph[5][adjacent_source].coords[2]], pen=(0.5, :orange))
                    else
                        GMT.lines!([graph[5][source].coords[1] graph[5][source].coords[2]; graph[5][adjacent_source].coords[1] graph[5][adjacent_source].coords[2]], pen=(0.5, :cyan))
                    end
                end
            end
        else 
            for site in eachindex(graph[1])
                for adjacent_site in graph[2][site+1:anz_sites, site].nzind .+ site
                    GMT.lines!([graph[1][site].coords[1] graph[1][site].coords[2]; graph[1][adjacent_site].coords[1] graph[1][adjacent_site].coords[2]], pen=(0.5, :black))
                end
            end
            
            for source in eachindex(graph[5])
                for adjacent_site in graph[2][1:anz_sites, source+anz_sites].nzind
                    GMT.lines!([graph[5][source].coords[1] graph[5][source].coords[2]; graph[1][adjacent_site].coords[1] graph[1][adjacent_site].coords[2]], pen=(0.5, :red))
                end
            end

            for source in eachindex(graph[5])
                for adjacent_source in graph[2][anz_sites+source+1:end, anz_sites+source].nzind .+ source
                    GMT.lines!([graph[5][source].coords[1] graph[5][source].coords[2]; graph[5][adjacent_source].coords[1] graph[5][adjacent_source].coords[2]], pen=(0.5, :cyan))
                end
            end
        end

        GMT.scatter!(long2, lat2, fmt=:pdf, marker=:cross, markeredge=(1.5, :red), size=0.4, markerfacecolor=:red)
        GMT.scatter!(long1, lat1, fmt=:pdf, marker=:circle, markeredgecolor=0, size=0.1, markerfacecolor=:blue)
        
        base_filename = joinpath(out_dir, string(time_period))
        
        if compare
            legend_entries = (
                vspace  = -0.25,
                header  = (text="Types of Connections", font=(12, "Times-Bold")),
                hline1  = (pen=2, offset=0.5),
                vline   = (pen=2, offset=0),
                symbol1 = (marker=:circle, size=0.4, dx_left=0.25, fill=:blue, dx_right=0.75, text="Site"),
                symbol2 = (marker=:cross,  size=0.4, dx_left=0.25, fill=:red, dx_right=0.75, text="Source", pen=2),
                hline2  = (pen=2, offset=0.5),
                symbol3 = (marker="-", size=0.6, dx_left=0.25, pen=(1, :black), dx_right=0.75, text="Default-Movement"),
                symbol4 = (marker="-", size=0.6, dx_left=0.25, pen=(1, :red), dx_right=0.75, text="Material Connections"),
                symbol5 = (marker="-", size=0.6, dx_left=0.25, pen=(1, :pink), dx_right=0.75, text="Source to Site"),
                symbol6 = (marker="-", size=0.6, dx_left=0.25, pen=(1, :orange), dx_right=0.75, text="Source to Source"),
                symbol7 = (marker="-", size=0.6, dx_left=0.25, pen=(1, :cyan), dx_right=0.75, text="Connectivity-Auxiliary edges")
            )
        else
            legend_entries = (
                vspace  = -0.25,
                header  = (text="Archaeological graphs with Source 2 Source connections", font=(12, "Times-Bold")),
                hline1  = (pen=1, offset=0.5),
                ncolumns = 2,
                vline   = (pen=1, offset=0),
                symbol1 = (marker=:circle, size=0.4, dx_left=3.5, fill=:blue, dx_right=3.75, text="Site from $time_period"),
                symbol2 = (marker=:cross,  size=0.4, dx_left=3.5, fill=:red, dx_right=3.75, text="Source from $time_period", pen=2),
                hline2  = (pen=1, offset=0.5),
                symbol3 = (marker="-", size=0.6, dx_left=0.25, pen=(0.5, :black), dx_right=0.75, text="Site 2 Site"),
                symbol4 = (marker="-", size=0.6, dx_left=0.25, pen=(0.5, :red), dx_right=0.75, text="Source 2 Site"),
                symbol7 = (marker="-", size=0.6, dx_left=0.25, pen=(0.5, :cyan), dx_right=0.75, text="Source 2 Source")
            )
        end
        
        GMT.legend!(
            legend_entries,
            pos = (offset=(5.25, -1.0), anchor=:MR, width=5.25, justify=:BL, spacing=1.2),
            clearance = (0.25, 0.25), box = (pen=0.5, fill=:azure1), font = "12p,Times-Bold,black",
            figsize = "10c", format = :pdf, show = false,
            figname = base_filename * ".pdf"
        )
    end
    @info "Source network plots generated successfully."
    return nothing
end
