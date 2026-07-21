"""
    compare_all_is_connected(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, material_weight_matrices::Dict, relevant_sites::Vector{ArchDataEntry})

Evaluates connections between sites based on spatial distance and incrementally added material data. 
Instead of a binary no-material/all-material comparison, it saves each incremental step of additional materials.

# Arguments
- `start_node::Int`: Index of the primary site.
- `end_node::Int`: Index of the secondary site.
- `distance_matrix::AbstractArray{Float64}`: Matrix containing spatial distances.
- `material_weight_matrices::Dict`: Dictionary of material matrices and distribution ranges.
- `relevant_sites::Vector{ArchDataEntry}`: List of valid archaeological sites.

# Output
- `Vector{Bool}`: A boolean vector representing connectivity at each stage. The first entry is the default connection; the second is the group-adjusted connection; subsequent entries represent added materials.
"""
function compare_all_is_connected(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, material_weight_matrices::Dict, relevant_sites::Vector{ArchDataEntry})
    r = rand()
    connectivity_vector = zeros(Bool, length(material_weight_matrices)+2)

    eff_max_distance = maximum_distance 
    
    Prob = RayleighSISO(distance_matrix[start_node,end_node], eff_max_distance; eta=2)
    connectivity_vector[1] = r < Prob

    if !isnothing(relevant_sites) && (relevant_sites[start_node].entry_classification == relevant_sites[end_node].entry_classification)
        eff_max_distance *= 2
        Prob = RayleighSISO(distance_matrix[start_node,end_node], eff_max_distance; eta=2)
    end
    connectivity_vector[2] = r < Prob

    for (i,(name,mat)) in enumerate(material_weight_matrices)
        Prob += (1-Prob)*RayleighSISO(distance_matrix[start_node,end_node], mat[2], eta=2+material_multiplier*mat[1][start_node, end_node])
        connectivity_vector[2+i] = r < Prob
    end
    return connectivity_vector
end


"""
    compare_all_connect_nodes(data::Arch_DataSet, MaterialData::MaterialDataList; prepared::Union{Dict{String,PreparedTimePeriod},Nothing} = nothing, output_dir::String = pwd(), overwrite::Bool = false)

Generates all sub-stage graphs using increasing amounts of materials to trace network evolution.

# Arguments
- `data::Arch_DataSet`: The archaeological dataset.
- `MaterialData::MaterialDataList`: The material properties list.

# Keywords
- `prepared::Union{Dict{String,PreparedTimePeriod},Nothing}`: Pre-indexed memory structures. Defaults to nothing.
- `output_dir::String`: Target export path. Defaults to the current working directory.
- `overwrite::Bool`: Grants file overwrite permissions. Defaults to false.

# Output
- `Vector{Any}`: A collection detailing adjacency matrices, time periods, and site vectors.
"""
function compare_all_connect_nodes(data::Arch_DataSet, MaterialData::MaterialDataList; prepared::Union{Dict{String,PreparedTimePeriod},Nothing} = nothing, output_dir::String = pwd(), overwrite::Bool = false)
    prep_data = isnothing(prepared) ? prepare_network_data(data, MaterialData) : prepared
    time_periods = collect(keys(prep_data))
    graphen = []

    for time_period in time_periods
        prep = prep_data[time_period]

        out_base = joinpath(output_dir, string(Dates.today()), "compare_all_Mat", time_period)
        if isdir(out_base) && !overwrite
            error("Output directory $out_base already exists. Set `overwrite=true` to proceed.")
        else
            mkpath(out_base)
        end

        material_connectivity_arr = adjacency_data(MaterialData, prep)
        relevant_sites = prep.sites
        anz_rel = length(relevant_sites)
        
        all_materials_adjacency = [spzeros(Int,anz_rel, anz_rel) for _ in 1:length(material_connectivity_arr)+2]
        distance_matrix = distance_allocation(relevant_sites)

        for i in 1:anz_rel
            for j in (i+1):anz_rel
                connec_vec_all_mats = compare_all_is_connected(i, j, distance_matrix, material_connectivity_arr, relevant_sites)
                for m in eachindex(connec_vec_all_mats)
                    if connec_vec_all_mats[m]
                        all_materials_adjacency[m][i, j] = 1.0
                        all_materials_adjacency[m][j, i] = 1.0
                    end
                end
            end
        end

        for m in eachindex(all_materials_adjacency)
            connect_components!(all_materials_adjacency[m], components(all_materials_adjacency[m])..., distance_matrix)
        end
        
        for (i, entry) in enumerate(all_materials_adjacency)
            writedlm(joinpath(out_base, "$(time_period)_stage$i.csv"), entry, ',')
        end

        push!(graphen, [all_materials_adjacency, time_period, relevant_sites, distance_matrix])
    end
    
    return graphen
end


"""
    compare_all_is_connected_sources(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, distance_matrix_sources::AbstractArray{Float64}, material_weight_matrices::Dict, sourcedata::Source_DataSet; relevant_sites::Vector{ArchDataEntry})

Source-aware version of `compare_all_is_connected`. Records incremental connectivity states 
including default distance, group-adjusted distance, and stages for each uniquely added material via specific sources.

# Output
- `Vector{Bool}`: Step-by-step connectivity flags.
"""
function compare_all_is_connected_sources(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, distance_matrix_sources::AbstractArray{Float64}, material_weight_matrices::Dict, sourcedata::Source_DataSet; relevant_sites::Vector{ArchDataEntry})
    r = rand()
    eff_max_distance = maximum_distance 
    
    if !isnothing(relevant_sites) && (relevant_sites[start_node].entry_classification == relevant_sites[end_node].entry_classification)
        eff_max_distance *= 2
    end

    Prob = RayleighSISO(distance_matrix[start_node,end_node], eff_max_distance; eta=2)
    connec_vec = [r < Prob]

    unique_materials = sort(unique([s.material_name for s in sourcedata.Data]))

    for mat_name in unique_materials
        if haskey(material_weight_matrices, mat_name)
            mat = material_weight_matrices[mat_name]
            if mat[1][start_node,end_node] > 0
                for (i,s) in enumerate(sourcedata.Data)
                    if s.material_name == mat_name
                        Prob += (1-Prob)*RayleighSISO(distance_matrix_sources[i,start_node]+distance_matrix_sources[i,end_node], mat[2], eta=2+material_multiplier*mat[1][start_node, end_node])
                    end
                end
            end
        end
        push!(connec_vec, r < Prob)
    end

    return connec_vec
end


"""
    compare_all_connect_nodes_sources(data::Arch_DataSet, MaterialData::MaterialDataList, SourceDataDict::Dict; double_distance::Bool = true, output_dir::String = pwd(), overwrite::Bool = false)

Creates source-inclusive spatial networks for all sub-stages with incremental materials.
Stage 1 serves as the spatial default, stage 2 includes group connectivity, and subsequent stages integrate material connections one by one.

# Output
- `Vector{Any}`: A collection of block matrices bridging sites and sources.
"""
function compare_all_connect_nodes_sources(data::Arch_DataSet, MaterialData::MaterialDataList, SourceDataDict::Dict; double_distance::Bool = true, output_dir::String = pwd(), overwrite::Bool = false)
    time_periods = unique([entry.time_period for entry in data.Data])
    graphen = []
    for time_period in time_periods
        @info "Evaluating incremental source connectivity for time period: $time_period"
        SourceData = SourceDataDict[time_period][1] 
        
        out_base = output_dir == pwd() ? joinpath(output_dir, string(Dates.today()), "compare_all_Mat_Sources", time_period) : output_dir
        if isdir(out_base) && !overwrite
            error("Output directory $out_base already exists. Set `overwrite=true` to proceed.")
        else
            mkpath(out_base)
        end

        source_adjacency, source_maxdistance = adjacency_sources(data, SourceData, MaterialData, time_period)
        material_connectivity_arr = adjacency_data(MaterialData, data, time_period)
        
        relevant_sites = data.Data[[entry.time_period==time_period for entry in data.Data]]
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

        # Create Source-to-Source connections
        shared_sites_mat = source_adjacency*source_adjacency'
        adj_source2source = spzeros(Bool, anz_sources, anz_sources)
        for i in 1:anz_sources
            for j in i+1:anz_sources
                adj_source2source[i,j] = Haversine_is_connected_source2source(i,j, distance_source_source, shared_sites_mat, SourceDataDict[time_period][2], eff_distance=500, beta=25, double_distance=double_distance)
            end
        end
        adj_source2source += adj_source2source'
        
        # Optional probabilistic reduction for non-double_distance mode
        if !double_distance
            for i in 1:anz_sources
                for j in 1:anz_rel
                    if adjacency_source2site[i,j] != 0
                        adjacency_source2site[i,j] = rand() < RayleighSISO(distance_source_site[i,j], source_maxdistance[i], eta=4, beta=20) ? 1 : 0
                    end
                end
            end
        end

        dummy_vec = compare_all_is_connected_sources(1, 2, distance_site_site, distance_source_site, material_connectivity_arr, SourceData; relevant_sites=relevant_sites)
        num_stages = length(dummy_vec)
        @debug "Number of evaluated material stages: $num_stages"
        
        all_site2site_adjacency = [spzeros(Int, anz_rel, anz_rel) for _ in 1:num_stages]

        for i in 1:anz_rel
            for j in (i+1):anz_rel
                connec_vec_all_mats = compare_all_is_connected_sources(i, j, distance_site_site, distance_source_site, material_connectivity_arr, SourceData; relevant_sites=relevant_sites)
                
                for m in 1:num_stages
                    if connec_vec_all_mats[m]
                        all_site2site_adjacency[m][i,j] = 1
                        all_site2site_adjacency[m][j,i] = 1
                    end
                end
            end
        end

        distances_total = [distance_site_site distance_source_site'; distance_source_site distance_source_source]
        stage_matrices = []

        for m in 1:num_stages
            Adj_total = [all_site2site_adjacency[m] adjacency_source2site'; adjacency_source2site adj_source2source]
            connect_components!(Adj_total, components(Adj_total)..., distances_total)
            writedlm(joinpath(out_base, "$(time_period)_stage$m.csv"), Adj_total, ',')
            push!(stage_matrices, Adj_total)
        end

        # Store the collection of block matrices for this time period
        push!(graphen, [stage_matrices, time_period, relevant_sites, SourceData.Data, distances_total, ids])
    end
    
    return graphen
end