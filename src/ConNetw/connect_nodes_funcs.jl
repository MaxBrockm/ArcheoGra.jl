"""
    adjacency_data(Material_data::MaterialDataList, data::Arch_DataSet, time_period::String)

Reads in the shared raw material distribution of vertices for a given time period. 
Connectivity is stored in a weighted adjacency matrix, where the weight corresponds to the cumulative amount of shared material presence.

# Arguments
- `Material_data::MaterialDataList`: List of all materials and their corresponding properties.
- `data::Arch_DataSet`: List of archaeological sites with their properties.
- `time_period::String`: Time period in which the material connectivity is compared.

# Output
- `Dict{String, Tuple{SparseMatrixCSC{Float64, Int}, Float64}}`: A dictionary mapping the material name to a tuple of its weighted adjacency matrix and its distribution distance.
"""
function adjacency_data(Material_data::MaterialDataList, data::Arch_DataSet, time_period::String)
    num_sites = num_entries(data, :time_period, time_period)
    timed_sites = data.Data[[getproperty(entry, :time_period)==time_period for entry in data.Data]]
    
    mat_adjacency_dict = Dict{String, Tuple{AbstractArray,Number}}()
    
    for (name,mat) in Material_data.Material_list 
        if mat.mat_prop[1] == time_period
            A = spzeros(num_sites, num_sites)
            relevant_sites = Tuple{Int,Int}[]

            for (i, entry) in enumerate(timed_sites)
                idx = findfirst(==(mat.material_id), entry.num_material_presence)
                if !isnothing(idx)
                    push!(relevant_sites, (i, entry.material_multiplicity[idx]))
                end
            end

            for a in eachindex(relevant_sites)
                for b in (a+1):length(relevant_sites)
                    si, mi = relevant_sites[a]
                    sj, mj = relevant_sites[b]
                    A[si, sj] = mi + mj
                end
            end
            mat_adjacency_dict[name] = (A+A', mat.mat_prop[2])
        end
    end

    return mat_adjacency_dict
end

"""
    adjacency_data(Material_data::MaterialDataList, prep::PreparedTimePeriod)

Optimized overload using pre-indexed data from `PreparedTimePeriod`.
Material-to-site mappings are read directly from `prep.material_map`, eliminating linear `findfirst` calls.

# Arguments
- `Material_data::MaterialDataList`: List of all materials and their corresponding properties.
- `prep::PreparedTimePeriod`: Pre-indexed views of the archaeological sites for a specific time period.

# Output
- `Dict{String, Tuple{SparseMatrixCSC{Float64, Int}, Float64}}`: A dictionary mapping the material name to a tuple of its weighted adjacency matrix and its distribution distance.
"""
function adjacency_data(Material_data::MaterialDataList, prep::PreparedTimePeriod)
    num_sites = length(prep.sites)
    mat_adjacency_dict = Dict{String, Tuple{AbstractArray,Number}}()

    for (name, mat) in Material_data.Material_list
        mat.mat_prop[1] == prep.time_period || continue
        mid = mat.material_id
        bucket = get(prep.material_map, mid, Tuple{Int,Int}[])
        isempty(bucket) && continue

        A = spzeros(num_sites, num_sites)
        for a in eachindex(bucket)
            for b in (a+1):length(bucket)
                si, mi = bucket[a]
                sj, mj = bucket[b]
                A[si, sj] = mi + mj
            end
        end
        mat_adjacency_dict[name] = (A + A', mat.mat_prop[2])
    end
    return mat_adjacency_dict
end


"""
    RayleighSISO(r::Number, r0::Number; beta::Number = 10.0, eta::Number = 2.0)

Evaluates the connectedness probability using a Rayleigh (SISO) mathematical model: 
\$R(x)=e^{-\\beta(\\frac{r}{r_0})^{\\eta}}\$

# Arguments
- `r::Number`: The current measured distance.
- `r0::Number`: The normalization or baseline maximum distance limit.

# Keywords
- `beta::Number`: The decay scaling parameter. Defaults to 10.0.
- `eta::Number`: The exponent modifying distance penalization. Defaults to 2.0.

# Output
- `Float64`: The computed connection probability between 0.0 and 1.0.
"""
@inline function RayleighSISO(r::Number, r0::Number; beta::Number = 10, eta::Number = 2)
    return exp(-beta*(r/r0)^eta)
end

"""
    connect_nodes(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, material_weight_matrices::Dict, relevant_sites::Vector{ArchDataEntry}; symmetric::Bool = true)

Evaluates whether two nodes should be connected based on their geographic distance, group classification, and shared materials.
Generates stacked probabilities by utilizing the `RayleighSISO` function, increasing the likelihood of connection for each shared material.

# Arguments
- `start_node::Int`: Index of the primary site.
- `end_node::Int`: Index of the secondary site.
- `distance_matrix::AbstractArray{Float64}`: Matrix containing distances between all nodes.
- `material_weight_matrices::Dict`: Dictionary containing matrices for each material tracking the amount of shared items.
- `relevant_sites::Vector{ArchDataEntry}`: Evaluated sites to access entry properties like group classifications.

# Keywords
- `symmetric::Bool`: Set to `true` to use symmetrical distances, or `false` to evaluate directional asymmetry. Defaults to `true`.

# Output
- `Tuple`: Returns boolean flags representing the material-enhanced connection and the default spatial connection. Returns 4 booleans if asymmetric, 2 if symmetric.
"""
function connect_nodes(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, material_weight_matrices::Dict, relevant_sites::Vector{ArchDataEntry}; symmetric :: Bool = true)
    r = rand()
    eff_max_distance = maximum_distance 
    
    if !isnothing(relevant_sites) && (relevant_sites[start_node].entry_classification == relevant_sites[end_node].entry_classification)
        eff_max_distance *= 2
    end

    ## Default connection
    ProbThere  = RayleighSISO(distance_matrix[start_node,end_node], eff_max_distance; eta=2)
    ProbThere_default = copy(ProbThere)

    
    for (_, mat) in material_weight_matrices
        mat[1][start_node, end_node] == 0.0 && continue
        ProbThere += (1.0 - ProbThere) * RayleighSISO(distance_matrix[start_node, end_node], mat[2], eta=2.0 + material_multiplier * mat[1][start_node, end_node])
    end

    if !symmetric
        ProbBack  = RayleighSISO(distance_matrix[end_node,start_node], eff_max_distance; eta=2)
        ProbBack_default = copy(ProbBack)

        for (_, mat) in material_weight_matrices
            (mat[1][end_node,start_node] == 0) && (continue)
            ProbBack += (1-ProbBack)*RayleighSISO(distance_matrix[end_node,start_node], mat[2], eta=2+material_multiplier*mat[1][end_node, start_node])
        end
        return  r < ProbThere, r < ProbThere_default, r < ProbBack, r < ProbBack_default
    else
        return  r < ProbThere, r < ProbThere_default
    end
end


"""
    connect_nodes(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, symmetric::Bool = true; distance_matrix_back::AbstractArray{Float64} = distance_matrix)

Variant of `connect_nodes` strictly designed to evaluate connections across time periods based exclusively on geographic limits without integrating material modifiers.

# Arguments
- `start_node::Int`: Index of the primary site.
- `end_node::Int`: Index of the secondary site.
- `distance_matrix::AbstractArray{Float64}`: Matrix containing spatial distances.
- `symmetric::Bool`: Set to `true` to use symmetrical distances. Defaults to `true`.

# Keywords
- `distance_matrix_back::AbstractArray{Float64}`: Required if `symmetric` is `false` to determine reverse distances. Defaults to `distance_matrix`.

# Output
- `Tuple`: Returns boolean flags identifying whether nodes satisfy the spatial probability threshold.
"""
function connect_nodes(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, symmetric::Bool = true; distance_matrix_back::AbstractArray{Float64} = distance_matrix)
    r = rand()
    eff_max_distance = Float64(maximum_distance)

    ProbThere  = RayleighSISO(distance_matrix[start_node,end_node], eff_max_distance; eta=2)

    if !symmetric
        ProbBack  = RayleighSISO(distance_matrix_back[end_node,start_node], eff_max_distance; eta=2)
        return  r < ProbThere, r < ProbBack
    else
        return  r < ProbThere
    end
end


"""
    adjacency_sources(data::Arch_DataSet, SourceData::Source_DataSet, MaterialData::MaterialDataList, time_period::String)

Creates an adjacency matrix quantifying the presence of materials supplied from known sources at specific archaeological sites.

# Arguments
- `data::Arch_DataSet`: List of archaeological sites.
- `SourceData::Source_DataSet`: Dataset containing verified raw material sources.
- `MaterialData::MaterialDataList`: List of valid materials.
- `time_period::String`: Specific historical period to evaluate.

# Output
- `Tuple{SparseMatrixCSC{Float64, Int}, Vector{Float64}}`: The source-to-site distribution matrix and the corresponding distribution distance array for each source.
"""
function adjacency_sources(data::Arch_DataSet, SourceData::Source_DataSet, MaterialData::MaterialDataList, time_period::String)
    timed_sites = data.Data[[getproperty(entry, :time_period)==time_period for entry in data.Data]]    
    SourceMat = spzeros(length(SourceData.Data), length(timed_sites))
    material_distance_per_source = zeros(length(SourceData.Data))
    for (s,source) in enumerate(SourceData.Data)
        mat_id = MaterialData.Material_list[source.material_name].material_id 
        for (v,site) in enumerate(timed_sites)
            idx = findfirst(==(mat_id), site.num_material_presence)
            SourceMat[s, v] = !isnothing(idx) && site.material_multiplicity[idx]
        end
        material_distance_per_source[s] = MaterialData.Material_list[source.material_name].mat_prop[2]
    end
    return SourceMat, material_distance_per_source
end

"""
    adjacency_sources(prep::PreparedTimePeriod, SourceData::Source_DataSet, MaterialData::MaterialDataList)

Optimized overload of `adjacency_sources` using pre-indexed data structures to process material distributions.
Per-site material lookups execute in O(1) performance using `site_material_dicts`.

# Arguments
- `prep::PreparedTimePeriod`: Pre-indexed archaeological dataset views.
- `SourceData::Source_DataSet`: Dataset containing physical sources.
- `MaterialData::MaterialDataList`: List of associated materials.

# Output
- `Tuple{SparseMatrixCSC{Float64, Int}, Vector{Float64}}`: The source-to-site distribution matrix and the corresponding distribution distance array for each source.
"""
function adjacency_sources(prep::PreparedTimePeriod, SourceData::Source_DataSet, MaterialData::MaterialDataList)
    num_sites = length(prep.sites)
    SourceMat = spzeros(length(SourceData.Data), num_sites)
    material_distance_per_source = zeros(length(SourceData.Data))
    for (s, source) in enumerate(SourceData.Data)
        mat_id = MaterialData.Material_list[source.material_name].material_id
        for (v, smd) in enumerate(prep.site_material_dicts)
            mult = get(smd, mat_id, 0)
            if mult != 0
                SourceMat[s, v] = mult
            end
        end
        material_distance_per_source[s] = MaterialData.Material_list[source.material_name].mat_prop[2]
    end
    return SourceMat, material_distance_per_source
end