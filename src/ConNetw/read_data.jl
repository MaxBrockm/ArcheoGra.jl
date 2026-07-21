"""
    num_entries(dataset::Arch_DataSet)

Returns the total number of archaeological entries in the given dataset.

# Arguments
- `dataset::Arch_DataSet`: The archaeological dataset to evaluate.

# Output
- `Int`: The number of entries.
"""
function num_entries(dataset::Arch_DataSet)
    return length(dataset.Data)
end
@inline Base.length(dataset::Arch_DataSet) = length(dataset.Data)

"""
    num_entries(dataset::Arch_DataSet, property::Symbol, keyvalue::Any)

Returns the number of archaeological entries in the given dataset that match a specific property value.

# Arguments
- `dataset::Arch_DataSet`: The archaeological dataset to evaluate.
- `property::Symbol`: The property symbol to check (e.g., `:time_period`).
- `keyvalue::Any`: The specific value the property must match.

# Output
- `Int`: The number of matching entries.
"""
function num_entries(dataset::Arch_DataSet, property::Symbol, keyvalue::Any)
    return count(entry -> getproperty(entry, property) == keyvalue, dataset.Data)
end


"""
    add_data!(dataset::Arch_DataSet, data::ArchDataEntry)

Appends a new archaeological entry to an existing dataset.

# Arguments
- `dataset::Arch_DataSet`: The existing archaeological dataset.
- `data::ArchDataEntry`: The new archaeological entry to append.

# Output
- `Vector{ArchDataEntry}`: The updated underlying data vector.
"""
function add_data!(dataset::Arch_DataSet, data::ArchDataEntry)
    push!(dataset.Data, data)
end


"""
    get_coordinates(dataset::Arch_DataSet)

Retrieves a strictly typed list of coordinates and their corresponding time periods from the dataset.

# Arguments
- `dataset::Arch_DataSet`: The dataset containing archaeological entries.

# Output
- `Vector{Tuple{Tuple{Float64, Float64}, String}}`: A vector of tuples containing coordinates and time periods.
"""
function get_coordinates(dataset::Arch_DataSet)
    return Tuple{Tuple{Float64, Float64}, String}[(entry.coords, entry.time_period) for entry in dataset.Data]
end


"""
    check_presence(coords::Tuple{Float64, Float64}, time::String, dataset::Arch_DataSet)

Checks if a specific pair of coordinates for a given time period is already present in the dataset.

# Arguments
- `coords::Tuple{Float64, Float64}`: The coordinates to check as (longitude, latitude).
- `time::String`: The time period string.
- `dataset::Arch_DataSet`: The dataset to search.

# Output
- `Bool`: `true` if the coordinate and time period combination exists, `false` otherwise.
"""
function check_presence(coords::Tuple{Float64, Float64}, time::String, dataset::Arch_DataSet)
    return (coords, time) in get_coordinates(dataset)
end

"""
    get_material_id(material_presence::Vector{String}, materiallist::MaterialDataList)

Maps a list of material names to their corresponding internal IDs based on the material dataset.

# Arguments
- `material_presence::Vector{String}`: Vector of material names.
- `materiallist::MaterialDataList`: The dataset defining the materials and their IDs.

# Output
- `Vector{Int}`: Vector of corresponding material IDs.
"""
function get_material_id(material_presence::Vector{String}, materiallist::MaterialDataList)
    return [materiallist.Material_list[material].material_id for material in material_presence]
end


"""
    add_material!(entry::ArchDataEntry, material_presence::AbstractVector{<:AbstractString}, material_id::AbstractVector{<:Any}, material_multiplicity::AbstractVector{<:Any})

Functions adds material presence to an already existing entry in a dataset of archeological sites. 
Also adds material multiplicity of the corresponding materials. Function is required whenever an archeological sites 
is added at an already existing pair of coordinates.

# Arguments:
- `entry::ArchDataEntry`: Entry of archeological site where material presence will be added.
- `material_presence::AbstractVector{<:AbstractString}`: Vector of material names which are added to existing site.
- `material_id::AbstractVector{<:Any}`: List of IDs of materials which are added to existing site.
- `material_multiplicity::AbstractVector{<:Any}`: Vector of multiplicity of the materials which are added to the entry.
"""
function add_material!(entry::ArchDataEntry, material_presence::Vector{String}, material_id::Vector{Int}, material_multiplicity::Vector{Int})
    entry.material_presence = unique([entry.material_presence; material_presence]) ## update presence list with new materials
    for i in eachindex(material_id)
        if material_id[i] in entry.num_material_presence ## add more findings to existing ones
            entry.material_multiplicity[findfirst(x->x==material_id[i], entry.num_material_presence)] += material_multiplicity[i]
        else
            push!(entry.num_material_presence, material_id[i]) ## add new materials
            push!(entry.material_multiplicity, material_multiplicity[i])
        end
    end
end


"""
    read_data(file::file_props; data::Arch_DataSet, MaterialData::MaterialDataList)

Reads archaeological dataset from an external file (XLSX, CSV, or TXT) into the internal structures.
Updates the provided dataset and material lists dynamically.

# Arguments
- `file::file_props`: Configuration object defining the file path, type, and column mappings.

# Keywords
- `data::Arch_DataSet`: Existing dataset to append to (defaults to an empty dataset).
- `MaterialData::MaterialDataList`: Existing materials list (defaults to an empty list).

# Output
- `Tuple{Arch_DataSet, MaterialDataList}`: The updated archaeological dataset and material list.
"""
function read_data(file::file_props;
    data::Arch_DataSet=Arch_DataSet(Vector{ArchDataEntry}(undef, 0)), MaterialData::MaterialDataList=MaterialDataList(Dict{String, MaterialEntry}())
    )
    local data_arr::Matrix{Any}

    if file.file_type == "xlsx"
        data_arr = XLSX.readxlsx(file.filePath)[file.sheetname][:][:,:]
        replace!(data_arr, missing=>0)
        replace!(data_arr, " "=>0)
    elseif file.file_type == "txt"
        data_arr = DelimitedFiles.readdlm(file.filePath, file.delim, Any, '\n')
    elseif file.file_type == "csv"
        data_arr = DelimitedFiles.readdlm(file.filePath, file.delim)
    else 
        error("Invalid filetype specified: $(file.file_type). Supported types are 'xlsx', 'txt', and 'csv'.")
    end
    
    # Build list of materials 
    for (m, material) in enumerate(data_arr[1, file.material_cols])
        mat_key = material * file.material_age
        get!(MaterialData.Material_list, mat_key, MaterialEntry(length(MaterialData.Material_list) + 1, (file.material_age, file.material_distance[m])))
    end

    # Build a Set of existing (coords, time_period) for O(1) lookup
    existing_coords = Set{Tuple{Tuple{Float64,Float64}, String}}(
        (entry.coords, entry.time_period) for entry in data.Data
    )

    # Read data from file into storage
    for row in eachrow(data_arr[2:end,:])
        coords_key = (Tuple{Float64, Float64}(row[[file.long_ind, file.lat_ind]]), file.material_age)
        if coords_key ∉ existing_coords 
            # Only add new entries
            material_presence = Vector{String}(string.(data_arr[1,file.material_cols][row[file.material_cols].!=0]) .* file.material_age)
            name = !isnothing(file.name_ind) ? string(row[file.name_ind]) : ""
            class = !isnothing(file.class_ind) ? row[file.class_ind] : nothing
            mat_mult = Vector{Int}(row[file.material_cols][row[file.material_cols] .!= 0])

            add_data!(data, ArchDataEntry(
                id = length(data) + 1, 
                name = name, 
                coords = coords_key[1], 
                time_period = file.material_age, 
                material_presence = material_presence, 
                num_material_presence = get_material_id(material_presence, MaterialData), 
                material_multiplicity = mat_mult, 
                entry_classification = class, 
                grid_id = row[file.grid_id]
            ))
            push!(existing_coords, coords_key)
        else 
            # Update existing entries if found multiple times
            @debug "Found identical location at coordinates $(coords_key[1]) for period $(file.material_age). Updating entry." 
            ind = findfirst(x -> getproperty(x, :coords) == coords_key[1] && getproperty(x, :time_period) == file.material_age, data.Data)
            if isnothing(ind)
                # Safety fallback if set is desynced
                material_presence = Vector{String}(string.(data_arr[1, file.material_cols][row[file.material_cols] .!= 0]))
                name = !isnothing(file.name_ind) ? string(row[file.name_ind]) : ""
                class = !isnothing(file.class_ind) ? row[file.class_ind] : nothing
                mat_mult = Vector{Int}(row[file.material_cols][row[file.material_cols] .!= 0])
                add_data!(data, ArchDataEntry(
                    id = length(data) + 1, 
                    name = name, 
                    coords = coords_key[1], 
                    time_period = file.material_age, 
                    material_presence = material_presence, 
                    num_material_presence = get_material_id(material_presence, MaterialData), 
                    material_multiplicity = mat_mult,
                    entry_classification = class, 
                    grid_id = row[file.grid_id]
                ))
                push!(existing_coords, coords_key)
            else
                material_presence = Vector{String}(string.(data_arr[1, file.material_cols][row[file.material_cols] .!= 0]) .* file.material_age)
                add_material!(data.Data[ind], material_presence, get_material_id(material_presence, MaterialData), Vector{Int}(row[file.material_cols][row[file.material_cols] .!= 0]))
            end
        end
    end
    
    return data, MaterialData
end


"""
    read_data_TOML(config_file::String)

Reads a setup TOML file and creates the corresponding cumulative archaeological and material datasets.

When reading in data from the same time period across multiple files, it is recommended to use the identical `time_period` tag.

# Arguments
- `config_file::String`: The file path to the setup `.toml` file.

# Output
- `Tuple{Arch_DataSet, MaterialDataList}`: The fully initialized datasets parsed from the configuration.
"""
function read_data_TOML(config_file::String)
    ArchData = Arch_DataSet(Vector{ArchDataEntry}(undef, 0))
    MatData = MaterialDataList(Dict{String, MaterialEntry}())
    config_dir = dirname(abspath(config_file))

    for cfg in TOML.parsefile(config_file)
        file_path_raw = cfg[2]["filePath"]
        file_path = ""
        if isfile(joinpath(config_dir, file_path_raw))
            file_path = joinpath(config_dir, file_path_raw)
            @info "Using file: $file_path"
        elseif isfile(joinpath(dirname(@__DIR__), "src", "DataFolder", file_path_raw))
            file_path = joinpath(dirname(@__DIR__), "src", "DataFolder", file_path_raw)
            @info "Using relative file path: $file_path"
        else
            error("File not found: $file_path_raw. Please check the path in the TOML configuration.")
        end
        ArchData, MatData = read_data(
            file_props(
                filePath = file_path,
                file_type = cfg[2]["file_type"],
                sheetname = get(cfg[2], "sheetname", ""),
                delim = get(cfg[2], "delim", ",")[1],
                material_cols = cfg[2]["material_cols_start"]:cfg[2]["material_cols_end"],
                material_age = cfg[2]["material_age"],
                material_distance = cfg[2]["material_distance"] isa Number ? Float64(cfg[2]["material_distance"]) * ones(length(cfg[2]["material_cols_start"]:cfg[2]["material_cols_end"])) : Vector{Float64}(cfg[2]["material_distance"]),
                grid_id = cfg[2]["grid_id"],
                long_ind = cfg[2]["long_ind"],
                lat_ind = cfg[2]["lat_ind"],
                name_ind = get(cfg[2], "name_ind", nothing),
                class_ind = get(cfg[2], "class_ind", nothing)
            );
            data = ArchData, 
            MaterialData = MatData
        )
    end
    return ArchData, MatData
end


"""
    read_source_data(toml_file::String, MatData::MaterialDataList)

Parses source raw material data from a TOML configuration and updates the global material dictionary.

# Arguments
- `toml_file::String`: Path to the `.toml` configuration file.
- `MatData::MaterialDataList`: The existing material dataset to be updated.

# Output
- `Dict{String, Tuple{Source_DataSet, Float64}}`: Dictionary mapping time periods to their source dataset and multiplier.
"""
function read_source_data(toml_file::String, MatData::MaterialDataList)
    SourceDict = Dict{String, Tuple{Source_DataSet, Float64}}()
    config_dir = dirname(abspath(toml_file))
    
    for cfg in TOML.parsefile(toml_file)
        SourceData = Source_DataSet(Vector{SourceDataEntry}(undef, 0), 0.0)
        
        mult_global = get(cfg[2], "source_data_mult", nothing)
        file_type = get(cfg[2], "source_file_type", "txt")
        source_path_raw = cfg[2]["source_data_path"]
        source_path = isabspath(source_path_raw) ? source_path_raw : normpath(joinpath(config_dir, source_path_raw))
        
        local data_arr::Matrix{Any}
        if file_type == "txt"
            data_arr = DelimitedFiles.readdlm(source_path, cfg[2]["source_data_delim"][1], Any, '\n')
        else
            data_arr = XLSX.readxlsx(source_path)[cfg[2]["source_sheetname"]][:][:,:]
        end
        replace!(data_arr, missing => 0)
        replace!(data_arr, " " => 0)

        for row in eachrow(data_arr[2:end, :])
            source_mult = isnothing(mult_global) ? 1.0 : Float64(mult_global)

            mat_name = string(row[cfg[2]["source_data_name_col"]]) * string(cfg[2]["material_age"])
            
            push!(SourceData.Data, SourceDataEntry(
                id = length(SourceData.Data) + 1,
                material_name = mat_name,
                coords = (Float64(row[cfg[2]["source_data_long_col"]]), Float64(row[cfg[2]["source_data_lat_col"]])),
                grid_id = row[cfg[2]["source_grid_id"]]
            ))
            
            # Update Material Dictionary if needed
            if !haskey(MatData.Material_list, mat_name)
                @debug "Adding new source material: $mat_name with distance $(cfg[2]["material_distance"]) to material list."
                MatData.Material_list[mat_name] = MaterialEntry(length(MatData.Material_list) + 1, (string(cfg[2]["material_age"]), 0.0))
            end
            SourceDict[cfg[2]["material_age"]] = (SourceData, source_mult)
        end
    end
    return SourceDict
end


# ==============================================================================
# Plotting Functions
# ==============================================================================

"""
    plot_sites(toml_file::String; output_dir::String = pwd(), overwrite::Bool = false)

Generates a map plot of ONLY the archaeological sites given in the corresponding TOML file.

# Arguments
- `toml_file::String`: Name or path of the setup `.toml` file.

# Keywords
- `output_dir::String`: Base directory to save the plots (defaults to current working directory).
- `overwrite::Bool`: If `false` (default), throws an error if the directory exists.
"""
function plot_sites(toml_file::String; output_dir::String = pwd(), overwrite::Bool = false)
    # Check if we should use the package's DataFolder or a user-provided absolute path
    file_path = isfile(toml_file) ? toml_file : joinpath(@__DIR__, "DataFolder", toml_file)
    data, mats = read_data_TOML(file_path)
    prep_data = prepare_network_data(data, mats)
    
    for (time_period, prep) in prep_data
        @debug "Plotting sites for time period: $time_period"
        
        out_path = joinpath(output_dir, string(Dates.today()), time_period)
        if isdir(out_path)
            if !overwrite
                error("Folder '$out_path' already exists! To overwrite files, use the argument `overwrite=true`.")
            else
                @info "Folder already exists. Plots will be overwritten (overwrite=true)."
            end
        else
            mkpath(out_path)
        end
        
        relevant_sites = prep.sites
        long = [entry.coords[1] for entry in relevant_sites]
        lat = [entry.coords[2] for entry in relevant_sites]

        left_border, right_border = extrema([p.coords[1] for p in relevant_sites])
        bottom_border, top_border = extrema([p.coords[2] for p in relevant_sites])

        Δlong = 2.0
        Δlat = 1.0

        GMT.coast(
            region="$(left_border - Δlong)/$(bottom_border - Δlat)/$(right_border + 2.0 * Δlong)/$(top_border + Δlat)+r", 
            proj=(name=:laea, center=[10, 50]), frame=:ag, res=:low, area=500, shore=:thin, borders=1,
            water="deepskyblue3@50", land="palegreen4@50"
        )
        GMT.scatter!(
            long, lat, mt=:png, marker=:circle,
            markeredgecolor=0, size=0.1, markerfacecolor=:blue, show=true, fmt=:png, dpi=1200,
            figname=joinpath(out_path, "$(time_period)_sites")
        )
    end
    return nothing
end

"""
    plot_material(toml_file::String; output_dir::String = pwd(), overwrite::Bool = false)

Generates a map plot visualizing the distribution and multiplicity of materials across sites.

# Arguments
- `toml_file::String`: Name or path of the setup `.toml` file.

# Keywords
- `output_dir::String`: Base directory to save the plots (defaults to current working directory).
- `overwrite::Bool`: If `false` (default), throws an error if the directory exists.
"""
function plot_material(toml_file::String; output_dir::String = pwd(), overwrite::Bool = false)
    file_path = isfile(toml_file) ? toml_file : joinpath(@__DIR__, "DataFolder", toml_file)
    data, mats = read_data_TOML(file_path)
    prep_data = prepare_network_data(data, mats)
    
    for (time_period, prep) in prep_data
        @debug "Plotting materials for time period: $time_period"
        
        out_path = joinpath(output_dir, string(Dates.today()), time_period)
        if isdir(out_path)
            if !overwrite
                error("Folder '$out_path' already exists! To overwrite files, use the argument `overwrite=true`.")
            else
                @info "Folder already exists. Plots will be overwritten (overwrite=true)."
            end
        else
            mkpath(out_path)
        end
        
        relevant_sites = prep.sites
        
        # Default map constraints
        (left_border, right_border, bottom_border, top_border) = (1.89, 29.0, 44.9, 53.7)
        Δlong = 2.0
        Δlat = 1.0

        GMT.coast(
            region="$(left_border - Δlong)/$(bottom_border - Δlat)/$(right_border + 2.0 * Δlong)/$(top_border + Δlat)+r", 
            proj=(name=:laea, center=[10, 50]), frame=:ag, res=:low, area=500, shore=:thin, borders=1,
            water="deepskyblue3@50", land="palegreen4@50",
        )
        
        for site in relevant_sites[1:end-1]
            num_mats = sum(site.material_multiplicity)
            cmcolor = num_mats == 0 ? "black" : "red" 
            GMT.scatter!(
                site.coords[1], site.coords[2], marker=:circle,
                markeredgecolor=0, size=(0.1 + num_mats * 0.01), markerfacecolor=cmcolor
            )
        end
        
        # Process the final point and trigger the show/save behavior
        num_mats = sum(relevant_sites[end].material_multiplicity)
        cmcolor = num_mats == 0 ? "black" : "red" 
        GMT.scatter!(
            relevant_sites[end].coords[1], relevant_sites[end].coords[2], marker=:circle,
            markeredgecolor=0, size=(0.1 + num_mats * 0.01), markercolor=cmcolor, show=true,
            figname=joinpath(out_path, "$(time_period)_materials")
        )
    end
    return nothing
end