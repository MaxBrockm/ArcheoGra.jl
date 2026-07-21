"""
    archeoLCPA(data::AbstractArray{ArchDataEntry}, k::Int, ℓ::Int; PlotOut::Bool=false, conn_func::ConnectionFunction = Tobler_Hiking_adjusted)

Function to setup LCPA for the current dataset. Takes the correct part of the topographic data for the datapoints and compresses it to the desires size (`k` and `ℓ`). The corresponding adjacency matrix for the Dijkstra Algorithm is created and connections removed for waterways. Afterwards the LCPA-algorithm is called and all paths and their lengths between nodes are calculated and returned.

# Arguments:
- `data::AbstractArray{ArchDataEntry}`: List of all archeological sites with their properties
- `k::Int`: Number of topography-pixels which are compressed into one -- width
- `l::Int`: Number of topography-pixels which are compressed into one -- heigth
- `PlotOut::Bool=false`: Parameter to decide wether to plot all least cost paths between all nodes. Might take very ong if enabled.
- `conn_func::ConnectionFunction = Tobler_Hiking_adjusted`: Choice of cost-function in LCPA. Standard is modified Tobler Hiking function.
"""
function archeoLCPA(data::AbstractArray{ArchDataEntry}, k::Int, ℓ::Int; PlotOut::Bool=false, conn_func::ConnectionFunction = Tobler_Hiking_adjusted)
    ## load topographic data
    top_data = tiff_to_matrix(joinpath("src", "DataFolder", "-9_35_32_54.tiff"))
    
    ## find correct map section for given data
    left_border, right_border = extrema([p.coords[1] for p in data])
    bottom_border, top_border = extrema([p.coords[2] for p in data])

    long_p_px = 41/2460 ## specific to this topo. dataset
    lat_p_px = 19/1140

    lx = floor(Int, (left_border+9)/long_p_px)
    ux = ceil(Int, (right_border+9)/long_p_px)
    ly = floor(Int, (bottom_border-35)/lat_p_px)
    uy = ceil(Int, (top_border-35)/lat_p_px)
    
    top_data = top_data[ly-10:uy+10, lx-10:ux+10]
    
    new_lb_corner = (-9+lx*long_p_px, 35+ly*lat_p_px); new_rt_corner = (-9+ux*long_p_px, 35+uy*lat_p_px)
  
    useTopdata = TopographicalData(top_data, size(top_data), new_lb_corner, new_rt_corner, (1,1))

    ## compress topographic data into smaller size for LCPA
    compressed_topData = compress(useTopdata, k, ℓ)

    ## remove unwanted connections
    ind_list = lexi2ind.([coords2lexi(data[i].coords, compressed_topData) for i in eachindex(data)], [compressed_topData])
    Ad_water = weighted_adjacency_matrix_haversine(compressed_topData, conn_func.f)

    Ad = delConnections(Ad_water, lexiMat2Vec(WaterLandBool(compressed_topData, ind_list)))

    ## determine LCPA 
    if conn_func.symmetric ### Currently placeholder for when all_pairs_LCPA_GMT gets the option for symmetry
        return all_pairs_LCPA_GMT(compressed_topData, Ad, data, PlotOutput=PlotOut)
    else 
        return all_pairs_LCPA_GMT(compressed_topData, Ad, data, PlotOutput=PlotOut)
    end
end

"""
    all_pairs_LCPA_GMT(TopData::TopographicalData, Ad::AbstractMatrix, sites::Vector{ArchDataEntry}; PlotOutput = false)

Slightly modified version for the LCPA calculation between all archaeological sites. Adjustments made for different data types and additional option to choose whether to plot or not.
"""
function all_pairs_LCPA_GMT(TopData::TopographicalData, Ad::AbstractMatrix, sites::Vector{ArchDataEntry}; PlotOutput = false)
    start_lexi_list = [coords2lexi(sites[i].coords, TopData) for i in eachindex(sites)]
    
    if PlotOutput
    GMT.coast(
        region="$(TopData.lBot[1])/$(TopData.lBot[2])/$(TopData.rTop[1])/$(TopData.rTop[2])+r6",
        proj=(name=:laea, center=[10,50]), 
        frame=:ag, 
        res=:low,
        area=500, 
        shore=:thin, 
        borders=1,
        rivers=:deepskyblue3,
        # show = true,
        water="deepskyblue3@50", land="palegreen4@50",
    )
    end
    length_list = length(start_lexi_list)

    path_save = Array{Any}(nothing, length_list,length_list)
    length_save = zeros(length_list,length_list)

    ## create lists of waypoints along the LCPA-path
    for (count_i,i) in enumerate(start_lexi_list)

        _, pred = dijkstra_shortest_path(Ad, i, targets = start_lexi_list)
        for (count_j,j) in enumerate(start_lexi_list)
            if count_i != count_j
                path_arr = []
                if i!=j
                    # connect site coordinates to cell-midpoint cooridnates of cells that are moved along in LCPA-path
                    path = get_path(j, i, pred)
                    path_coords = ind2coords.(lexi2ind.(path, [TopData]), [TopData])
                    path_arr = zeros(length(path_coords),2); for i in eachindex(path_coords) path_arr[i,:] .= path_coords[i] end
                    path_save[count_i,count_j] = path_arr
                    path_save[count_i,count_j][1,:] = [sites[count_j].coords...]; path_save[count_i,count_j][end,:] = [sites[count_i].coords...]; 
                    length_save[count_i,count_j] = sum([distance(path_arr[i,:], path_arr[i+1,:]) for i=1:length(path_coords)-1]) 
                else 
                    # if start-index == end-index only connect site coordinates and no cell-midpoint coordinates 
                    path_arr = zeros(2,2)
                    path_arr[1,:] = [sites[count_j].coords...]; path_arr[2,:] = [sites[count_i].coords...];
                    path_save[count_i,count_j] = path_arr
                    length_save[count_i,count_j] = distance(path_arr[1,:], path_arr[2,:])
                end
                if PlotOutput GMT.lines!(path_arr, pen=(1.5,:red)) end
            end
        end
    end
    if PlotOutput
    GMT.scatter!(
            float.([sites[i].coords[1] for i in eachindex(sites)]), float.([sites[i].coords[2] for i in eachindex(sites)]),
            fmt=:png, 
            marker=:circle,
            markeredgecolor=0, 
            size=0.1, 
            markerfacecolor=:black, 
            show=true,
            figname=joinpath("src", "OutputFolder", string(Dates.today()), string(TopData.kl))
    )
    end
    printstyled("DONE!", color=:green)
    return path_save, length_save
end

"""
    connect_nodes_LCPA(data::Arch_DataSet, MaterialData::MaterialDataList; k :: Int =100, ℓ::Int =100, PlotOut::Bool=false, conn_func :: ConnectionFunction = Tobler_Hiking_adjusted)

Creates networks for all time periods from the set of archaeological sites and the list of given materials.
If the resulting network from this consists of more than one components, additional edges are added to 
connect all components. Uses the distance determined by the LCPA. Has optimization for if connection functions are symmetric or not. Returns Plots of connected networks with comparison, to what the networks would look like without material presence. Return value is a vector of all created graphs and their adjacency matrices. Contains different version of adjacency matrices for comparison. Two entries each of return vector corresponds to one time period. The first of two entries corresponds to the networks created with material. The values in the first position of that network is a list of sites, the second of an unweighted adjacency matrix, only showing the presence/absence of edges. The third entry is the adjacency matrix, as it was before ensuring that the reusling network is connected. The last entry is the weighted adjacency matrix, i.e., giving the presence and abesence of edges and also their length as calculated by the LCPA.
The second entry of the given time period of the returned vector contains the same information, but without the presence of material.
The parameters 'k' and 'ℓ' are choosen to compress the topographic data into grid cells of desired size. For reference, in the original data each pixel acounts for a 1.5km × 1.5km area. Using the inbuild compression the final grid cells used for the LCPA are of size (k⋅1.5km)×(ℓ⋅1.5km).
# Arguments:
- `data::Arch_DataSet`: List of all archeological sites with their properties
- `MaterialData::MaterialDataList`: List of all materials and their properties
- `k`: Number of topography-pixels which are compressed into one -- width; in preprocessing of LCPA
- `ℓ`: Number of topography-pixels which are compressed into one -- heigth; in preprocessing of LCPA
- `PlotOut::Bool=false`: Parameter to decide wether to plot all least cost paths between all nodes. Might take very ong if enabled.
- `conn_func::ConnectionFunction = Tobler_Hiking_adjusted`: Choice of cost-function in LCPA. Standard is modified Tobler Hiking function.
"""
function connect_nodes_LCPA(data::Arch_DataSet, MaterialData::MaterialDataList; k :: Int =100, ℓ::Int =100, PlotOut::Bool=false, conn_func :: ConnectionFunction = Tobler_Hiking_adjusted, prepared::Union{Dict{String,PreparedTimePeriod},Nothing}=nothing)
    prep_data = isnothing(prepared) ? prepare_network_data(data, MaterialData) : prepared
    time_periods = collect(keys(prep_data))
    graphen = []
    for time_period in time_periods
        prep = prep_data[time_period]
        @show time_period
        ## Create output path
        if conn_func.symmetric
            try
                mkpath(joinpath("src", "OutputFolder", string(Dates.today()), time_period, "symm"))
            catch 
                @info "Possibly overwriting files"
            end
        else
            try
                mkpath(joinpath("src", "OutputFolder", string(Dates.today()), time_period, "non_symm"))
            catch
                @info "Possibly overwriting files"
            end
        end
        
        # preprocessing data for network connections
        material_connectivity_arr = adjacency_data(MaterialData, prep)

        relevant_sites = prep.sites

        paths, dists = archeoLCPA(relevant_sites, k, ℓ, PlotOut=PlotOut, conn_func = conn_func)

        anz_rel = length(relevant_sites)
        adjacency_matrix_daten = spzeros(anz_rel, anz_rel)
        adjacency_matrix_daten_compare = spzeros(anz_rel, anz_rel)
        
        # symmetric case
        if conn_func.symmetric
            for i in 1:anz_rel
                for j in (i+1):anz_rel
                    adjacency_matrix_daten[i,j], adjacency_matrix_daten_compare[i,j] = connect_nodes(i,j,dists, material_connectivity_arr, relevant_sites, symmetric = true)
                    adjacency_matrix_daten[j,i] = adjacency_matrix_daten[i,j]
                    adjacency_matrix_daten_compare[j,i] = adjacency_matrix_daten_compare[i,j]
                end
            end
        else ## asymmetric case
            for i in 1:anz_rel
                for j in (i+1):anz_rel
                    adjacency_matrix_daten[i,j], adjacency_matrix_daten_compare[i,j],adjacency_matrix_daten[j,i], adjacency_matrix_daten_compare[j,i] = connect_nodes(i,j,dists, material_connectivity_arr, relevant_sites, symmetric = false)
                end
            end
        end
        ## copies for comparison of outputs
        non_conn_adj = copy(adjacency_matrix_daten); non_conn_adj_comp = copy(adjacency_matrix_daten_compare)

        ## Saving needs to be exchanged by function 

        string_help = conn_func.symmetric ? "symm" : "non_symm"
        writedlm(joinpath("src", "OutputFolder", string(Dates.today()), time_period, string_help, "non_conn_$(time_period).csv"), adjacency_matrix_daten, ',')
        writedlm(joinpath("src", "OutputFolder", string(Dates.today()), time_period, string_help, "no_mat_non_conn_$(time_period).csv"), adjacency_matrix_daten_compare, ',')

        ### Material edges only
        writedlm(joinpath("src", "OutputFolder", string(Dates.today()), time_period, string_help, "material_diff_$(time_period).csv"), adjacency_matrix_daten-adjacency_matrix_daten_compare, ',')
        
        ## Connect all components by means of shortest distances between clusters
        connect_components!(adjacency_matrix_daten, components(adjacency_matrix_daten)..., dists)
        connect_components!(adjacency_matrix_daten_compare, components(adjacency_matrix_daten_compare)..., dists)
        

        writedlm(joinpath("src", "OutputFolder", string(Dates.today()), time_period, string_help, "$(time_period).csv"), adjacency_matrix_daten, ',')
        writedlm(joinpath("src", "OutputFolder", string(Dates.today()), time_period, string_help, "no_mat_$(time_period).csv"), adjacency_matrix_daten_compare, ',')

        weighted_adjacency_matrix = adjacency_matrix_daten .* dists
         # write adjacency matrix to csv
        writedlm(joinpath("src", "OutputFolder", string(Dates.today()), time_period, string_help, "weight_$(time_period).csv"), weighted_adjacency_matrix, ',')

        ## Push final graph
        push!(graphen, [relevant_sites, adjacency_matrix_daten, non_conn_adj, weighted_adjacency_matrix])

        weighted_adjacency_matrix_comp = adjacency_matrix_daten_compare .* dists
        writedlm(joinpath("src", "OutputFolder", string(Dates.today()), time_period, string_help, "weight_no_mat$(time_period).csv"), weighted_adjacency_matrix_comp, ',')

        ## Push comparison of final graph
        push!(graphen, [relevant_sites, adjacency_matrix_daten_compare, non_conn_adj_comp, weighted_adjacency_matrix_comp])

        long = [entry.coords[1] for entry in relevant_sites]
        lat = [entry.coords[2] for entry in relevant_sites]

        left_border, right_border = extrema([p.coords[1] for p in relevant_sites])
        bottom_border, top_border = extrema([p.coords[2] for p in relevant_sites])

        Δlong = 2; Δlat=1;
        if conn_func.symmetric
            GMT.coast(
            region="$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+2*Δlong)/$(top_border+Δlat)+r", 
            proj=(name=:laea, center=[10,50]), frame=:ag, res=:low, area=500, shore=:thin, borders=1,
            water="deepskyblue3@50", land="palegreen4@50",
        )
        @showprogress for knoten in axes(adjacency_matrix_daten,1)
            for adjacent_knoten in adjacency_matrix_daten[:,knoten].nzind
                if ((adjacent_knoten in non_conn_adj_comp[:,knoten].nzind) && (adjacent_knoten in non_conn_adj[:,knoten].nzind)) && .!isnothing(paths[knoten, adjacent_knoten])
                    GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:black))
                elseif adjacent_knoten in non_conn_adj[:,knoten].nzind && .!isnothing(paths[knoten, adjacent_knoten])
                    GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:red))
                elseif .!isnothing(paths[knoten, adjacent_knoten])
                    GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:cyan))
                end
            end
        end
        GMT.scatter!(
            long, lat, mt=:png, marker=:circle,
            markeredgecolor=0, size=0.1, markerfacecolor=:blue
        )
        GMT.legend!((
            vspace=-0.25,
            header=(text="Archaeological graphs", font=(12,"Times-Roman")),
            hline1=(pen=1, offset=0.5),
            ncolumns=2,
            symbol1=(marker=:circle, size=0.4, dx_left=0.25, fill=:blue, dx_right=0.75, text="Site from $time_period"),
            symbol2=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:black), dx_right=0.75, text="Default movement"),
            vline=(pen=1, offset=0),
            symbol3=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:red), dx_right=0.75, text="Material Connectivity"),
            symbol4=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:cyan), dx_right=0.75, text="Connected graph")
            ),
            pos=(anchor=:TL, width=10, justify=:BL, spacing=1.2),
            clearance=(0.25,0.25), box=(pen=0.5, fill=:azure1),
            figsize=10, show=true,
            figname=joinpath("src", "OutputFolder", string(Dates.today()), "$(time_period)_symm")
            )
        else
            GMT.coast(
                region="$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+2*Δlong)/$(top_border+2*Δlat)+r",
                proj=(name=:laea, center=[10,50]), frame=:ag, res=:low, area=500, shore=:thin, 
                borders=1, water="deepskyblue3@50", land="palegreen4@50",
            )

            @showprogress for knoten in axes(adjacency_matrix_daten,1)
                adjacent_vec = adjacency_matrix_daten[:,knoten].nzind
                ## There ...  
                for adjacent_knoten in adjacent_vec[adjacent_vec.<knoten]
                    if ((adjacent_knoten in non_conn_adj_comp[:,knoten].nzind) && (adjacent_knoten in non_conn_adj[:,knoten].nzind)) && .!isnothing(paths[knoten, adjacent_knoten])
                        GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:black))
                    elseif adjacent_knoten in non_conn_adj[:,knoten].nzind && .!isnothing(paths[knoten, adjacent_knoten])
                        GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:lightcoral))
                    elseif .!isnothing(paths[knoten, adjacent_knoten])
                        GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:cyan))
                    end
                end
                ## and back again...
                for adjacent_knoten in adjacent_vec[adjacent_vec.>knoten]
                    if ((adjacent_knoten in non_conn_adj_comp[:,knoten].nzind) && (adjacent_knoten in non_conn_adj[:,knoten].nzind)) && .!isnothing(paths[knoten, adjacent_knoten])
                        GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:white))
                    elseif adjacent_knoten in non_conn_adj[:,knoten].nzind && .!isnothing(paths[knoten, adjacent_knoten])
                        GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:red))
                    elseif .!isnothing(paths[knoten, adjacent_knoten])
                        GMT.lines!(paths[knoten, adjacent_knoten], pen=(0.5,:darkblue))
                    end
                end
                ### -- A hobbits tale, by
            end
            GMT.scatter!(
                long, lat, 
                fmt=:png, 
                marker=:circle,
                markeredgecolor=0, 
                size=0.1, 
                markerfacecolor=:blue, 
            )
            ### Bilbo Baggins
            GMT.legend!((
                vspace=-0.25,
                header=(text="Archaeological graphs", font=(12,"Times-Roman")),
                hline1=(pen=1, offset=0.5),
                symbol1=(marker=:circle, size=0.4, dx_left=3.5, fill=:blue, dx_right=3.75, text="Site from $time_period"),
                hline2=(pen=1, offset=0.5),
                ncolumns=2,
                vline=(pen=1, offset=0),
                symbol2=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:black), dx_right=0.75, text="Default movement || There"),
                symbol3=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:white), dx_right=0.75, text="Default movement || Back"),
                symbol4=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:lightcoral), dx_right=0.75, text="Material Connectivity || There"),
                symbol5=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:red), dx_right=0.75, text="Material Connectivity || Back"),
                symbol6=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:cyan), dx_right=0.75, text="Connected graph || There"),
                symbol7=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:darkblue), dx_right=0.75, text="Connected graph || Back")
                ),
                pos=(offset=(0.25, -0.25), anchor=:TL, width=10, justify=:BL, spacing=1.2),
                clearance=(0.25,0.25), box=(pen=0.5, fill=:azure1),
                figsize=10, show=true, figname=joinpath("src", "OutputFolder", string(Dates.today()), "$(time_period)_nonsymm")
                )
        end
    end
    return graphen
end