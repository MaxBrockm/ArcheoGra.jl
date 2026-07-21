### Functions to connect sites from consecutive time periods to model the transition between time slices
"""
    create_connections(archdata::Arch_DataSet, ordered_times::Vector{String}; k::Int=16, ℓ::Int=16, conn_func::ConnectionFunction=Tobler_Hiking_mod)

This function creates connections between vertices of different time periods. The connections are intended to support the transitions between time periods, when further analysing the networks. 
Connections are created similar to the networks within one time period, however only connecting existing sites to new sites of the next time period. The order of time periods is entered as a parameter of the function. The connections are made regardeless of materials or groupings of the archaeological data. 

# Arguments:
- `archdata`: Archaeological data, containing data from all entered time periods
- `ordered_times`: List of sorted times, each to consecutive time periods are connected
- `k`, `ℓ`: Parameters for compression of topographic data in LCPA  
"""
function create_connections(archdata::Arch_DataSet, ordered_times::Vector{String}; k::Int=16, ℓ::Int=16, conn_func::ConnectionFunction=Tobler_Hiking_mod, prepared::Union{Dict{String,PreparedTimePeriod},Nothing}=nothing)
    ## Build a quick time_period → sites index (no material data needed here)
    tp_sites = if !isnothing(prepared)
        Dict(tp => p.sites for (tp, p) in prepared)
    else
        _idx = Dict{String, Vector{ArchDataEntry}}()
        for entry in archdata.Data
            push!(get!(Vector{ArchDataEntry}, _idx, entry.time_period), entry)
        end
        _idx
    end

    sites_time1 = nothing
    sites_time2 = tp_sites[ordered_times[1]]

    for times in 1:length(ordered_times)-1
        time_1 = ordered_times[times]
        time_2 = ordered_times[times+1]

        try
            mkpath(joinpath("src", "OutputFolder", string(Dates.today()), "cross_time"))
        catch
            @info "Possibly overwriting files!"
        end
        sites_time1 = sites_time2
        Set_coords_time1 = Set([sites_time1[i].coords for i in eachindex(sites_time1)])

        # Determine which sites need to be connected to cross time period
        sites_time2 = tp_sites[time_2]
        new_sites = [sites_time2[j] for j in eachindex(sites_time2) if sites_time2[j].coords ∉ Set_coords_time1]
        
        long_sites_time1 = [p.coords[1] for p in sites_time1]
        lat_sites_time_1 = [p.coords[2] for p in sites_time1]

        long_new_sites = [p.coords[1] for p in new_sites]
        lat_new_sites = [p.coords[2] for p in new_sites]
        # Setup cost grid
        top_data = tiff_to_matrix(joinpath("src", "DataFolder", "-9_35_32_54.tiff"))

        left_border, right_border = extrema([extrema(long_sites_time1)..., extrema(long_new_sites)...])
        bottom_border, top_border = extrema([extrema(lat_sites_time_1)..., extrema(lat_new_sites)...])

        long_p_px = 41/2460 ## specific to this topo. dataset
        lat_p_px = 19/1140
    
        lx = floor(Int, (left_border+9)/long_p_px)
        ux = ceil(Int, (right_border+9)/long_p_px)
        ly = floor(Int, (bottom_border-35)/lat_p_px)
        uy = ceil(Int, (top_border-35)/lat_p_px)
        
        top_data = top_data[ly-10:uy+10, lx-10:ux+10]
        
        new_lb_corner = (-9+lx*long_p_px, 35+ly*lat_p_px); new_rt_corner = (-9+ux*long_p_px, 35+uy*lat_p_px)
      
        useTopdata = TopographicalData(top_data, size(top_data), new_lb_corner, new_rt_corner, (1,1))
        compressed_topData = compress(useTopdata, k, ℓ)
        
        ind_list = lexi2ind.([coords2lexi(archdata.Data[i].coords, compressed_topData) for i in eachindex(archdata.Data)], [compressed_topData])
        Ad_water = weighted_adjacency_matrix_haversine(compressed_topData, conn_func.f)

        LCPA_Ad = delConnections(Ad_water, lexiMat2Vec(WaterLandBool(compressed_topData, ind_list)))

        # Calculate distances between each one old vertex to all new vertices
        cross_time_distances_there = zeros(length(sites_time1), length(new_sites))
        cross_time_paths_there = Array{Any}(undef, length(sites_time1), length(new_sites))
        for i in eachindex(sites_time1)
            cross_time_paths_there[i,:], cross_time_distances_there[i,:] = one_to_rest_LCPA(compressed_topData, LCPA_Ad, sites_time1[i], new_sites) 
        end

        cross_time_distances_back = zeros(length(new_sites), length(sites_time1))
        cross_time_paths_back = Array{Any}(undef, length(new_sites), length(sites_time1))
        if !conn_func.symmetric
            for i in eachindex(new_sites)
                cross_time_paths_back[i,:], cross_time_distances_back[i,:] = one_to_rest_LCPA(compressed_topData, LCPA_Ad, new_sites[i], sites_time1) 
            end
        end

        cross_time_adjacency_there = spzeros(length(sites_time1), length(new_sites))
        cross_time_adjacency_back = spzeros(length(new_sites), length(sites_time1))
        # Create connections between old and new vertices
        if conn_func.symmetric 
            @showprogress for i in eachindex(sites_time1)
                for j in eachindex(new_sites)
                    cross_time_adjacency_there[i,j] = connect_nodes(i, j, cross_time_distances_there)
                end
            end
            writedlm(joinpath("src", "OutputFolder", string(Dates.today()), "cross_time", "cross_time_adj_$(time_1)$(time_2).csv"), cross_time_adjacency_there, ',')
            @show nnz(cross_time_adjacency_there)/(length(sites_time1)*length(new_sites))
        else 
            @showprogress for i in eachindex(sites_time1)
                for j in eachindex(new_sites)
                    cross_time_adjacency_there[i,j], cross_time_adjacency_back[j,i] = connect_nodes(i, j, cross_time_distances_there,conn_func.symmetric, distance_matrix_back = cross_time_distances_back)
                end
            end
            writedlm(joinpath("src", "OutputFolder", string(Dates.today()), "cross_time", "cross_time_adj_$(time_1)$(time_2)_there.csv"), cross_time_adjacency_there, ',')
            writedlm(joinpath("src", "OutputFolder", string(Dates.today()), "cross_time", "cross_time_adj_$(time_1)$(time_2)_back.csv"), cross_time_adjacency_back, ',')
            @show (nnz(cross_time_adjacency_there)+nnz(cross_time_adjacency_back))/(2*length(sites_time1)*length(new_sites))
        end

        # Plotten
        Δlong = 2; Δlat=1;
        if conn_func.symmetric
            GMT.coast(
            region="$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+2*Δlong)/$(top_border+Δlat)+r", 
            proj=(name=:laea, center=[10,50]), frame=:ag, res=:low, area=500, shore=:thin, borders=1,
            water="deepskyblue3@50", land="palegreen4@50",
        )
        @showprogress for knoten in axes(cross_time_adjacency_there,1)
            for adjacent_knoten in cross_time_adjacency_there[knoten,:].nzind
                if .!isnothing(cross_time_paths_there[knoten, adjacent_knoten])
                    GMT.lines!(cross_time_paths_there[knoten, adjacent_knoten], pen=(0.5,:black))
                end
            end
        end
        GMT.scatter!(
            long_sites_time1, lat_sites_time_1, mt=:png, marker=:circle,
            markeredgecolor=0, size=0.1, markerfacecolor=:blue
        )
        GMT.scatter!(long_new_sites, lat_new_sites, mt=:png, marker=:circle,
        markeredgecolor=0, size=0.1, markerfacecolor=:red)
        GMT.legend!((
            vspace=-0.25,
            header=(text="Cross-time-connections", font=(12,"Times-Roman")),
            hline1=(pen=1, offset=0.5),
            ncolumns=2,
            symbol1=(marker=:circle, size=0.4, dx_left=0.25, fill=:blue, dx_right=0.75, text="Site from $time_1"),
            symbol2=(marker=:circle, size=0.4, dx_left=0.25, fill=:red, dx_right=0.75, text="Site from $time_2"),
            vline=(pen=1, offset=0),
            symbol3=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:black), dx_right=0.75, text="Connections")
            ),
            pos=(anchor=:TL, width=10, justify=:BL, spacing=1.2),
            clearance=(0.25,0.25), box=(pen=0.5, fill=:azure1),
            figsize=10, show=true,
            figname=joinpath("src", "OutputFolder", string(Dates.today()), "cross_time", "$(time_1)$(time_2)_non_symm")
            )
        else
            GMT.coast(
            region="$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+2*Δlong)/$(top_border+Δlat)+r", 
            proj=(name=:laea, center=[10,50]), frame=:ag, res=:low, area=500, shore=:thin, borders=1,
            water="deepskyblue3@50", land="palegreen4@50",
        )
        @showprogress for knoten in axes(cross_time_adjacency_there,1)
            for adjacent_knoten in cross_time_adjacency_there[knoten,:].nzind
                if .!isnothing(cross_time_paths_there[knoten, adjacent_knoten])
                    GMT.lines!(cross_time_paths_there[knoten, adjacent_knoten], pen=(0.5,:black))
                end
            end
        end
        @showprogress for knoten in axes(cross_time_adjacency_back,1)
            for adjacent_knoten in cross_time_adjacency_back[knoten,:].nzind
                if .!isnothing(cross_time_paths_back[knoten, adjacent_knoten])
                    GMT.lines!(cross_time_paths_back[knoten, adjacent_knoten], pen=(0.5,:white))
                end
            end
        end
        GMT.scatter!(
            long_sites_time1, lat_sites_time_1, mt=:png, marker=:circle,
            markeredgecolor=0, size=0.1, markerfacecolor=:blue
        )
        GMT.scatter!(long_new_sites, lat_new_sites, mt=:png, marker=:circle,
        markeredgecolor=0, size=0.1, markerfacecolor=:red)
        GMT.legend!((
            vspace=-0.25,
            header=(text="Cross-time-connections", font=(12,"Times-Roman")),
            hline1=(pen=1, offset=0.5),
            ncolumns=2,
            symbol1=(marker=:circle, size=0.4, dx_left=0.25, fill=:blue, dx_right=0.75, text="Site from $time_1"),
            symbol2=(marker=:circle, size=0.4, dx_left=0.25, fill=:red, dx_right=0.75, text="Site from $time_2"),
            vline=(pen=1, offset=0),
            symbol3=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:black), dx_right=0.75, text="Connections There"),
            symbol4=(marker="-", size=0.6,  dx_left=0.25, pen=(0.5,:white), dx_right=0.75, text="Connections Back")
            ),
            pos=(anchor=:TL, width=10, justify=:BL, spacing=1.2),
            clearance=(0.25,0.25), box=(pen=0.5, fill=:azure1),
            figsize=10, show=true,
            figname=joinpath("src", "OutputFolder", string(Dates.today()), "cross_time", "$(time_1)$(time_2)_non_symm")
            )
        end
    end
end


"""
    one_to_rest_LCPA(TopData::TopographicalData, Ad::AbstractMatrix, one_site::ArchDataEntry, rest_sites::Vector{ArchDataEntry})

Calculates the least cost paths from one specific site to all other provided sites.

Required in cross-time connections, when vertices of the same time period are not 
intended to be connected among themselves, but only to other time periods.

# Arguments
- `TopData::TopographicalData`: The topographical map data.
- `Ad::AbstractMatrix`: The adjacency matrix of the graph.
- `one_site::ArchDataEntry`: The starting site (must contain a `.coords` tuple).
- `rest_sites::Vector{ArchDataEntry}`: A list of target sites.

# Output
- `paths_save::Vector{Matrix{Float64}}`: A vector of coordinate matrices, each representing a path.
- `lengths_save::Vector{Float64}`: A vector containing the total geographic length of each path.
"""
function one_to_rest_LCPA(TopData::TopographicalData, Ad::AbstractMatrix, one_site::ArchDataEntry, rest_sites::Vector{ArchDataEntry})
    rest_lexi_list = [coords2lexi(site.coords, TopData) for site in rest_sites]
    one_lexi = coords2lexi(one_site.coords, TopData)
    
    # Type-stable initialization
    paths_save = Vector{Matrix{Float64}}(undef, length(rest_lexi_list))
    lengths_save = zeros(length(rest_lexi_list))
    
    # Use the new, highly optimized Min-Heap Dijkstra with early stopping
    _, pred = dijkstra_shortest_path(Ad, one_lexi, targets=rest_lexi_list)
    
    for j in eachindex(rest_lexi_list)
        if one_lexi != rest_lexi_list[j]
            # Get path nodes and convert to coordinates
            path_nodes = get_path(rest_lexi_list[j], one_lexi, pred)
            path_coords = [ind2coords(lexi2ind(node, TopData), TopData) for node in path_nodes]
            
            # Initialize Matrix for this path
            paths_save[j] = zeros(length(path_coords), 2)
            
            # Fill intermediate coordinates
            for i in 2:(length(path_coords)-1)
                paths_save[j][i, :] .= path_coords[i]
            end
            
            # Explicitly set start and end to the exact ArchDataEntry coordinates
            paths_save[j][1, :] .= [rest_sites[j].coords...]
            paths_save[j][end, :] .= [one_site.coords...]
            
            # Calculate path geographic length
            lengths_save[j] = sum(distance((paths_save[j][i, 1], paths_save[j][i, 2]), 
                                           (paths_save[j][i+1, 1], paths_save[j][i+1, 2])) 
                                  for i in 1:(length(path_coords)-1))
        else
            # Edge case: Start and target are the same
            paths_save[j] = zeros(2, 2)
            paths_save[j][1, :] .= [rest_sites[j].coords...]
            paths_save[j][2, :] .= [one_site.coords...]
            lengths_save[j] = distance((paths_save[j][1, 1], paths_save[j][1, 2]), 
                                       (paths_save[j][2, 1], paths_save[j][2, 2]))
        end
    end
    
    return paths_save, lengths_save
end