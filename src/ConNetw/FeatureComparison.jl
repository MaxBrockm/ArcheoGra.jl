"""
    extract_degree_distribution(Ad::SparseMatrixCSC)

Extracts the standard degree and the two-hop degree distribution for all vertices in the graph.

# Arguments
- `Ad::SparseMatrixCSC`: The adjacency matrix of the graph.

# Output
- `Tuple{Vector{Int}, Vector{Float64}}`: A tuple containing the degree vector and the two-hop degree vector.
"""
function extract_degree_distribution(Ad::SparseMatrixCSC)
    degree_vector = vec(count(>(0), Ad, dims=2))
    two_hop_degree_vector = [1 / degree_vector[i] * sum(degree_vector[Ad[:, i].nzind]) for i in axes(Ad, 1)]
    return degree_vector, two_hop_degree_vector
end 

"""
    extract_clustering_coefficient(Ad::SparseMatrixCSC)

Calculates the local clustering coefficient and the two-hop clustering coefficient for each vertex.

# Arguments
- `Ad::SparseMatrixCSC`: The adjacency matrix of the graph.

# Output
- `Tuple{Vector{Float64}, Vector{Float64}}`: A tuple containing the clustering coefficients and two-hop clustering coefficients.
"""
function extract_clustering_coefficient(Ad::SparseMatrixCSC)
    num_verts = size(Ad,1)
    coeff_vector = zeros(num_verts)
    two_hop_coeff_vector = zeros(num_verts)
    for i in 1:num_verts
        neigh = Ad[:,i].nzind
        # edges counted twice in adjacency matrix so factor 2 is suppressed
        coeff_vector[i] = length(neigh) == 1 ? 0 : length(Ad[neigh,neigh].nzval)/(length(neigh)*(length(neigh)-1)) 
    end
    for i in 1:num_verts
        neigh = Ad[:,i].nzind
        two_hop_coeff_vector[i] = sum(coeff_vector[neigh]) / length(neigh)
    end
    return coeff_vector, two_hop_coeff_vector
end

"""
    extract_egonet_features(Ad::SparseMatrixCSC)

Extracts features related to the ego-network of each vertex, specifically the number of outgoing edges 
and the number of unique neighbors adjacent to the ego-network.

# Arguments
- `Ad::SparseMatrixCSC`: The adjacency matrix of the graph.

# Output
- `Tuple{Vector{Int}, Vector{Int}}`: Outgoing edges from the ego-network and the number of unique neighbors.
"""
function extract_egonet_features(Ad::SparseMatrixCSC)
    anz_verts = size(Ad,1)
    ego_net_out_edges = zeros(anz_verts)
    ego_net_neigh = zeros(anz_verts)
    for i in 1:anz_verts
        ego_net = vcat(i, Ad[:,i].nzind)
        ego_net_incicator = zeros(Bool,anz_verts); ego_net_incicator[ego_net] .= true 
        ego_net_out_edges[i] = length(Ad[ego_net_incicator, .!ego_net_incicator].nzval)
        ego_net_neigh[i] = length(unique(Ad[.!ego_net_incicator,ego_net_incicator].rowval))
    end
    return ego_net_out_edges, ego_net_neigh
end


"""
    extract_path_features(Ad::SparseMatrixCSC, max_distance::Number)

Calculates the betweenness centrality, closeness centrality, and normalized diameter of the graph 
using an optimized Dijkstra's algorithm with a binary min-heap.

# Arguments
- `Ad::SparseMatrixCSC`: The weighted adjacency matrix.
- `max_distance::Number`: The maximum possible distance for normalization.

# Output
- `Tuple{Vector{Float64}, Vector{Float64}, Float64}`: Betweenness centrality, closeness centrality, and normalized diameter.
"""
function extract_path_features(Ad::SparseMatrixCSC, max_distance::Number)
    anz_verts = size(Ad, 1)
    betweenness_centrality_vector = zeros(anz_verts) 
    closeness_centrality_vector = zeros(anz_verts)
    diameter  = 0

    for s in 1:anz_verts ## calculates δₛ(v) for all v simultaneously
        
        stack = Int[] ## visited notes
        P = [Int[] for _ in 1:anz_verts] ## predecessors
        sigma = zeros(Int, anz_verts) ## number of shortest paths from s to other vertices 
        dist = fill(Inf, anz_verts) ## distances from s to other vertices 
        visited = falses(anz_verts) ## list of unsivited nodes

        sigma[s] = 1
        dist[s] = 0

        ### Dijkstra 
        queue = collect(1:anz_verts)  

        while !isempty(queue)
            # Find node u with minimal dist[u] in queue; starts with s
            u_idx = argmin(dist[queue])
            u = queue[u_idx]
            deleteat!(queue, u_idx)

            push!(stack, u)

            for v in Ad[:,u].nzind ## go through every node that is connected to u
                alt = dist[u] + Ad[u, v] ## cost to go through u to v 
                if alt < dist[v]
                    dist[v] = alt
                    sigma[v] = sigma[u]
                    P[v] = [u]
                elseif alt == dist[v] ## keep track of multiple shortest paths
                    sigma[v] += sigma[u]
                    push!(P[v], u)
                end
            end
        
        end

        # Update δₛ(v) for all v; needs to go in reverse from most distante end points of paths from s to closest neighbours of s 
        delta = zeros(anz_verts)
        while !isempty(stack)
            v = pop!(stack)
            for w in P[v]
                delta[w] += (sigma[w] / sigma[v]) * (1 + delta[v])
            end
            if v != s
                betweenness_centrality_vector[v] += delta[v]
            end
        end

        closeness_centrality_vector[s] = (anz_verts - 1)/sum(dist) ## uses already calculated shorest paths 
        diameter  = maximum([diameter; dist])
    end

    # normalizing for undirected graphs 
    betweenness_centrality_vector ./= 2
    betweenness_centrality_vector .*= 2/((anz_verts-1)*(anz_verts-2))

    return betweenness_centrality_vector, closeness_centrality_vector, diameter/max_distance
end


"""
    extract_density(Ad::SparseMatrixCSC)

Calculates the edge density of the given graph.

# Arguments
- `Ad::SparseMatrixCSC`: The adjacency matrix.

# Output
- `Float64`: The computed graph density.
"""
function extract_density(Ad::SparseMatrixCSC)
    anz_verts = size(Ad,1)
    return length(Ad.nzval)/(anz_verts*(anz_verts-1))
end


"""
    _heappush!(heap::Vector{Tuple{Float64, Int}}, item::Tuple{Float64, Int})

Pushes an item onto the binary min-heap while maintaining the heap invariant.
"""
function _heappush!(heap::Vector{Tuple{Float64, Int}}, item::Tuple{Float64, Int})
    push!(heap, item)
    # sift up
    i = length(heap)
    @inbounds while i > 1
        p = i >> 1
        heap[i][1] < heap[p][1] || break
        heap[i], heap[p] = heap[p], heap[i]
        i = p
    end
end

"""
    _heappop!(heap::Vector{Tuple{Float64, Int}})

Pops the minimum element from a binary min-heap of `(cost, node)` tuples.
"""
function _heappop!(heap::Vector{Tuple{Float64, Int}})
    n = length(heap)
    heap[1], heap[n] = heap[n], heap[1]
    result = pop!(heap)
    n -= 1
    # sift down
    i = 1
    @inbounds while true
        smallest = i
        l = 2i; r = 2i + 1
        if l <= n && heap[l][1] < heap[smallest][1]; smallest = l; end
        if r <= n && heap[r][1] < heap[smallest][1]; smallest = r; end
        smallest == i && break
        heap[i], heap[smallest] = heap[smallest], heap[i]
        i = smallest
    end
    return result
end


"""
    feature_aggregation(feature_vector::Vector{Float64})

Aggregates a feature vector into four statistical categories: mean, median, coefficient of variation, and skewness. 
This summarizes vectors of different lengths into easily comparable metrics.

# Arguments
- `feature_vector::Vector{Float64}`: The raw vector of features to aggregate.

# Output
- `Vector{Float64}`: A 4-element vector containing the aggregated statistics.
"""
function feature_aggregation(feature_vector::Vector{Float64})
    aggregated_vector = zeros(4)
    ## mean
    aggregated_vector[1] = mean(feature_vector)
    ## median
    aggregated_vector[2] = median(feature_vector)

    σ = std(feature_vector, mean=aggregated_vector[1])
    aggregated_vector[3] = aggregated_vector[1] == 0 ? 0.0 : σ / aggregated_vector[1]
    aggregated_vector[4] = σ == 0 ? 0.0 : sum(((feature_vector .- aggregated_vector[1]) ./ σ).^3) / length(feature_vector)
    return aggregated_vector
end

"""
    network_features(Ad::SparseMatrixCSC, max_distance::Number)

Extracts and aggregates all primary network features from the given adjacency matrix.

# Arguments
- `Ad::SparseMatrixCSC`: The weighted adjacency matrix.
- `max_distance::Number`: Maximum geographic distance used for normalization.

# Output
- `Tuple{Matrix{Float64}, Float64, Float64}`: Aggregated feature matrix (8x4), normalized diameter, and graph density.
"""
function network_features(Ad::SparseMatrixCSC, max_distance; print_results::Bool = false)
    
    degree_dist, two_hop_degree_dist = extract_degree_distribution(Ad)
    coeff_vector, two_hop_coeff_vector = extract_clustering_coefficient(Ad)
    ego_net_out_edges, ego_net_neigh = extract_egonet_features(Ad)
    betweenness_centrality_vector, closeness_centrality_vector, diameter = extract_path_features(Ad, max_distance)
    density = extract_density(Ad)

    aggregated_features_mat = zeros(8,4)
    for (i,feature) in enumerate([degree_dist, two_hop_degree_dist, coeff_vector, two_hop_coeff_vector, ego_net_out_edges, ego_net_neigh, betweenness_centrality_vector, closeness_centrality_vector])
        aggregated_features_mat[i,:] = feature_aggregation(feature)
    end
    if print_results
        colnames = ["mean", "median", "Coefficient of Variation", "skewness"]
        rownames = ["degree", "two-hop-degree", "clustering coeff.", "two-hop-clustering coeff.", "#outgoing edges of ego-network", "#neighbours of ego-network", "betweenness centrality", "closeness centrality"]
        pretty_table(aggregated_features_mat, header=colnames, row_labels=rownames, backend=Val(:latex))
        println("Normalized diameter of graph: $diameter.\\\\")
        println("Density of graph: $density.\\\\")
    end
    return aggregated_features_mat, diameter, density
end


"""
    mahalanobis_pinv(x::AbstractVector, μ::AbstractVector, Σ::AbstractMatrix; rtol::Float64=1e-10)

Calculates the squared Mahalanobis distance using the pseudo-inverse via eigendecomposition, ensuring numerical stability.

# Arguments
- `x::AbstractVector`: The target vector.
- `μ::AbstractVector`: The mean vector.
- `Σ::AbstractMatrix`: The covariance matrix.

# Keywords
- `rtol::Float64`: Relative tolerance to filter out small eigenvalues.

# Output
- `Float64`: The squared Mahalanobis distance.
"""
function mahalanobis_pinv(x, μ, Σ; rtol=1e-10)
    d = x .- μ
    F = eigen(Hermitian(Σ))
    keep = F.values .> rtol * maximum(F.values, init=1.0)
    Q = F.vectors[:, keep]; Λ = F.values[keep]
    y = Q' * d
    return sum((y.^2) ./ Λ)                   # d²
end

"""
    mahalanobis_distance(x::AbstractVector, μ::AbstractVector, cov_mat::AbstractArray)

Wrapper function to calculate the Mahalanobis distance robustly.
"""
@inline function mahalanobis_distance(x::AbstractVector, μ::AbstractVector, cov_mat::AbstractArray) 
    return mahalanobis_pinv(x, μ, cov_mat)
end

"""
    standardized_distance_squared(x::AbstractVector, μ::AbstractVector, σ::AbstractVector)

Calculates the squared standardized distance (Z-score squared) between vectors.
"""
function standardized_distance_squared(x::AbstractVector, μ::AbstractVector, σ::AbstractVector)
    d = 0
    for i in eachindex(x)
        d += (σ[i]>0 ? ((x[i]-μ[i])/σ[i])^2 : 0)
    end
    return d
end


"""
    empirical_quantile(sorted::AbstractVector{<:Real}, p::Real)

Calculates the empirical `p`-quantile of a sorted vector with linear interpolation.
"""
function empirical_quantile(sorted::AbstractVector{<:Real}, p::Real)
    @assert 0.0 <= p <= 1.0
    n = length(sorted)
    n == 1 && return sorted[1]
    rank = p * (n - 1) + 1
    i = floor(Int, rank); f = rank - i
    return i >= n ? sorted[end] : (1 - f) * sorted[i] + f * sorted[i + 1]
end

"""
    empirical_cdf_interpolated(sorted::AbstractVector{<:Real}, x::Real)

Evaluates the empirical cumulative distribution function (CDF) with linear interpolation.
"""
function empirical_cdf_interpolated(sorted::AbstractVector{<:Real}, x::Real)
    s = (sorted[1] <= sorted[end]) ? sorted : @view(sorted[end:-1:1])
    n = length(sorted)
    if x <= s[1];  return 1.0; end
    if x >= s[end]; return 0.0; end
    i = searchsortedlast(s, x)       # 1 ≤ i < n, und sorted[i] ≤ x < sorted[i+1]
    i >= n && return 0.0             # safety guard (e.g. x = NaN/Inf bypasses the >= check above)
    a = s[i]; b = s[i+1]
    if a == b
        return 1-(i - 1) / (n - 1)          
    else
        t = (x - a) / (b - a)             
        return 1-(i - 1 + t) / (n - 1)
    end
end


"""
    create_graph_statistic_Haversine(toml_file::String; num_graphs::Int = 1000, double_distance::Bool = false, output_dir::String = pwd(), overwrite::Bool = false)

Generates statistical features over a large collection of randomized graphs to establish a baseline for the given parameters.
Calculates covariance matrices and Mahalanobis distances to allow testing against standard random network expectations.
Statistics can be created via multiple threads to accelerate the process. Uses all threads available; however the number of threads needs to be preset via the `JULIA_NUM_THREADS` environment variable by the user.

# Arguments
- `toml_file::String`: Path to the configuration TOML.

# Keywords
- `num_graphs::Int`: Number of simulated graphs to generate for the distribution. Defaults to 1000.
- `double_distance::Bool`: Controls the distance penalty calculation. Defaults to false.
- `output_dir::String`: Target export folder. Defaults to the active user working location.
- `overwrite::Bool`: Permits file deletion if directories exist. Defaults to false.

# Output
- `Nothing`: Persists statistics, covariances, and distances as text files to the local file system.
"""
function create_graph_statistic_Haversine(toml_file::String; num_graphs::Int = 1000, double_distance::Bool = false, output_dir::String = pwd(), overwrite::Bool = false)
    data_path = isfile(toml_file) ? toml_file : joinpath(@__DIR__, "src", "DataFolder", toml_file)

    ArData, MaData = read_data_TOML(data_path)
    SoData = nothing; incl_sources = false; src_string = ""
     
    try 
        SoData = read_source_data(data_path, MaData)
        incl_sources = true
        src_string = "source"
    catch
        @info "Failed to read in Source Data. Will continue creating graphs using direct distances without source-nodes"
    end

    ## Pre-compute indexed data once
    prep_data = prepare_network_data(ArData, MaData)
    time_periods = unique(entry.time_period for entry in ArData.Data)

    out_base = joinpath(output_dir, "GraphStatistics", "Haversine")
    if !isdir(out_base)
        mkpath(out_base)
    end

    colnames = ["mean", "median", "Coefficient of Variation", "skewness"]
    rownames = ["degree", "two-hop-degree", "clustering coeff.", "two-hop-clustering coeff.", "#outgoing edges of ego-network", "#neighbours of ego-network", "betweenness centrality", "closeness centrality"]
    num_features = length(colnames)*length(rownames)
    
    feature_dict = Dict{String, AbstractArray}()
    max_distance_dict = Dict{String,Number}()
    distances_dict = Dict{String, AbstractArray}()

    for time_period in time_periods
        target_file = joinpath(out_base, "$(time_period)$(src_string).txt")

        if isfile(target_file)
            if !overwrite
                error("Statistic for the graph from $time_period already exists at $target_file. Set `overwrite=true` to proceed.")
            else
                @warn "Overwriting existing statistical data for $time_period."
            end
        elseif isfile(joinpath(@__DIR__, "src", "OutputFolder", "GraphStatistics", "Haversine", "$(time_period)$(src_string).txt"))
            error("Statistic for the graph from $time_period already exists in the default output folder. Set `overwrite=true` to proceed.")
        end

        feature_dict[time_period * src_string] = zeros(num_features+2, num_graphs)
        feature_dict[time_period * src_string * "compare"] = zeros(num_features+2, num_graphs)

        distances_dict[time_period * src_string] = zeros(num_features+2, num_graphs)
        distances_dict[time_period * src_string * "compare"] = zeros(num_features+2, num_graphs)

        max_distance_dict[time_period * src_string] = maximum(distance_allocation(prep_data[time_period].sites))
    end

    cp(data_path, joinpath(out_base, "Setup_file.toml"), force=true)
    open(joinpath(out_base, "Setup_file.toml"), "a") do io
        write(io, "\nmaximum_distance = $maximum_distance\nmaterial_multiplier=$material_multiplier")
    end

    rng = MersenneTwister(Threads.threadid())
    Random.seed!(rng)

    @info "Initiating parallel graph generation sequence..."
    Threads.@threads for j in 1:num_graphs
        
        if j % 50 == 0
            @debug "Processed graph generation milestone: $j / $num_graphs"
        end

        graphen = (incl_sources ? Haversine_connect_nodes_sources(ArData, MaData, SoData, plot_results=false, verbose = false, double_distance=double_distance, prepared=prep_data) : Haversine_connect_nodes(ArData, MaData, plot_results = false, verbose = false, prepared=prep_data))
        for i in 1:2:length(graphen)
            graph = graphen[i]
            comparison_graph = graphen[i+1]
            time_period = graph[1][1].time_period
            Ad_weighted = graph[4]
            comparison_Ad_weighted = comparison_graph[4]
            features, diameter, density = network_features(Ad_weighted, max_distance_dict[time_period*src_string], print_results=false)
            comp_features, comp_diameter, comp_density = network_features(comparison_Ad_weighted, max_distance_dict[time_period*src_string], print_results=false)
            feature_dict[time_period * src_string][:,j] = [reshape(features', num_features); diameter; density]
            feature_dict[time_period * src_string * "compare"][:,j] = [reshape(comp_features', num_features); comp_diameter; comp_density]
        end

        # Second run for distance calibration calculations
        graphen = (incl_sources ? Haversine_connect_nodes_sources(ArData, MaData, SoData, plot_results=false, verbose = false, double_distance=double_distance, prepared=prep_data) : Haversine_connect_nodes(ArData, MaData, plot_results = false, verbose = false, prepared=prep_data))
        
        for i in 1:2:length(graphen)
            graph = graphen[i]
            comparison_graph = graphen[i+1]
            time_period = graph[1][1].time_period
            Ad_weighted = graph[4]
            comparison_Ad_weighted = comparison_graph[4]
            features, diameter, density = network_features(Ad_weighted, max_distance_dict[time_period*src_string], print_results=false)
            comp_features, comp_diameter, comp_density = network_features(comparison_Ad_weighted, max_distance_dict[time_period*src_string], print_results=false)
            distances_dict[time_period * src_string][:,j] = [reshape(features', num_features); diameter; density]
            distances_dict[time_period * src_string * "compare"][:,j] = [reshape(comp_features', num_features); comp_diameter; comp_density]
        end
    end

    @info "Graph generation complete. Exporting statistics to $out_base."

    for (time_period, feature_mat) in feature_dict
        results = []
        
        mean_vec = vec(mean(feature_mat, dims=2))
        std_vec = vec(std(feature_mat, dims=2))

        cov_mat = cov(feature_mat') 
        
        num_aggregates = length(colnames)
        for (i,feature) in enumerate(rownames)
            for (j,aggre) in enumerate(colnames)
                push!(results, (feature*"-"*aggre, (mean_vec[(i-1)*num_aggregates+j], std_vec[(i-1)*num_aggregates+j])))
            end
        end
        push!(results, ("normalized diameter", (mean_vec[num_aggregates*length(rownames)+1], std_vec[num_aggregates*length(rownames)+1])))
        push!(results, ("density", (mean_vec[num_aggregates*length(rownames)+2], std_vec[num_aggregates*length(rownames)+2])))
        distances = [mahalanobis_distance(distances_dict[time_period][:,j], vec(mean_vec), cov_mat) for j in 1:num_graphs]

        open(joinpath(out_base, "$(time_period).txt"), "w") do io
            println(io, "Key,mean,std")
            for (key, value) in results
                println(io, "$key,$(value[1]),$(value[2])")
            end
        end
        writedlm(joinpath(out_base, "$(time_period)_Covariance.txt"), cov_mat)
        writedlm(joinpath(out_base, "$(time_period)_distances.txt"), distances)
    end
end

"""
    is_average(Ad_weighted::SparseMatrixCSC, max_distance::Number, connection_method::String, time_period::String; incl_sources::Bool=false, comparison::Bool=false, p::Real=0.99, output_dir::String=pwd())

Checks if a given network is statistically typical relative to pre-computed network samples by evaluating its Mahalanobis distance.

# Arguments
- `Ad_weighted::SparseMatrixCSC`: Target network to analyze.
- `max_distance::Number`: Normalization baseline for spatial metrics.
- `connection_method::String`: Identifier for the directory storing model weights.
- `time_period::String`: Historical epoch.

# Keywords
- `incl_sources::Bool`: References the model utilizing physical source data.
- `comparison::Bool`: Toggles evaluation against the comparative baseline.
- `p::Real`: Acceptance threshold mapped against empirical CDF. Defaults to 0.99.
- `output_dir::String`: Reference to the base output location where statistics reside.

# Output
- `Tuple{Bool, Float64, Float64, Float64}`: Indicator boolean, calculated probability, raw Mahalanobis distance, and the threshold limit.
"""
function is_average(Ad_weighted::SparseMatrixCSC, max_distance::Number, connection_method::String, time_period::String; incl_sources::Bool=false, comparison::Bool=false, p::Real=0.99, output_dir::String=pwd())
    source_string = (incl_sources ? "source" : "")
    comp_string   = (comparison ? "compare" : "")
    base = joinpath(output_dir, "GraphStatistics", connection_method, "$(time_period)$(source_string)")
    
    stats_file = base * "$(comp_string).txt"
    cov_file   = base * "$(comp_string)_Covariance.txt"
    dists_file = base * "$(comp_string)_distances.txt"

    if !isfile(stats_file)
        stats_file = joinpath(dirname(@__DIR__), "OutputFolder", "GraphStatistics", connection_method, "$(time_period)$(source_string)$(comp_string).txt")
         if !isfile(stats_file)
            error("Statistical data unavailable for $(connection_method). Execute baseline analysis first.")
         end
    end

    if !isfile(cov_file)
        cov_file = joinpath(dirname(@__DIR__), "OutputFolder", "GraphStatistics", connection_method, "$(time_period)$(source_string)$(comp_string)_Covariance.txt")
         if !isfile(cov_file)
            error("Statistical data unavailable for $(connection_method). Execute baseline analysis first.")
         end
    end

    if !isfile(dists_file)
        dists_file = joinpath(dirname(@__DIR__), "OutputFolder", "GraphStatistics", connection_method, "$(time_period)$(source_string)$(comp_string)_distances.txt")
         if !isfile(dists_file)
            error("Statistical data unavailable for $(connection_method). Execute baseline analysis first.")
         end
    end

    data_arr = readdlm(stats_file, ',', Any, '\n')[2:end, :]              
    μ        = data_arr[:, 2]                                              
    Σ        = readdlm(cov_file)
    dists    = vec(readdlm(dists_file))                                    

    feats, diam, dens = network_features(Ad_weighted, max_distance; print_results=false)
    x  = [reshape(feats', prod(size(feats))); diam; dens]
    md = mahalanobis_distance(x, μ, Σ)                                     

    s = sort(dists)                                                                    
    p_hat = empirical_cdf_interpolated(s, md)                              
    thr   = empirical_quantile(s, p)                                      

    ok = p_hat >= p

    status = ok ? "ACCEPTED" : "REJECTED"
    @info "Mahalanobis \$d^2\$=$(md); P(Acceptance) ≈ $(round(p_hat*100; digits=2))%; Threshold at p=$(Int(round(p*100)))%: $(thr) → $(status)"
    
    return ok, p_hat, md, thr
end

"""
    graph_sensitivity(toml_file::String; num_graphs::Int = 1000, double_distance::Bool = false, output_dir::String = pwd(), overwrite::Bool = false)

Simulates network variations to ascertain robustness of physical connections over multiple statistical passes.
Statistics can be created via multiple threads to accelerate the process. Uses all threads available; however the number of threads needs to be preset via the `JULIA_NUM_THREADS` environment variable by the user.
# Arguments
- `toml_file::String`: Base setup parameters.

# Keywords
- `num_graphs::Int`: Number of simulated cycles.
- `double_distance::Bool`: Connection probability modifier.
- `output_dir::String`: Root path for visual rendering objects.
- `overwrite::Bool`: Grants IO operation overwrites.

# Output
- `Tuple`: Returns aggregated sparse matrices detailing interaction frequencies across epochs.
"""
function graph_sensitivity(toml_file::String; num_graphs::Int = 1000, double_distance::Bool = false, output_dir::String = pwd(), overwrite::Bool = false)
    data_path = isfile(toml_file) ? toml_file : joinpath(@__DIR__, "src", "DataFolder", toml_file)

    ArData, MaData = read_data_TOML(data_path)
    SoData = nothing; incl_sources = false; src_string = ""
     
    try 
        SoData = read_source_data(data_path, MaData)
        incl_sources = true
        src_string = "source"
    catch
        @info "Failed to read in Source Data. Will continue creating graphs using direct distances without source-nodes"
    end

    out_base = joinpath(output_dir, "GraphStatistics", "Haversine")
    if isdir(out_base) && !overwrite
        error("Output directory $out_base already exists. Set `overwrite=true` to proceed.")
    end
    mkpath(out_base)

    ## Pre-compute indexed data once
    prep_data = prepare_network_data(ArData, MaData)
    time_periods = unique(entry.time_period for entry in ArData.Data)

    cp(data_path, joinpath(out_base, "Setup_file.toml"), force=true)
    open(joinpath(out_base, "Setup_file.toml"), "a") do io
        write(io, "\nmaximum_distance = $maximum_distance\nmaterial_multiplier=$material_multiplier" )
    end

    rng = MersenneTwister(Threads.threadid())
    Random.seed!(rng)

    graphen = (incl_sources ? Haversine_connect_nodes_sources(ArData, MaData, SoData, plot_results=false, double_distance=double_distance, prepared=prep_data) : Haversine_connect_nodes(ArData, MaData, plot_results = false, prepared=prep_data))

    Ad_total_Magda = nothing
    Ad_total_oLBK = nothing
    Ad_total_fullLBK = nothing
    Ad_total_Magda_new = nothing
    
    for i in 1:2:length(graphen)
        graph = graphen[i]
        comparison_graph = graphen[i+1]
        time_period = graph[1][1].time_period
        Ad = graph[2]
        if time_period=="Magdafinal"
            Ad_total_Magda = zeros(size(Ad))
        elseif time_period=="oLBK"
            Ad_total_oLBK = zeros(size(Ad))
        elseif time_period=="fullLBK"
            Ad_total_fullLBK = zeros(size(Ad))
        elseif time_period=="Magdanew"
            Ad_total_Magda_new = zeros(size(Ad))
        else
            @error "unknown Time period"
        end
    end

    @info "Initiating parallel sensitivity simulation mapping..."
    Threads.@threads for j in 1:num_graphs
        
        if j % 50 == 0 println(j) end

        graphen = (incl_sources ? Haversine_connect_nodes_sources(ArData, MaData, SoData, plot_results=false, double_distance=double_distance, prepared=prep_data) : Haversine_connect_nodes(ArData, MaData, plot_results = false, prepared=prep_data))
        for i in 1:2:length(graphen)
            graph = graphen[i]
            comparison_graph = graphen[i+1]
            time_period = graph[1][1].time_period
            Ad = graph[2]
            if time_period=="Magdafinal"
                Ad_total_Magda += Ad
            elseif time_period=="oLBK"
                Ad_total_oLBK += Ad
            elseif time_period=="fullLBK"
                Ad_total_fullLBK += Ad
            elseif time_period=="Magdanew"
                Ad_total_Magda_new += Ad
            else
                @error "unknown Time period"
            end
        end
    end

    if !isnothing(Ad_total_Magda)
        SourceData = SoData["Magdafinal"][1] 
        relevant_sites = ArData.Data[[entry.time_period=="Magdafinal" for entry in ArData.Data]]
        relevant_ids = [site.grid_id for site in relevant_sites] 
        source_ids = [source.grid_id for source in SourceData.Data]
        ids = [relevant_ids; source_ids]
        Ad_total_Magda ./= num_graphs
        display(plot_network_robustness(Ad_total_Magda, name="Magdafinal", output_dir=out_base))
        writedlm(joinpath(out_base, "Magdafinal_adjacency_probabilities.csv"), [[""] ids'; ids Ad_total_Magda], ',')
    end
    if !isnothing(Ad_total_oLBK)
        relevant_sites = ArData.Data[[entry.time_period=="oLBK" for entry in ArData.Data]]
        relevant_ids = [site.grid_id for site in relevant_sites] 
        SourceData = SoData["oLBK"][1]
        source_ids = [source.grid_id for source in SourceData.Data]
        ids = [relevant_ids; source_ids]
        Ad_total_oLBK ./= num_graphs
        display(plot_network_robustness(Ad_total_oLBK, name="oLBK", output_dir=out_base))
        writedlm(joinpath(out_base, "oLBK_adjacency_probabilities.csv"), [[""] ids'; ids Ad_total_oLBK], ',')
    end
    if !isnothing(Ad_total_fullLBK)
        relevant_sites = ArData.Data[[entry.time_period=="fullLBK" for entry in ArData.Data]]
        relevant_ids = [site.grid_id for site in relevant_sites] 
        SourceData = SoData["fullLBK"][1]
        source_ids = [source.grid_id for source in SourceData.Data]
        ids = [relevant_ids; source_ids]
        Ad_total_fullLBK ./= num_graphs
        display(plot_network_robustness(Ad_total_fullLBK, name="LBK", output_dir=out_base))
        writedlm(joinpath(out_base, "fullLBK_adjacency_probabilities.csv"), [[""] ids'; ids Ad_total_fullLBK], ',')
    end
    if !isnothing(Ad_total_Magda_new)
        relevant_sites = ArData.Data[[entry.time_period=="Magdanew" for entry in ArData.Data]]
        relevant_ids = [site.grid_id for site in relevant_sites] 
        SourceData = SoData["Magdanew"][1]
        source_ids = [source.grid_id for source in SourceData.Data]
        ids = [relevant_ids; source_ids]
        Ad_total_Magda_new ./= num_graphs
        display(plot_network_robustness(Ad_total_Magda_new, name="Magdalenien", output_dir=out_base))
        writedlm(joinpath(out_base, "Magdanew_adjacency_probabilities.csv"), [[""] ids'; ids Ad_total_Magda_new], ',')
    end
    @info "Robustness analysis maps exported successfully."
    return Ad_total_Magda, Ad_total_oLBK, Ad_total_fullLBK, Ad_total_Magda_new
end


"""
    plot_network_robustness(prob_matrix::AbstractMatrix; name::String = "Graph", output_dir::String = pwd())

Plots a "spy plot" graphic interpreting edge probabilities derived from simulated network graphs.

# Arguments
- `prob_matrix::AbstractMatrix`: The density likelihood matrix evaluating historical connectivity points.

# Keywords
- `name::String`: Chart target label.
- `output_dir::String`: Root path for visual rendering objects.

# Output
- `Plots.Plot`: Returns the generated plot object.
"""
function plot_network_robustness(prob_matrix::AbstractMatrix; name::String = "Graph", output_dir::String = pwd())
    N = size(prob_matrix, 1)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    
    for i in 1:N
        for j in 1:N
            if prob_matrix[i, j] > 0.0
                push!(rows, i)
                push!(cols, j)
                push!(vals, prob_matrix[i, j])
            end
        end
    end
    
    p = scatter(
        cols, rows,
        marker_z = vals,
        color = cgrad(:brg),
        markersize = vals .^ 2 .* 6,  
        markerstrokewidth = 0,
        yflip = true,
        framestyle = :box,
        legend = false,
        colorbar = true,
        colorbar_title = "Edge Probability",
        xlabel = "Target Node",
        ylabel = "Source Node",
        title = "Network Robustness (Spy Plot) - $name",
        size = (800, 600)
    )
    
    savefig(p, joinpath(output_dir, "$(name)_Robustness_SpyPlot.png"))
    return p
end