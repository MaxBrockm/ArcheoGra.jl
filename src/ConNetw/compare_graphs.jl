"""
    matrix_spectral_distance(matrix1::AbstractMatrix, matrix2::AbstractMatrix, which::Symbol; k_first::Int = 0)

Calculates the spectral distance between two graphs using their corresponding matrices and the Euclidean norm.

# Arguments
- `matrix1::AbstractMatrix`: First matrix to compare.
- `matrix2::AbstractMatrix`: Second matrix to compare.
- `which::Symbol`: Specification for the eigenvalues to calculate (e.g., `:SM` for smallest magnitude).

# Keywords
- `k_first::Int`: Number of eigenvalues to compute. If `0`, all but one are computed.

# Output
- `Float64`: The computed spectral distance.
"""
function matrix_spectral_distance(matrix1::AbstractMatrix, matrix2::AbstractMatrix, which::Symbol; k_first::Int = 0)
    num_eigs = k_first != 0 ? k_first : size(matrix1, 1) - 1
    eigenvalues1 = sort(eigs(matrix1, which=which, nev=num_eigs)[1])
    eigenvalues2 = sort(eigs(matrix2, which=which, nev=num_eigs)[1])
    return sqrt(sum((eigenvalues1 .- eigenvalues2) .^ 2))
end



"""
    effective_resistance(w_lap_mat::SparseMatrixCSC; klargest::Union{Int, Nothing} = nothing)

Determines the effective resistance of an underlying weighted graph using the eigenvalues of its weighted Laplacian matrix.
Effective resistance serves as a measure of connectedness.

In order to achieve this, the weighted adjacency matrix uses the inverse edge lengths as its (nonzero) entries. 
Shorter or additional edges are favorable, while longer or missing edges are penalized.

# Arguments
- `w_lap_mat::SparseMatrixCSC`: Weighted Laplacian matrix of the graph.

# Keywords
- `klargest::Union{Int, Nothing}`: Optional number of largest eigenvalues to include. If the smallest eigenvalues are negligible, setting this can accelerate the calculation.

# Output
- `Float64`: The computed effective resistance.
"""
function effective_resistance(w_lap_mat::SparseMatrixCSC; klargest::Union{Int, Nothing} = nothing)
    n = size(w_lap_mat, 1)
    # Start at 2 to exclude the zero eigenvalue of a connected graph
    start_idx = isnothing(klargest) ? 2 : klargest
    EV = eigvals(Symmetric(Matrix(w_lap_mat)), start_idx:n)
    return n * sum(1 ./ EV)
end


"""
    feature_extraction(Ad_weighted::SparseMatrixCSC, edge_comparison::Bool)

Determines the degrees, two-hop degrees, and (optionally) edge lengths of a given weighted graph.
The weights must correspond to the edge lengths of the graph.

# Arguments
- `Ad_weighted::SparseMatrixCSC`: Weighted adjacency matrix of the graph.
- `edge_comparison::Bool`: If `true`, non-zero edge values are also extracted.

# Output
- `Tuple{Vector{Int}, Vector{Float64}, Vector{Float64}}`: Contains the degree vector, two-hop degree vector, and edge lengths (empty if `edge_comparison` is `false`).
"""
function feature_extraction(Ad_weighted::SparseMatrixCSC, edge_comparison::Bool)
    degree_vector = vec(count(>(0), Ad_weighted, dims=2))
    feature_two_hop_degree = [1.0 / degree_vector[i] * sum(degree_vector[Ad_weighted[:, i].nzind]) for i in axes(Ad_weighted, 1)]
    
    edge_vals = edge_comparison ? Ad_weighted.nzval : Float64[]
    return (degree_vector, feature_two_hop_degree, edge_vals)
end


"""
    canberra_distance(d1::AbstractVector, d2::AbstractVector)

Calculates the Canberra distance between two vectors. 
This method scales with the absolute value of the entries, making it highly sensitive to changes near zero.

# Arguments
- `d1::AbstractVector`: First feature vector.
- `d2::AbstractVector`: Second feature vector.

# Output
- `Float64`: The Canberra distance.
"""
@inline canberra_distance(d1::AbstractVector, d2::AbstractVector) = sum(abs.(d1 - d2) ./ (d1 + d2))


"""
    feature_vector_comparison(agg_vec1::AbstractVector, agg_vec2::AbstractVector)

Calculates different measures of distance (Canberra, Euclidean, Difference, Absolute Difference) between two vectors.

# Arguments
- `agg_vec1::AbstractVector`: First aggregated feature vector.
- `agg_vec2::AbstractVector`: Second aggregated feature vector.

# Output
- `Tuple{Float64, Float64, Vector{Float64}, Vector{Float64}}`: Distance metrics and difference vectors.
"""
function feature_vector_comparison(agg_vec1::AbstractVector, agg_vec2::AbstractVector)
    can_dist = canberra_distance(agg_vec1, agg_vec2)
    euclid_dist = norm(agg_vec1 .- agg_vec2, 2)
    diff = agg_vec1 .- agg_vec2
    abs_diff = abs.(agg_vec1 .- agg_vec2)
    return (can_dist, euclid_dist, diff, abs_diff)
end


"""
    feature_comparison(Ad_weighted1::SparseMatrixCSC, comp_graph::Function; num_graphs::Int = 1000, titleaddon::String = "", edge_comparison::Bool = true)

Applies a feature comparison, comparing edge length, degree, and two-hop degree against a generated reference graph.
Features are extracted and then aggregated into four statistical categories.

# Arguments
- `Ad_weighted1::SparseMatrixCSC`: Weighted adjacency matrix of the primary graph.
- `comp_graph::Function`: Generator function producing the comparison graph's adjacency matrix.

# Keywords
- `num_graphs::Int`: Number of reference graphs to generate for averaging.
- `titleaddon::String`: Suffix added to the plot titles.
- `edge_comparison::Bool`: Indicates whether to include edge length comparisons.

# Output
- `Nothing`: Displays visual plots and prints metrics directly.
"""
function feature_comparison(Ad_weighted1::SparseMatrixCSC, comp_graph :: Function; num_graphs::Int = 1000, titleaddon="", edge_comparison=true)
    fvec1 = feature_extraction(Ad_weighted1, edge_comparison)
    
    num_features = edge_comparison ? 3 : 2
    aggregated_mat1 = zeros(4, num_features)
    aggregated_mat2 = zeros(4, num_features)

    for i in 1:num_features
        aggregated_mat1[:, i] = feature_aggregation(fvec1[i])
    end

    bins_deg = 0:1:maximum(fvec1[1])
    bins_two_hop_deg = 0:1:ceil(Int, maximum(fvec1[2]))
    bins_edge_length = edge_comparison ? LinRange(0, maximum(fvec1[3]), 36) : LinRange(0, 1, 2)

    average_degree_freq = zeros(length(bins_deg))
    average_two_hop_freq = zeros(length(bins_two_hop_deg))
    average_edge_length_freq = edge_comparison ? zeros(length(bins_edge_length)) : Float64[]
    
    freq_deg1 = [[length(fvec1[1][bins_deg[i] .< fvec1[1] .<= bins_deg[i+1]]) for i in 1:length(bins_deg)-1]; length(fvec1[1][fvec1[1] .> bins_deg[end]])]
    freq_two_hop_deg1 = [[length(fvec1[2][bins_two_hop_deg[i] .< fvec1[2] .<= bins_two_hop_deg[i+1]]) for i in 1:length(bins_two_hop_deg)-1]; length(fvec1[2][fvec1[2] .> bins_two_hop_deg[end]])]
    freq_edge_length1 = edge_comparison ? [[length(fvec1[3][bins_edge_length[i] .< fvec1[3] .<= bins_edge_length[i+1]]) for i in 1:length(bins_edge_length)-1]; length(fvec1[3][fvec1[3] .> bins_edge_length[end]])] : Int[]

    @showprogress for _ in 1:num_graphs
        Ad_weighted2 = comp_graph()
        fvec2 = feature_extraction(Ad_weighted2, edge_comparison)

        average_degree_freq .+= [[length(fvec2[1][bins_deg[i] .< fvec2[1] .<= bins_deg[i+1]]) for i in 1:length(bins_deg)-1]; length(fvec2[1][fvec2[1] .> bins_deg[end]])]
        average_two_hop_freq .+= [[length(fvec2[2][bins_two_hop_deg[i] .< fvec2[2] .<= bins_two_hop_deg[i+1]]) for i in 1:length(bins_two_hop_deg)-1]; length(fvec2[2][fvec2[2] .> bins_two_hop_deg[end]])]
        
        if edge_comparison
            average_edge_length_freq .+= [[length(fvec2[3][bins_edge_length[i] .< fvec2[3] .<= bins_edge_length[i+1]]) for i in 1:length(bins_edge_length)-1]; length(fvec2[3][fvec2[3] .> bins_edge_length[end]])]
        end

        for i in 1:num_features
            aggregated_mat2[:, i] .+= feature_aggregation(fvec2[i])
        end
    end

    average_degree_freq ./= num_graphs
    average_two_hop_freq ./= num_graphs
    aggregated_mat2 ./= num_graphs

    histograms = Plots.Plot[]
    comparison_histograms = Plots.Plot[]

    push!(histograms, bar(bins_deg, [freq_deg1 average_degree_freq], xticks=bins_deg[1:2:end], fillcolor=[:blue :red], xlabel="degree", ylabel="frequency", bar_width=0.85, title="Degree frequency in both graphs: $titleaddon", label=["graph 1" "graph 2"], fillalpha=[1 0.5], leftmargin=8Plots.mm))
    push!(histograms, bar(bins_two_hop_deg, [freq_two_hop_deg1 average_two_hop_freq], xticks=bins_two_hop_deg[1:2:end], fillcolor=[:blue :red], xlabel="two-hop degree", bar_width=0.85, ylabel="frequency", title="Two hop degree frequency in both graph: $titleaddon", label=["graph 1" "graph 2"], fillalpha=[1 0.5], leftmargin=8Plots.mm))

    if edge_comparison
        average_edge_length_freq ./= num_graphs
        push!(histograms, bar(bins_edge_length, [freq_edge_length1 average_edge_length_freq], xticks=trunc.(bins_edge_length[1:5:end], digits=2), fillcolor=[:blue :red], xlabel="edge length", bar_width=bins_edge_length[2]-bins_edge_length[1], ylabel="frequency", title="Edge length frequency in both graphs: $titleaddon", label=["graph 1" "graph 2"], fillalpha=[1 0.5], leftmargin=8Plots.mm))

        h_edge_length_diff = bar(bins_edge_length, min.(freq_edge_length1 .- average_edge_length_freq, 0.0), xticks=trunc.(bins_edge_length[1:5:end], digits=2), bar_width=bins_edge_length[2]-bins_edge_length[1], fillcolor=:red, legend=true, title="Difference in edge length: $titleaddon", xlabel="degree", ylabel="difference in frequency", label="higher freq. in graph 2", leftmargin=8Plots.mm)
        bar!(h_edge_length_diff, bins_edge_length, max.(freq_edge_length1 .- average_edge_length_freq, 0.0), fillcolor=:blue, label="higher freq. in graph 1")
        push!(comparison_histograms, h_edge_length_diff)
    end

    h_deg_diff = bar(bins_deg, min.(freq_deg1 .- average_degree_freq, 0.0), xticks=bins_deg[1:2:end], bar_width=bins_deg[2] - bins_deg[1], fillcolor=:red, legend=:topright, title="Difference in degree $titleaddon", xlabel="degree", ylabel="difference in frequency", label="higher freq. in graph 2", leftmargin=8Plots.mm)
    bar!(h_deg_diff, bins_deg, max.(freq_deg1 .- average_degree_freq, 0.0), fillcolor=:blue, label="higher freq. in graph 1")

    h_two_hop_deg_diff = bar(bins_two_hop_deg, min.(freq_two_hop_deg1 .- average_two_hop_freq, 0.0), xticks=bins_two_hop_deg[1:2:end], bar_width=bins_two_hop_deg[2] - bins_two_hop_deg[1], fillcolor=:red, legend=true, title="Difference in two hop degree: $titleaddon", xlabel="degree", ylabel="difference in frequency", label="higher freq. in graph 2", leftmargin=8Plots.mm)
    bar!(h_two_hop_deg_diff, bins_two_hop_deg, max.(freq_two_hop_deg1 .- average_two_hop_freq, 0.0), fillcolor=:blue, label="higher freq. in graph 1")
    
    push!(comparison_histograms, h_deg_diff, h_two_hop_deg_diff)
    
    display(plot(histograms..., layout=Plots.grid(length(histograms), 1, heights=ones(length(histograms))/length(histograms)), size=(750, 450*length(histograms))))
    display(plot(comparison_histograms..., layout=Plots.grid(length(histograms), 1, heights=ones(length(histograms))/length(histograms)), size=(750, 450*length(histograms))))
    
    pretty_table(edge_comparison ? aggregated_mat1 : hcat(aggregated_mat1, fill("-", 4)), header=(["Graph 1", "", ""], ["degree", "two-hop degree", "edge length"]), row_labels=["mean", "median", "standard deviation", "maximum"], crop=:none)
    pretty_table(edge_comparison ? aggregated_mat2 : hcat(aggregated_mat2, fill("-", 4)), header=(["Graph 2", "", ""], ["degree", "two-hop degree", "edge length"]), row_labels=["mean", "median", "standard deviation", "maximum"], crop=:none)

    comparison_vec = feature_vector_comparison(vec(aggregated_mat1), vec(aggregated_mat2))

    @info "Canberra distance: $(comparison_vec[1])"
    @info "Euclidean distance: $(comparison_vec[2])"

    difference_mat1 = reshape(comparison_vec[3], (4, num_features))
    difference_mat2 = reshape(comparison_vec[4], (4, num_features))
    
    pretty_table(edge_comparison ? difference_mat1 : hcat(difference_mat1, fill("-", 4)), header=(["Difference ", "(of aggregated values)", ""], ["degree", "two-hop degree", "edge length"]), row_labels=["mean", "median", "standard deviation", "maximum"], crop=:none)
    pretty_table(edge_comparison ? difference_mat2 : hcat(difference_mat2, fill("-", 4)), header=(["Absolute difference ", "(of aggregated values)", ""], ["degree", "two-hop degree", "edge length"]), row_labels=["mean", "median", "standard deviation", "maximum"], crop=:none)
    return nothing
end


###########################################################
###########################################################
##### Service functions; possibly outdated
###########################################################
###########################################################

"""
    SRGG(n::Int, beta::Real)

Creates a Soft Random Geometric Graph (SRGG) with `n` vertices. 
Vertices are uniformly distributed in [0,1]^2 and connected using a Rayleigh connection function with parameter `beta`. 
The function ensures the final graph is fully connected.

# Arguments
- `n::Int`: Number of vertices.
- `beta::Real`: Parameter for the connection probability.

# Output
- `SparseMatrixCSC`: The weighted adjacency matrix of the resulting graph.
"""
function SRGG(n::Int, beta::Real)
    x_coords = rand(n)
    y_coords = rand(n)
    distance_mat = zeros(n, n)
    
    I = Int[]
    J = Int[]
    sizehint!(I, n * 4) 
    sizehint!(J, n * 4)
    
    @inbounds for i in 1:n
        for j in (i+1):n
            sq_dist = (x_coords[i] - x_coords[j])^2 + (y_coords[i] - y_coords[j])^2
            dist = sqrt(sq_dist)
            
            distance_mat[i, j] = dist
            distance_mat[j, i] = dist
            
            prob = exp(-beta * sq_dist)
            if prob > rand()
                push!(I, i, j)
                push!(J, j, i)
            end
        end
    end
    
    adjacency_mat = sparse(I, J, ones(length(I)), n, n)
    connect_components!(adjacency_mat, components(adjacency_mat)..., distance_mat)
    
    return distance_mat .* adjacency_mat
end

"""
    get_matrices(adjacency_matrix::AbstractMatrix)

Determines the Laplacian and normalized Laplacian matrices of a graph from its adjacency matrix.

# Arguments
- `adjacency_matrix::AbstractMatrix`: The adjacency matrix of the graph.

# Output
- `Tuple{SparseMatrixCSC, SparseMatrixCSC}`: The Laplacian matrix and the normalized Laplacian matrix.
"""
function get_matrices(adjacency_matrix::AbstractMatrix)
    degrees = vec(sum(adjacency_matrix, dims=1))
    laplacian_matrix = spdiagm(degrees) - adjacency_matrix
    normalized_laplacian_matrix = I - spdiagm(1 ./ degrees) * adjacency_matrix
    return (laplacian_matrix, normalized_laplacian_matrix)
end

### Flagged for CleanUp
"""
    produce_results(graph_adj_mat::SparseMatrixCSC; output_dir::String = pwd(), overwrite::Bool = false)

Applies functions for graph comparison: spectral plots, spectral distances, cluster analysis, feature comparison, and effective resistance.
The provided graph is compared against standard random models (SRGG, Barabási-Albert, and Watts-Strogatz).

# Arguments
- `graph_adj_mat::SparseMatrixCSC`: Weighted adjacency matrix of the graph to be analyzed.

# Keywords
- `output_dir::String`: Directory to save the spectral plot.
- `overwrite::Bool`: If `false` (default), prevents overwriting existing directories.
"""
function produce_results(graph_adj_mat::SparseMatrixCSC; output_dir::String = pwd(), overwrite::Bool = false)
    n = size(graph_adj_mat, 1)
    m = Int(length(graph_adj_mat.nzval) / 2)

    adj_I, adj_J, _ = findnz(graph_adj_mat)
    non_weighted_adj = sparse(adj_I, adj_J, ones(2*m))

    lap_mat, norm_lap_mat = get_matrices(graph_adj_mat)
    norm_lap_EV = real.(eigvals(Matrix(norm_lap_mat)))
    lap_EV = real.(eigvals(Matrix(lap_mat)))
    adj_EV = real.(eigvals(Matrix(graph_adj_mat)))

    nw_lap_mat, nw_norm_lap_mat = get_matrices(non_weighted_adj)
    nw_norm_lap_EV = real.(eigvals(Matrix(nw_norm_lap_mat)))
    nw_lap_EV = real.(eigvals(Matrix(nw_lap_mat)))
    nw_adj_EV = real.(eigvals(Matrix(non_weighted_adj)))

    @info "Spectral plots in progress."
    spec_plot, averaged_eigenvalues_norm_lap = spectral_plots_comparison(n, num_edges=m)
    xx = 0:0.001:2
    plot!(spec_plot, xx, smoothed_eigenvalue_dist.(xx, [norm_lap_EV]), label=LaTeXString("Compared graph"), lw=4)
    display(spec_plot)
    
    # OS-independent output directory creation and check
    out_path = joinpath(output_dir, string(Dates.today()), "CompareGraphs")
    if isdir(out_path)
        if !overwrite
            error("Folder '$out_path' already exists! To overwrite files, use the argument `overwrite=true`.")
        else
            @info "Folder already exists. Plot will be overwritten (overwrite=true)."
        end
    else
        mkpath(out_path)
    end
    savefig(spec_plot, joinpath(out_path, "spectral_plot.pdf"))
    
    ### Spectral distances & effective resistance for SRGG
    @info "Spectral distances in progress."
    effective_resistance_SRGG = 0.0
    average_EV_SRGG_adj = zeros(n)
    average_EV_SRGG_lap = zeros(n)
    average_EV_BA_adj = zeros(n)
    average_EV_BA_lap = zeros(n)
    average_EV_WS_adj = zeros(n)
    average_EV_WS_lap = zeros(n)
    
    k_BA = ceil(Int, n / 2 - sqrt((n / 2)^2 - m))
    k_WS = ceil(Int, 2 * m / n)
    
    @showprogress for _ in 1:100
        adj_SRGG = SRGG(n, 15)
        lap_mat_SRGG, _ = get_matrices(adj_SRGG)
        
        average_EV_SRGG_adj .+= eigvals(Matrix(adj_SRGG))
        eigvals_lapmat_SRGG = eigvals(Matrix(lap_mat_SRGG))
        average_EV_SRGG_lap .+= eigvals_lapmat_SRGG    
        
        effective_resistance_SRGG += n * sum(1 ./ eigvals_lapmat_SRGG)
        
        adj_BA = adjacency_matrix(barabasi_albert(n, k_BA))
        lap_mat_BA, _ = get_matrices(adj_BA)
        average_EV_BA_adj .+= eigvals(Matrix(adj_BA))
        average_EV_BA_lap .+= eigvals(Matrix(lap_mat_BA))
        
        adj_WS = adjacency_matrix(watts_strogatz(n, k_WS, 0.5))
        lap_mat_WS, _ = get_matrices(adj_WS)
        average_EV_WS_adj .+= eigvals(Matrix(adj_WS))
        average_EV_WS_lap .+= eigvals(Matrix(lap_mat_WS))
    end
   
    average_EV_BA_adj ./= 100
    average_EV_BA_lap ./= 100
    average_EV_SRGG_adj ./= 100
    average_EV_SRGG_lap ./= 100
    average_EV_WS_adj ./= 100
    average_EV_WS_lap ./= 100
    
    spec_dist_mat = zeros(3, 3)
    spec_dist_mat[1, :] = [norm(adj_EV .- average_EV_SRGG_adj, 2), norm(nw_adj_EV .- average_EV_BA_adj, 2), norm(nw_adj_EV .- average_EV_WS_adj, 2)]
    spec_dist_mat[2, :] = [norm(lap_EV .- average_EV_SRGG_lap, 2), norm(nw_lap_EV .- average_EV_BA_lap, 2), norm(nw_lap_EV .- average_EV_WS_lap, 2)]
    spec_dist_mat[3, :] = [norm(norm_lap_EV .- averaged_eigenvalues_norm_lap[:, 1], 2), norm(nw_norm_lap_EV .- averaged_eigenvalues_norm_lap[:, 2], 2), norm(nw_norm_lap_EV .- averaged_eigenvalues_norm_lap[:, 3], 2)]
    
    @info "Spectral distances to averaged eigenvalue distributions:"
    pretty_table(spec_dist_mat, header=["SRGG", "Barabási-Albert", "Watts-Strogatz"], row_labels=["adjacency matrix", "laplacian matrix", "normalized laplacian matrix"], crop=:none)

    ### Cluster Analysis
    @info "Cluster analysis in progress."
    lap_EV = eigvals(Matrix(lap_mat))
    @debug "Eigenvalues computed successfully."
    
    num_clusters = identify_spectral_gap(lap_EV)
    @debug "Number of clusters computed: $num_clusters."
    
    cluster = 1:n
    if num_clusters != 1
        cluster = spectral_kmeans(lap_mat, num_clusters)
    end

    ### Feature Comparison
    @info "Feature comparison to SRGG in progress."
    beta = 15
    comp_SRGG() = SRGG(n, beta)
    feature_comparison(graph_adj_mat, comp_SRGG, titleaddon="SRGG, β=$beta")

    @info "Feature comparison to Barabási-Albert graph in progress."
    comp_barabasi_albert() = adjacency_matrix(barabasi_albert(n, k_BA))
    feature_comparison(graph_adj_mat, comp_barabasi_albert, titleaddon="Barabási-Albert", edge_comparison=false)

    @info "Feature comparison to Watts-Strogatz in progress."
    comp_watts_strogatz() = adjacency_matrix(watts_strogatz(n, k_WS, 0.5))
    feature_comparison(graph_adj_mat, comp_watts_strogatz, titleaddon="Watts-Strogatz", edge_comparison=false)

    ### Effective Resistance
    effective_resistance_SRGG /= 100
    graph_effective_resistance = effective_resistance(lap_mat)
    
    @info "Effective resistance computed:"
    @info "-> SRGG reference: $effective_resistance_SRGG"
    @info "-> Analyzed graph: $graph_effective_resistance"
    
    return nothing
end

