"""
    smoothed_eigenvalue_dist(x::Number, eigenvalues::Vector{<:Number}; sigma::Number = 0.015)

Smoothed eigenvalue distribution for spectral plots using a Gaussian kernel. 
Based on "The Laplacian spectrum of neural networks" (S. de Lange, M. de Reus, M. van den Heuvel).

# Arguments
- `x::Number`: The evaluation point.
- `eigenvalues::Vector{<:Number}`: The array of graph eigenvalues.

# Keywords
- `sigma::Number`: The bandwidth parameter. Defaults to 0.015.

# Output
- `Float64`: The computed smoothed density value at `x`.
"""
@inline smoothed_eigenvalue_dist(x::Number, eigenvalues::Vector{<:Number}; sigma::Number=0.015) = sum(exp.(-(x .- eigenvalues) .^ 2 / (2 * sigma^2))) / sqrt(2 * pi * sigma^2)

"""
    Lorentz_kernel(x::Number, eigenvalues::Vector{<:Number}; sigma::Number = 0.03)

Calculates the Lorentz (Cauchy) kernel density for a given set of eigenvalues.

# Arguments
- `x::Number`: The evaluation point.
- `eigenvalues::Vector{<:Number}`: The array of graph eigenvalues.

# Keywords
- `sigma::Number`: The bandwidth parameter. Defaults to 0.03.

# Output
- `Float64`: The computed density value at `x`.
"""
@inline Lorentz_kernel(x::Number, eigenvalues::Vector{<:Number};sigma::Number=0.03) = sum(sigma./((eigenvalues .-x).^2 .+ sigma^2))

"""
    spectral_plots_comparison(num_vert::Int; num_edges::Int = 0, num_graphs::Int = 1000, kernel::Function = smoothed_eigenvalue_dist, sigma::Number = 0.015)

Creates spectral plots using a smoothed eigenvalue distribution. Eigenvalues are averaged over multiple iterations of simulated random graphs.
Compared baseline models include Soft Random Geometric Graphs (SRGG), Barabási-Albert, and Watts-Strogatz.

# Arguments
- `num_vert::Int`: Number of vertices in the simulated graphs.

# Keywords
- `num_edges::Int`: Target number of edges for the baseline graphs. Defaults to 0.
- `num_graphs::Int`: Number of graph instances to simulate and average. Defaults to 1000.
- `kernel::Function`: Kernel density function to apply. Defaults to `smoothed_eigenvalue_dist`.
- `sigma::Number`: Bandwidth parameter for the kernel. Defaults to 0.015.

# Output
- `Tuple{Plots.Plot, Matrix{Float64}}`: The generated plot object and a matrix containing the averaged eigenvalues for each model.
"""
function spectral_plots_comparison(num_vert::Int; num_edges::Int = 0, num_graphs::Int=1000, kernel::Function=smoothed_eigenvalue_dist, sigma::Number=0.015)
    xx = 0:0.001:2
    p1 = plot(legend=:topright, size=(1200, 600), title=LaTeXString("Spectral Plots"), xlabel=LaTeXString("Eigenvalues"), ylabel=LaTeXString("relative frequency"), left_margin = [8Plots.mm 0Plots.mm], bottom_margin = [4Plots.mm 0Plots.mm])

    averaged_eigs_mat = zeros(num_vert, 3)

    # SRGG Base
    beta = floor(Int, num_vert * (num_vert - 1) / 2.0 * pi / num_edges)
    @debug "Calculated SRGG beta coefficient: $beta"
    
    avg_edges_count = 0.0
    eigensum = zeros(Float64, num_vert)
    @info "Simulating SRGG networks..."
    for j in 1:num_graphs 
        ad_g = SRGG(num_vert, beta)
        avg_edges_count += length(ad_g.nzval) / 2.0
        _, norm_lap_mat = get_matrices(ad_g)
        eigensum .+= real.(eigvals(Matrix(norm_lap_mat)))
        (j % 250 == 0) && @debug "Processed SRGG graph $j / $num_graphs"
    end
    avg_edges_int = ceil(Int, avg_edges_count / num_graphs)
    averaged_eigenvalues = eigensum ./ num_graphs
    plot!(p1, xx, 1.0 .+ kernel.(xx, [averaged_eigenvalues], sigma=sigma), label=LaTeXString("SRGG, \$\\beta=$beta\$"), lw=6)
    averaged_eigs_mat[:, 1] = averaged_eigenvalues

    @debug "Average edges simulated in SRGG: $avg_edges_int"
    num_edges == 0 && (num_edges = avg_edges_int)
   
    # Barabási-Albert Base
    k_ba = ceil(Int, num_vert / 2.0 - sqrt((num_vert / 2.0)^2 - num_edges)) 
    eigensum .= 0.0
    
    @info "Simulating Barabási-Albert networks..."
    for j in 1:num_graphs 
        g = barabasi_albert(num_vert, k_ba)
        _, norm_lap_mat = get_matrices(adjacency_matrix(g))
        eigensum .+= real.(eigvals(Matrix(norm_lap_mat)))
        (j % 250 == 0) && @debug "Processed BA graph $j / $num_graphs"
    end
    
    averaged_eigenvalues = eigensum ./ num_graphs
    plot!(p1, xx, 1.0 .+ kernel.(xx, [averaged_eigenvalues], sigma=sigma), label=LaTeXString("Barabási-Albert"), lw=6)
    averaged_eigs_mat[:, 2] = averaged_eigenvalues

    # Watts-Strogatz Base
    k_ws = ceil(Int, 2.0 * num_edges / num_vert)
    eigensum .= 0.0
    
    @info "Simulating Watts-Strogatz networks..."
    for j in 1:num_graphs
        g = watts_strogatz(num_vert, k_ws, 0.5)
        _, norm_lap_mat = get_matrices(adjacency_matrix(g))
        eigensum .+= real.(eigvals(Matrix(norm_lap_mat)))
        (j % 250 == 0) && @debug "Processed WS graph $j / $num_graphs"
    end
    
    averaged_eigenvalues = eigensum ./ num_graphs
    plot!(p1, xx, 1.0 .+ kernel.(xx, [averaged_eigenvalues], sigma=sigma), label=LaTeXString("Watts-Strogatz"), lw=6)
    averaged_eigs_mat[:, 3] = averaged_eigenvalues

    return p1, averaged_eigs_mat
end

"""
    archaeological_spectral_plots(source_toml::String; num_graphs::Int = 100, kernel_and_par::Tuple{Function, Number} = (smoothed_eigenvalue_dist, 0.015), seed::Any = nothing, incl_sources::Bool = false, double_distance::Bool = false, output_dir::String = pwd(), overwrite::Bool = false)

Generates and evaluates the spectral evolution of archaeological networks as material-driven edges are added incrementally.

# Arguments
- `source_toml::String`: Filename or path to the configuration TOML file.

# Keywords
- `num_graphs::Int`: Number of baseline graphs to simulate for null model comparisons. Defaults to 100.
- `kernel_and_par::Tuple{Function, Number}`: Tuple containing the kernel function and its bandwidth parameter. Defaults to `(smoothed_eigenvalue_dist, 0.015)`.
- `seed::Any`: Random seed for reproducibility. Defaults to nothing.
- `incl_sources::Bool`: Indicates whether to include raw material source nodes. Defaults to false.
- `double_distance::Bool`: Determines connection probability thresholds. Defaults to false.
- `output_dir::String`: Target export folder. Defaults to the active user working directory.
- `overwrite::Bool`: Grants file replacement privileges. Throws an error if set to false and the directory already exists. Defaults to false.

# Output
- `Nothing`: Persists generated plots directly to the local file system.
"""
function archaeological_spectral_plots(source_toml::String; num_graphs::Int = 100, kernel_and_par::Tuple{Function, Number} = (smoothed_eigenvalue_dist, 0.015), seed::Any = nothing, incl_sources::Bool = false, double_distance::Bool = false, output_dir::String = pwd(), overwrite::Bool = false)
    
    out_dir = joinpath(output_dir, string(Dates.today()), "SpectralPlots")
    if isdir(out_dir)
        if !overwrite
            error("Output directory $out_dir already exists. Set `overwrite=true` to proceed.")
        else
            @warn "Overwriting existing spectral plots in $out_dir"
        end
    else
        mkpath(out_dir)
    end

    if seed !== nothing
        Random.seed!(seed)
    end
    
    data_path = isfile(source_toml) ? source_toml : joinpath(dirname(@__DIR__), "DataFolder", source_toml)
    ArData, MaData = read_data_TOML(data_path)
    
    graphen = if incl_sources
        SoData = read_source_data(data_path, MaData)
        compare_all_connect_nodes_sources(ArData, MaData, SoData, double_distance=double_distance, output_dir=out_dir, overwrite=overwrite)
    else
        compare_all_connect_nodes(ArData, MaData)
    end

    # Copy TOML to output for tracking parameters
    target_toml = joinpath(out_dir, basename(data_path))
    cp(data_path, target_toml, force=true)
    open(target_toml, "a") do io
        write(io, "\nmaximum_distance = $maximum_distance\nmaterial_multiplier=$material_multiplier")
    end

    kernel = kernel_and_par[1]
    sigma = Float64(kernel_and_par[2])

    for i in eachindex(graphen)
        time_period = graphen[i][2]
        @info "Evaluating Spectral Evolution for Time Period: $time_period"
        
        adj_mats = graphen[i][1]
        num_verts = size(adj_mats[1], 1)
        
        # 1. Filter out redundant matrices (e.g., when default == group-default)
        valid_indices = Int[]
        for j in 1:length(adj_mats)
            if j == 2 && isapprox(adj_mats[1], adj_mats[2])
                continue
            end
            push!(valid_indices, j)
        end

        # 2. Determine indices for the 3 snapshot plots (Start, Middle, End)
        idx_start = valid_indices[min(3, length(valid_indices))] 
        idx_mid   = valid_indices[ceil(Int, length(valid_indices) / 2.0)]
        idx_end   = valid_indices[end]
        snapshot_indices = [idx_start, idx_mid, idx_end]

        # Tracking arrays for distance trajectories
        euclidean_SRGG = Float64[]
        euclidean_BA   = Float64[]
        euclidean_WS   = Float64[]
        kl_SRGG        = Float64[]
        kl_BA          = Float64[]
        kl_WS          = Float64[]
        distances_SRGG = Float64[]
        distances_BA   = Float64[]
        distances_WS   = Float64[]
        wasserstein_SRGG = Float64[]
        wasserstein_BA   = Float64[]
        wasserstein_WS   = Float64[]
        x_labels       = String[]

        # 3. Iterate over valid matrices
        @info "Calculating Spectral Evolution Stages..."
        for j_mat in valid_indices
            adj_mat = adj_mats[j_mat]
            dropzeros!(adj_mat) 
            num_edges = Int(nnz(adj_mat) / 2)
            @debug "Number of edges in archaeological network at stage $j_mat: $num_edges"
            
            # Generate labels for X-axis
            mat_string = if j_mat == 1
                "Default Dist."
            elseif j_mat == 2
                "Group Dist."
            else
                "+ $(j_mat-2) Mat."
            end
            push!(x_labels, mat_string)

            # Compute null models
            plot_Mat, avg_eigs_mat = spectral_plots_comparison(num_verts, num_edges=num_edges, num_graphs=num_graphs, kernel=kernel, sigma=sigma)
            
            # Compute eigenvalues of the archaeological graph
            _, norm_lap_Mat = get_matrices(adj_mat)
            arch_eigs = real.(eigvals(Matrix(norm_lap_Mat)))

            # Euclidean Distances
            euclid_SRGG = norm(arch_eigs .- avg_eigs_mat[:, 1])
            euclid_BA   = norm(arch_eigs .- avg_eigs_mat[:, 2])
            euclid_WS   = norm(arch_eigs .- avg_eigs_mat[:, 3])

            # Continuous PDF calculations for Divergence & Overlap
            step_size = 0.001
            xx = 0.0:step_size:2.0
            
            pdf_arch = kernel.(xx, [arch_eigs], sigma=sigma) ./ num_verts
            pdf_SRGG = kernel.(xx, [avg_eigs_mat[:, 1]], sigma=sigma) ./ num_verts
            pdf_BA   = kernel.(xx, [avg_eigs_mat[:, 2]], sigma=sigma) ./ num_verts
            pdf_WS   = kernel.(xx, [avg_eigs_mat[:, 3]], sigma=sigma) ./ num_verts

            # Kullback-Leibler Divergence
            eps_val = 1e-10 
            D_KL_SRGG = sum((pdf_arch .+ eps_val) .* log.((pdf_arch .+ eps_val) ./ (pdf_SRGG .+ eps_val))) * step_size
            D_KL_BA   = sum((pdf_arch .+ eps_val) .* log.((pdf_arch .+ eps_val) ./ (pdf_BA .+ eps_val))) * step_size
            D_KL_WS   = sum((pdf_arch .+ eps_val) .* log.((pdf_arch .+ eps_val) ./ (pdf_WS .+ eps_val))) * step_size

            # Spectral Overlap Distance (1 - Cosine Similarity)
            overlap_SRGG = sum(pdf_arch .* pdf_SRGG) / (norm(pdf_arch) * norm(pdf_SRGG))
            overlap_BA   = sum(pdf_arch .* pdf_BA)   / (norm(pdf_arch) * norm(pdf_BA))
            overlap_WS   = sum(pdf_arch .* pdf_WS)   / (norm(pdf_arch) * norm(pdf_WS))
            
            D_Cos_SRGG = 1.0 - overlap_SRGG
            D_Cos_BA   = 1.0 - overlap_BA
            D_Cos_WS   = 1.0 - overlap_WS

            # Wasserstein-1 Distance (Area between CDFs)
            cdf_arch = cumsum(pdf_arch) .* step_size
            cdf_SRGG = cumsum(pdf_SRGG) .* step_size
            cdf_BA   = cumsum(pdf_BA)   .* step_size
            cdf_WS   = cumsum(pdf_WS)   .* step_size

            W1_SRGG = sum(abs.(cdf_arch .- cdf_SRGG)) * step_size
            W1_BA   = sum(abs.(cdf_arch .- cdf_BA))   * step_size
            W1_WS   = sum(abs.(cdf_arch .- cdf_WS))   * step_size

            # Store computed metrics
            push!(euclidean_SRGG, euclid_SRGG)
            push!(euclidean_BA,   euclid_BA)
            push!(euclidean_WS,   euclid_WS)
            push!(kl_SRGG, D_KL_SRGG)
            push!(kl_BA,   D_KL_BA)
            push!(kl_WS,   D_KL_WS)
            push!(distances_SRGG, D_Cos_SRGG)
            push!(distances_BA,   D_Cos_BA)
            push!(distances_WS,   D_Cos_WS)
            push!(wasserstein_SRGG, W1_SRGG)
            push!(wasserstein_BA,   W1_BA)
            push!(wasserstein_WS,   W1_WS)
            
            # 4. Save Snapshots
            if j_mat in snapshot_indices
                snapshot_title = if j_mat == idx_start
                    "Baseline Spatial Network ($mat_string)"
                elseif j_mat == idx_mid
                    "Intermediate Network ($mat_string)"
                else
                    "Final Material Network ($mat_string)"
                end

                plot!(plot_Mat, xx, 1.0 .+ kernel.(xx, [arch_eigs], sigma=sigma), label=LaTeXString("Arch. graph"), lw=6, color=:black, title=snapshot_title)
                savefig(plot_Mat, joinpath(out_dir, "$(time_period)_snapshot_stage_$(j_mat).png"))
            end
        end

        # 5. Export Trajectory Plots at the end of the epoch
        @info "Exporting Evolution Plots for Time Period: $time_period"

        # Euclidean Evolution
        traj_plot_euclidean = plot(title="Euclidean Evolution ($time_period)", xlabel="Added Materials", ylabel=LaTeXString("Euclidean Distance"), legend=:topright, size=(1000, 600), bottom_margin=8Plots.mm, left_margin=4Plots.mm)
        plot!(traj_plot_euclidean, 1:length(valid_indices), euclidean_SRGG, label="Euclid: SRGG", marker=:circle, lw=3, color=:blue)
        plot!(traj_plot_euclidean, 1:length(valid_indices), euclidean_BA,   label="Euclid: Barabási-Albert", marker=:square, lw=3, color=:red)
        plot!(traj_plot_euclidean, 1:length(valid_indices), euclidean_WS,   label="Euclid: Watts-Strogatz", marker=:diamond, lw=3, color=:green)
        xticks!(traj_plot_euclidean, 1:length(valid_indices), x_labels, rotation=45)
        savefig(traj_plot_euclidean, joinpath(out_dir, "$(time_period)_euclidean_evolution.png"))

        # Wasserstein Evolution
        traj_plot_w1 = plot(title="Wasserstein Evolution ($time_period)", xlabel="Added Materials", ylabel=LaTeXString("Wasserstein Distance \$W_1\$"), legend=:topright, size=(1000, 600), bottom_margin=8Plots.mm, left_margin=4Plots.mm)
        plot!(traj_plot_w1, 1:length(valid_indices), wasserstein_SRGG, label="W1: SRGG", marker=:circle, lw=3, color=:blue)
        plot!(traj_plot_w1, 1:length(valid_indices), wasserstein_BA,   label="W1: Barabási-Albert", marker=:square, lw=3, color=:red)
        plot!(traj_plot_w1, 1:length(valid_indices), wasserstein_WS,   label="W1: Watts-Strogatz", marker=:diamond, lw=3, color=:green)
        xticks!(traj_plot_w1, 1:length(valid_indices), x_labels, rotation=45)
        savefig(traj_plot_w1, joinpath(out_dir, "$(time_period)_wasserstein_evolution.png"))

        # Kullback-Leibler Evolution
        traj_plot_kl = plot(title="Kullback-Leibler Evolution ($time_period)", xlabel="Added Materials", ylabel=LaTeXString("KL Divergence"), legend=:topright, size=(1000, 600), bottom_margin=8Plots.mm, left_margin=4Plots.mm)
        plot!(traj_plot_kl, 1:length(valid_indices), kl_SRGG, label="KL: SRGG", marker=:circle, lw=3, color=:blue)
        plot!(traj_plot_kl, 1:length(valid_indices), kl_BA,   label="KL: Barabási-Albert", marker=:square, lw=3, color=:red)
        plot!(traj_plot_kl, 1:length(valid_indices), kl_WS,   label="KL: Watts-Strogatz", marker=:diamond, lw=3, color=:green)
        xticks!(traj_plot_kl, 1:length(valid_indices), x_labels, rotation=45)
        savefig(traj_plot_kl, joinpath(out_dir, "$(time_period)_kullback_leibler_evolution.png"))

        # Spectral Overlap Distance (Cosine)
        traj_plot_overlap = plot(title="Evolution of Spectral Distance ($time_period)", xlabel="Added Materials", ylabel=LaTeXString("Spectral Distance \$D_{spec}\$"), legend=:topright, size=(1000, 600), bottom_margin=8Plots.mm, left_margin=4Plots.mm)
        plot!(traj_plot_overlap, 1:length(valid_indices), distances_SRGG, label="Cosine: SRGG", marker=:circle, lw=3, color=:blue)
        plot!(traj_plot_overlap, 1:length(valid_indices), distances_BA, label="Cosine: Barabási-Albert", marker=:square, lw=3, color=:red)
        plot!(traj_plot_overlap, 1:length(valid_indices), distances_WS, label="Cosine: Watts-Strogatz", marker=:diamond, lw=3, color=:green)
        xticks!(traj_plot_overlap, 1:length(valid_indices), x_labels, rotation=45)
        savefig(traj_plot_overlap, joinpath(out_dir, "$(time_period)_spectral_overlap_evolution.png"))
    end
    
    @info "Spectral analysis complete. All plots saved to $out_dir"
    return nothing
end


