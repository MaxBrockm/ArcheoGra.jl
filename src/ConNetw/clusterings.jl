###############################################################
###### MCL Utility
###############################################################

"""
    extract_clusters(M::SparseMatrixCSC)

Given a converged MCL/R-MCL stochastic matrix, reads off cluster assignments from the attractor (first nonzero) in each column.

# Arguments
- `M::SparseMatrixCSC`: The converged stochastic matrix.

# Output
- `Tuple{Vector{Vector{Int}}, Vector{Int}}`: A tuple containing the `clusters` (vector of vectors of vertex indices) and `node_labels`.
"""
function extract_clusters(M::SparseMatrixCSC)
    n = size(M, 2)
    node_labels = Vector{Int}(undef, n)
    for j in 1:n
        node_labels[j] = M[:, j].nzind[1]
    end
    unique_labels = unique(node_labels)
    label_map = Dict(l => i for (i, l) in enumerate(unique_labels))
    for j in 1:n
        node_labels[j] = label_map[node_labels[j]]
    end
    k = length(unique_labels)
    clusters = [Int[] for _ in 1:k]
    for j in 1:n
        push!(clusters[node_labels[j]], j)
    end
    return clusters, node_labels
end


"""
    mcl_inflate_prune!(M::SparseMatrixCSC, r::Int; prune_ratio::Float64 = 0.25)

In-place inflate (element-wise power `r`), column-normalise, prune small entries, and drop zeros.

# Arguments
- `M::SparseMatrixCSC`: The matrix to be manipulated.
- `r::Int`: The inflation power.

# Keywords
- `prune_ratio::Float64`: Ratio used to determine the pruning threshold based on average column values.

# Output
- `SparseMatrixCSC`: The mutated matrix `M`.
"""
function mcl_inflate_prune!(M::SparseMatrixCSC, r::Int; prune_ratio::Float64 = 0.25)
    M .= M .^ r
    col_sums = vec(sum(M; dims=1))
    for j in axes(M, 2)
        col_sums[j] != 0 && (@views M[:, j] ./= col_sums[j])
    end
    for j in axes(M, 2)
        nzv = M[:, j].nzval
        isempty(nzv) && continue
        avg_val = sum(nzv) / length(nzv)
        threshold = prune_ratio * avg_val
        for idx in eachindex(M[:, j].nzind)
            M[M[:, j].nzind[idx], j] < threshold && (M[M[:, j].nzind[idx], j] = 0.0)
        end
    end
    dropzeros!(M)
    return M
end

"""
    column_normalise!(M::SparseMatrixCSC)

Normalise each column of a sparse matrix to sum to 1.

# Arguments
- `M::SparseMatrixCSC`: The matrix to normalise.

# Output
- `SparseMatrixCSC`: The mutated matrix `M`.
"""
function column_normalise!(M::SparseMatrixCSC)
    for j in axes(M, 2)
        s = sum(M[:, j])
        s != 0 && (@views M[:, j] ./= s)
    end
    return M
end


"""
    build_coarse_adjacency(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})

Constructs a coarsened adjacency matrix representing the connectivity between distinct clusters.

# Arguments
- `Ad::SparseMatrixCSC`: Original adjacency matrix.
- `clusters::Vector{Vector{Int}}`: Vertex assignments for each cluster.

# Output
- `SparseMatrixCSC`: The normalised, coarsened matrix `M`.
"""
function build_coarse_adjacency(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})
    k = length(clusters)
    M = spzeros(k, k)
    for i in 1:k, j in (i+1):k
        M[i, j] = length(Ad[clusters[i], clusters[j]].nzval)
    end
    for i in 1:k
        M[i, i] = sum(M[i, :])
    end
    M .+= M'
    column_normalise!(M)
    return M
end

###############################################################
###### Spectral gap detection
###############################################################

"""
    identify_spectral_gap(eigenvalues::AbstractVector{<:Number}; start::Int = 3)

Returns the number of clusters suggested by the largest relative gap among the sorted Laplacian eigenvalues.

# Arguments
- `eigenvalues::AbstractVector{<:Number}`: List of sorted eigenvalues.

# Keywords
- `start::Int`: Starting index for gap evaluation.

# Output
- `Int`: Estimated optimal cluster count `k`.
"""
function identify_spectral_gap(eigenvalues::AbstractVector{<:Number}; start::Int = 3)
    n = length(eigenvalues)
    n < start && return 1
    best_idx = start
    best_gap = -Inf
    for i in start:n
        gap = eigenvalues[i] / eigenvalues[i-1]
        if gap > best_gap
            best_gap = gap
            best_idx = i
        end
    end
    return best_idx
end

###############################################################
###### Spectral k-means
###############################################################

"""
    spectral_kmeans(lap_mat::SparseMatrixCSC, k::Int; max_iter::Int = 300)

Spectral k-means mapping via Lloyd's algorithm on Laplacian eigenvectors.

# Arguments
- `lap_mat::SparseMatrixCSC`: The Laplacian matrix.
- `k::Int`: Number of target clusters.

# Keywords
- `max_iter::Int`: Maximum iterations to process.

# Output
- `Vector{Vector{Int}}`: Resulting clustered vertex assignments.
"""
function spectral_kmeans(lap_mat::SparseMatrixCSC, k::Int; max_iter::Int = 300)
    n = size(lap_mat, 1)
    V = eigvecs(Matrix(lap_mat))[:, 1:k]

    centroids = [V[Int(ceil(t)), :] for t in LinRange(1, n, k)]
    assignments = zeros(Int, n)
    distances = Vector{Float64}(undef, k)

    for _ in 1:max_iter
        changed = false
        for i in 1:n
            for c in 1:k
                distances[c] = norm(V[i, :] .- centroids[c])
            end
            best = argmin(distances)
            if assignments[i] != best
                assignments[i] = best
                changed = true
            end
        end
        !changed && break
        for c in 1:k; centroids[c] = zeros(k); end
        counts = zeros(Int, k)
        for i in 1:n
            c = assignments[i]
            centroids[c] .+= V[i, :]
            counts[c] += 1
        end
        for c in 1:k
            counts[c] > 0 && (centroids[c] ./= counts[c])
        end
    end

    clusters = [Int[] for _ in 1:k]
    for i in 1:n
        push!(clusters[assignments[i]], i)
    end
    return clusters
end

###############################################################
###### Fluid Communities (FluidC)
###############################################################

"""
    FluidC_clustering(Ad::SparseMatrixCSC, max_k::Int; max_iters::Int = 350)

Extracts graph clusters based on the Fluid Communities algorithm by Parés et al. (2017).

# Arguments
- `Ad::SparseMatrixCSC`: Unweighted adjacency matrix.
- `max_k::Int`: Maximum number of communities.

# Keywords
- `max_iters::Int`: Termination threshold for evaluation cycles.

# Output
- `Vector{Vector{Int}}`: Processed community clusters.
"""
function FluidC_clustering(Ad::SparseMatrixCSC, max_k::Int, max_iters::Int=350)
    n = size(Ad, 1)
    seeds = sample(1:n, max_k; replace=false)

    clusters = [Int[v] for v in seeds]
    vertex_to_cluster = zeros(Int, n)
    for (c, v) in enumerate(seeds)
        vertex_to_cluster[v] = c
    end

    prev_clusters = [Int[] for _ in 1:max_k]
    pot = Vector{Float64}(undef, max_k)

    for iter in 1:max_iters
        prev_clusters .= deepcopy(clusters)

        for v in sample(1:n, n; replace=false)
            inv_sizes = 1.0 ./ length.(clusters)
            cur = vertex_to_cluster[v]
            neighbours = Ad[:, v].nzind

            for c in 1:max_k
                count_in_c = count(w -> vertex_to_cluster[w] == c, neighbours)
                cur == c && (count_in_c += 1)
                pot[c] = count_in_c * inv_sizes[c]
            end

            max_pot = maximum(pot)
            max_pot == 0 && continue
            candidates = findall(==(max_pot), pot)

            if cur != 0 && candidates != [cur]
                has_nbr = any(w -> vertex_to_cluster[w] == cur, neighbours)
                !has_nbr && filter!(!=(cur), candidates)
            end

            isempty(candidates) && continue
            cur in candidates && continue

            new_c = rand(candidates)
            cur != 0 && filter!(!=(v), clusters[cur])
            push!(clusters[new_c], v)
            vertex_to_cluster[v] = new_c
        end

        clusters == prev_clusters && break
    end
    return clusters
end

"""
    FluidC_clustering_modified(Ad::SparseMatrixCSC, max_k::Int; min_vertices::Int = 5, max_iters::Int = 350)

Modified implementation of the Fluid Communities algorithm incorporating sub-graph division criteria.

# Arguments
- `Ad::SparseMatrixCSC`: Density adjacency matrix.
- `max_k::Int`: Target maximum community count.

# Keywords
- `min_vertices::Int`: Restricts splits resulting in elements possessing fewer valid assignments.
- `max_iters::Int`: Termination threshold to prevent infinite loops.

# Output
- `Vector{Vector{Int}}`: Community lists.
"""
function FluidC_clustering_modified(Ad::SparseMatrixCSC, max_k::Int; min_vertices::Int = 5, max_iters::Int = 350)
    num_vertices = size(Ad, 1)
    new_cluster = [[v] for v in sample(1:num_vertices, max_k, replace=false)]
    vertex_list = zeros(Int, num_vertices)
    for cluster_ind in eachindex(new_cluster)
        vertex_list[new_cluster[cluster_ind][1]] = cluster_ind
    end
    cluster_volumes = ones(max_k)
    last_cluster = [Int[] for _ in 1:max_k]
    iters = 0
    cont = true

    # Add iteration limit check to the while condition
    while cont && iters < max_iters 
        last_cluster = deepcopy(new_cluster)
        
        for v in sample(1:num_vertices, num_vertices, replace=false)
            dens_cluster = cluster_volumes ./ length.(new_cluster)
            current_cluster = findfirst(x -> v in x, new_cluster)

            neighbours_v = Ad[:, v].nzind
            pot_cluster = zeros(length(new_cluster))
            for c in eachindex(pot_cluster)
                pot_cluster[c] = sum([w in new_cluster[c] ? dens_cluster[c] / (Ad + I)[v, w] : 0 for w in [neighbours_v..., v]])
            end
            C_prime = findall(x -> x == maximum(pot_cluster) && x != 0, pot_cluster)
            isempty(C_prime) && continue

            if !(current_cluster in C_prime) 
                new_cluster_ind = rand(C_prime)
                !isnothing(current_cluster) && filter!(!=(v), new_cluster[current_cluster]) 
                push!(new_cluster[new_cluster_ind], v)
                vertex_list[v] = new_cluster_ind
            end
        end
        
        dens_cluster = cluster_volumes ./ length.(new_cluster)
        non_split_clusters = deepcopy(new_cluster)
        for (cluster_ind, cluster) in enumerate(non_split_clusters)
            isempty(cluster) || length(cluster) < min_vertices && continue
            sub_adj = Ad[cluster, cluster]
            sub_components = components(sub_adj)[2]
            proper_indeces = [cluster[sub_co] for sub_co in sub_components]
            for sub_co in proper_indeces[2:end]
                push!(cluster_volumes, cluster_volumes[cluster_ind] * length(sub_co) / length(cluster))
                cluster_volumes[cluster_ind] *= (1.0 - length(sub_co) / length(cluster))
                setdiff!(cluster, sub_co)
                push!(new_cluster, sub_co)
                new_cluster_ind = length(new_cluster)
                for v in sub_co
                    vertex_list[v] = new_cluster_ind
                end
            end
        end
        
        del_list = ones(Bool, length(new_cluster))
        for cluster_ind in eachindex(new_cluster)
            if cluster_ind <= length(last_cluster)
                if isempty(new_cluster[cluster_ind])
                    del_list[cluster_ind] = false
                    continue
                end
                if new_cluster[cluster_ind] == last_cluster[cluster_ind] && ((length(new_cluster[cluster_ind]) < min_vertices) || cluster_volumes[cluster_ind] < min_vertices / num_vertices)
                    adjacent_clusters = Int[]
                    for v in new_cluster[cluster_ind]
                        for w in Ad[:, v].nzind
                            pot_clust = findfirst(x -> w in x, new_cluster)
                            if !isnothing(pot_clust)
                                del_list[pot_clust] && push!(adjacent_clusters, pot_clust)
                            end
                        end
                    end

                    adjacent_clusters = filter(!=(cluster_ind), unique(adjacent_clusters))
                    if !isempty(adjacent_clusters)
                        new_c = rand(adjacent_clusters)
                        del_list[cluster_ind] = false
                        new_cluster[new_c] = [new_cluster[new_c]; new_cluster[cluster_ind]]
                        cluster_volumes[new_c] += cluster_volumes[cluster_ind]
                    end
                end
            end
        end
        
        cluster_volumes = cluster_volumes[del_list]
        new_cluster = new_cluster[del_list]
        cont = false
        if length(new_cluster) == length(last_cluster)
            for i in eachindex(new_cluster)
                (!isempty(setdiff(new_cluster[i], last_cluster[i])) || !isempty(setdiff(last_cluster[i], new_cluster[i]))) && (cont = true)
            end
        end
        
        iters += 1
    end
    
    # Optional: Warn if the loop terminated because it hit the iteration limit
    if iters == max_iters && cont
        @warn "FluidC_clustering_modified reached max_iters ($max_iters) without fully converging. Output may be sub-optimal."
    end

    return new_cluster
end


###############################################################
###### MCL
###############################################################

"""
    MCL(Ad::SparseMatrixCSC; r::Int = 2, max_iter::Int = 500)

Executes the Markov Cluster Algorithm (MCL) on a given graph adjacency matrix to identify clusters.

The algorithm simulates random walks (expansion) followed by inflation to strengthen strong currents and weaken weak ones, eventually converging to a set of distinct clusters.

# Arguments
- `Ad::SparseMatrixCSC`: The unweighted adjacency matrix of the graph.

# Keywords
- `r::Int`: The inflation power parameter. Higher values lead to finer, more granular clusters. Defaults to 2.
- `max_iter::Int`: The maximum number of iterations allowed for convergence. Defaults to 500.

# Output
- `Tuple{Vector{Vector{Int}}, Vector{Int}}`: A tuple containing a vector of cluster arrays (each array holding the indices of vertices in that cluster) and a vector of node labels.
"""
function MCL(Ad::SparseMatrixCSC; r::Int = 2, max_iter::Int = 500)
    n = size(Ad, 1)
    A = Ad + spdiagm(vec(sum(Ad; dims=2))) / (2n)
    M = A * spdiagm(1.0 ./ vec(sum(A; dims=2)))
    M_old = spzeros(n, n)

    for iter in 1:max_iter
        M_old .= M
        M = M^2
        mcl_inflate_prune!(M, r)
        
        if M == M_old
            @debug "MCL converged after \$iter iterations."
            break
        end
    end
    return extract_clusters(M)
end

# ─────────────────────────────────────────────────────────────
# 5. Regularised MCL
# ─────────────────────────────────────────────────────────────

"""
    R_MCL(Ad::SparseMatrixCSC; r::Int = 2, tol::Float64 = 1e-2, max_iter::Int = 500)

Executes the Regularised Markov Cluster Algorithm (R-MCL) to address potential over-fragmentation.

This variant incorporates a regularisation step that penalises overly small clusters, leading to more balanced community detection by mixing the expanded matrix with the original transition probabilities.

# Arguments
- `Ad::SparseMatrixCSC`: The unweighted adjacency matrix of the graph.

# Keywords
- `r::Int`: The inflation power parameter. Defaults to 2.
- `tol::Float64`: The tolerance limit used to evaluate matrix convergence. Defaults to 1e-2.
- `max_iter::Int`: The maximum number of iterations allowed. Defaults to 500.

# Output
- `Tuple{Vector{Vector{Int}}, Vector{Int}}`: A tuple containing the identified clusters and their corresponding node labels.
"""
function R_MCL(Ad::SparseMatrixCSC; r::Int = 2, tol::Float64 = 1e-2, max_iter::Int = 500)
    n = size(Ad, 1)
    A2 = Ad + spdiagm(vec(sum(Ad; dims=2))) / (2n)
    inv_d = spdiagm(1.0 ./ vec(sum(A2; dims=2)))
    A = inv_d * A2 + A2 * inv_d
    M_G = A * spdiagm(1.0 ./ vec(sum(A; dims=2)))
    M = copy(M_G)
    M_old = spzeros(n, n)
    M_old2 = spzeros(n, n)

    for iter in 1:max_iter
        M_old2 .= M_old
        M_old .= M
        M = M * M_G
        mcl_inflate_prune!(M, r)
        column_normalise!(M)
        
        if norm(M - M_old) ≤ tol || norm(M - M_old2) ≤ tol
            @debug "R-MCL converged after \$iter iterations."
            break
        end
    end
    return extract_clusters(M)
end

"""
    restarted_MCL(Ad::SparseMatrixCSC; r::Int = 2, max_iter::Int = 350, max_clusters::Int = size(Ad, 1))

Executes a restarted, hierarchical variant of the Markov Cluster Algorithm.

This approach iteratively coarsens the graph structure by treating detected clusters as super-nodes in a new adjacency matrix, applying MCL recursively until a target maximum number of clusters is reached.

# Arguments
- `Ad::SparseMatrixCSC`: The starting adjacency matrix of the graph.

# Keywords
- `r::Int`: The inflation power parameter. Defaults to 2.
- `max_iter::Int`: The maximum number of hierarchical coarsening levels/restarts. Defaults to 350.
- `max_clusters::Int`: The convergence condition representing the maximum acceptable number of clusters. Defaults to the total number of vertices in `Ad`.

# Output
- `Vector{Vector{Int}}`: A collection of identified clusters mapped back to the original vertex indices.
"""
function restarted_MCL(Ad::SparseMatrixCSC; r::Int = 2, max_iter::Int = 350, max_clusters::Int = size(Ad, 1))
    n = size(Ad, 1)
    A = Ad + spdiagm(vec(sum(Ad; dims=2))) / (2n)
    M = A * spdiagm(1.0 ./ vec(sum(A; dims=2)))
    M_old = spzeros(size(M)...)
    clusters = [Int[i] for i in 1:n]

    for iter in 1:max_iter
        for _ in 1:500
            M_old .= M
            M = M^2
            mcl_inflate_prune!(M, r)
            M == M_old && break
        end

        level_clusters, _ = extract_clusters(M)
        new_clusters = [Int[] for _ in level_clusters]
        for (j, lc) in enumerate(level_clusters)
            for idx in lc
                append!(new_clusters[j], clusters[idx])
            end
        end
        clusters = new_clusters
        
        @debug "Restarted MCL level \$iter: reduced to \$(length(clusters)) clusters."
        length(clusters) ≤ max_clusters && break

        M = build_coarse_adjacency(Ad, clusters)
        M_old = spzeros(size(M)...)
    end
    return clusters
end

"""
    restarted_R_MCL(Ad::SparseMatrixCSC; r::Int = 2, tol::Float64 = 1e-2, max_iter::Int = 350, max_clusters::Int = size(Ad, 1))

Executes a restarted, hierarchical variant of the Regularised Markov Cluster Algorithm.

Combines the coarsening approach of restarted MCL with the regularisation logic of R-MCL to prevent over-fragmentation during the hierarchical reduction process.

# Arguments
- `Ad::SparseMatrixCSC`: The starting adjacency matrix of the graph.

# Keywords
- `r::Int`: The inflation power parameter. Defaults to 2.
- `tol::Float64`: The tolerance limit used to evaluate matrix convergence. Defaults to 1e-2.
- `max_iter::Int`: The maximum number of hierarchical coarsening levels/restarts. Defaults to 350.
- `max_clusters::Int`: The convergence condition representing the maximum acceptable number of clusters. Defaults to the total number of vertices in `Ad`.

# Output
- `Vector{Vector{Int}}`: A collection of identified clusters mapped back to the original vertex indices.
"""
function restarted_R_MCL(Ad::SparseMatrixCSC; r::Int = 2, tol::Float64 = 1e-2, max_iter::Int = 350, max_clusters::Int = size(Ad, 1))
    n = size(Ad, 1)
    A2 = Ad + spdiagm(vec(sum(Ad; dims=2))) / (2n)
    inv_d = spdiagm(1.0 ./ vec(sum(A2; dims=2)))
    A = inv_d * A2 + A2 * inv_d
    M_G = A * spdiagm(1.0 ./ vec(sum(A; dims=2)))
    M = copy(M_G)
    M_old = spzeros(size(M)...)
    M_old2 = spzeros(size(M)...)
    clusters = [Int[i] for i in 1:n]

    for iter in 1:max_iter
        for _ in 1:500
            M_old2 .= M_old
            M_old .= M
            M = M * M_G
            mcl_inflate_prune!(M, r)
            column_normalise!(M)
            norm(M - M_old) ≤ tol && break
            norm(M - M_old2) ≤ tol && break
        end

        level_clusters, _ = extract_clusters(M)
        new_clusters = [Int[] for _ in level_clusters]
        for (j, lc) in enumerate(level_clusters)
            for idx in lc
                append!(new_clusters[j], clusters[idx])
            end
        end
        clusters = new_clusters
        
        @debug "Restarted R-MCL level \$iter: reduced to \$(length(clusters)) clusters."
        length(clusters) ≤ max_clusters && break

        M = build_coarse_adjacency(Ad, clusters)
        M_G = copy(M)
        M_old = spzeros(size(M)...)
        M_old2 = spzeros(size(M)...)
    end
    return clusters
end

###############################################################
###### Recursive Normalised Cut
###############################################################

"""
    normalized_cut(Ad::SparseMatrixCSC, kmax::Int; eig_smooth_tol::Real = 0.05)

Sequentially bisects the graph utilizing recursive Fiedler partitions to identify components.

# Arguments
- `Ad::SparseMatrixCSC`: Unweighted adjacency matrix.
- `kmax::Int`: Terminating budget count representing maximum generated clusters.

# Keywords
- `eig_smooth_tol::Real`: Tolerance cutoff required for further partitioning a sub-graph.

# Output
- `Vector{Vector{Int}}`: Resulting clustered vertex assignments.
"""
function normalized_cut(Ad::SparseMatrixCSC, kmax::Int; eig_smooth_tol::Real = 0.05)
    Ad_work = copy(Ad)
    (maximum(Ad_work) > 1) && (Ad_work.nzval .= 1.0 ./ Ad_work.nzval)
    Dinv = spdiagm(vec(1.0 ./ sqrt.(sum(Ad_work; dims=1))))
    L_sym = I - Dinv * Ad_work * Dinv

    clusters = [collect(1:size(Ad_work, 1))]
    ncut_vals = [1.0]

    for _ in 1:kmax
        length(clusters) ≥ kmax && break
        new_clusters = Vector{Int}[]
        new_ncut = Float64[]

        budget = min(kmax - length(clusters), length(clusters))
        split_order = sortperm(ncut_vals; rev=true)
        to_split = Set(split_order[1:budget])

        for (i, clust) in enumerate(clusters)
            if i ∉ to_split || length(clust) < 2
                push!(new_clusters, clust)
                push!(new_ncut, ncut_vals[i])
                continue
            end

            part_a, part_b, nval = _ncut_bisect(L_sym, clust, eig_smooth_tol)
            if isnothing(part_a)
                push!(new_clusters, clust)
                push!(new_ncut, ncut_vals[i])
            else
                push!(new_clusters, part_a, part_b)
                push!(new_ncut, nval, nval)
            end
        end

        new_clusters == clusters && break
        clusters = new_clusters
        ncut_vals = new_ncut
    end
    return clusters
end

"""
    _ncut_bisect(L_sym::SparseMatrixCSC, indices::Vector{Int}, smooth_tol::Real)

Calculates optimal sparse sub-array eigenvectors executing a single partition bisection based on Fiedler approximation logic.

# Arguments
- `L_sym::SparseMatrixCSC`: Sub-Laplacian to process.
- `indices::Vector{Int}`: Global vertex representations of mapped Laplacian input scope.
- `smooth_tol::Real`: Required numeric variance margin necessary to authorize continuous cut splitting.

# Output
- `Tuple`: Tuple containing divided index partition paths `(part_a, part_b)` and best evaluated trial indicator. Return format evaluates to `(nothing, nothing, nothing)` upon boundary trigger.
"""
function _ncut_bisect(L_sym::SparseMatrixCSC, indices::Vector{Int}, smooth_tol::Real)
    Mat = L_sym[indices, indices]
    n_mat = size(Mat, 1)
    
    # Fallback to dense solver (LAPACK) since ARPACK fails to converge on extremely small matrices.
    ev = Real.(eigvecs(Matrix(Mat))[:, 2])

    ratios = abs.(ev[1:end-1] ./ ev[2:end]) .- 1
    maximum(abs, ratios) < smooth_tol && return nothing, nothing, nothing

    perm = sortperm(ev)
    best_val, best_idx = Inf, 1
    for k in LinRange(0, length(perm), 10)[2:end-1]
        ki = round(Int, k)
        indicator = ones(length(perm))
        indicator[perm[1:ki]] .= -1
        trial = indicator' * Mat * indicator
        if trial < best_val
            best_val, best_idx = trial, ki
        end
    end
    return indices[perm[1:best_idx]], indices[perm[best_idx+1:end]], best_val
end

###############################################################
###### Cluster matching — Hungarian algorithm
###############################################################

"""
    jaccard_distance_matrix(c1::Vector{<:AbstractVector}, c2::Vector{<:AbstractVector})

Returns the normalised Jaccard distance between any two input clusters. Preferable to calculating the raw symmetric difference because overall cluster dimensions will not skew analytical representations.

# Arguments
- `c1::Vector{<:AbstractVector}`: Vector of input reference clusters.
- `c2::Vector{<:AbstractVector}`: Evaluation target clusters.

# Output
- `Matrix{Float64}`: Symmetrical relational distance mappings.
"""
function jaccard_distance_matrix(c1::Vector{<:AbstractVector}, c2::Vector{<:AbstractVector})
    s1 = [Set(c) for c in c1]
    s2 = [Set(c) for c in c2]
    D  = zeros(length(s1), length(s2))
    for i in eachindex(s1), j in eachindex(s2)
        u = length(s1[i]) + length(s2[j]) - length(intersect(s1[i], s2[j]))
        D[i, j] = u == 0 ? 0.0 : 1.0 - length(intersect(s1[i], s2[j])) / u
    end
    return D
end

"""
    hungarian_match(cost::Matrix{Float64})

Generates minimum bipartite matching coordinates evaluating `O(n^3)` dimensions without external Munkres package dependencies. Internally pads matrices to square parameters guaranteeing non-quadradic functional availability.

# Arguments
- `cost::Matrix{Float64}`: Numeric representation of relationship constraints.

# Output
- `Vector{Int}`: The vector representing corresponding assignments (where `row_assignment[i] = j` identifies Row `i` assigned to Column `j`).
"""
function hungarian_match(cost::Matrix{Float64})
    m, n = size(cost)
    sz   = max(m, n)
    C    = fill(2.0, sz, sz)
    C[1:m, 1:n] .= cost

    u   = zeros(sz + 1)
    v   = zeros(sz + 1)
    p   = zeros(Int, sz + 1)
    way = zeros(Int, sz + 1)

    for i in 1:sz
        p[0+1] = i
        j0     = 0
        minVal = fill(Inf, sz + 1)
        used   = falses(sz + 1)
        while true
            used[j0 + 1] = true
            i0 = p[j0 + 1]
            δ  = Inf
            j1 = -1
            for j in 1:sz
                used[j + 1] && continue
                cur = C[i0, j] - u[i0] - v[j + 1]
                if cur < minVal[j + 1]
                    minVal[j + 1] = cur
                    way[j + 1]    = j0
                end
                if minVal[j + 1] < δ
                    δ  = minVal[j + 1]
                    j1 = j
                end
            end
            for j in 0:sz
                if used[j + 1]
                    u[p[j + 1]] += δ
                    v[j + 1]    -= δ
                else
                    minVal[j + 1] -= δ
                end
            end
            j0 = j1
            p[j0 + 1] == 0 && break
        end
        while j0 != 0
            p[j0 + 1] = p[way[j0 + 1] + 1]
            j0 = way[j0 + 1]
        end
    end

    assignment = zeros(Int, sz)
    for j in 1:sz
        p[j + 1] != 0 && (assignment[p[j + 1]] = j)
    end
    return assignment[1:m]
end

"""
    match_clusters(reference::Vector{<:AbstractVector}, target::Vector{<:AbstractVector})

Executes optimal cluster configuration routing between historical stages and present parameters ensuring visual consistency. Unmatched terminal paths are retained as `Int[]` mapping objects.

# Arguments
- `reference::Vector{<:AbstractVector}`: Baseline parameter arrays.
- `target::Vector{<:AbstractVector}`: Incoming mapped elements undergoing re-alignment.

# Output
- `Vector{Vector{Int}}`: Resulting clustered vertex array containing structured path associations extending over `max(n_ref, n_tgt)`.
"""
function match_clusters(reference::Vector{<:AbstractVector}, target::Vector{<:AbstractVector})
    isempty(reference) || isempty(target) && return target

    D          = jaccard_distance_matrix(reference, target)
    assignment = hungarian_match(D)

    n_ref = length(reference)
    n_tgt = length(target)
    out   = [Int[] for _ in 1:max(n_ref, n_tgt)]

    used = falses(n_tgt)
    for i in 1:n_ref
        j = assignment[i]
        if 1 <= j <= n_tgt
            out[i]  = target[j]
            used[j] = true
        end
    end
    
    slot = n_ref + 1
    for j in 1:n_tgt
        if !used[j] && slot <= length(out)
            out[slot] = target[j]
            slot += 1
        end
    end
    return out
end


###############################################################
###### Quality measures
###############################################################

"""
    modularity(Ad::SparseMatrixCSC, clusters::Vector)

Calculates the network modularity scalar to determine internal connection concentration thresholds separating assigned domains.

# Arguments
- `Ad::SparseMatrixCSC`: Adjacency representation object.
- `clusters::Vector`: Grouped collection allocations.

# Output
- `Float64`: Final computed fractional strength parameter.
"""
function modularity(Ad::SparseMatrixCSC, clusters::Vector)
    k = length(clusters)
    twom = sum(Ad.nzval)
    twom == 0 && return 0.0
    E = zeros(k, k)
    for i in 1:k, j in 1:k
        E[i, j] = sum(Ad[clusters[i], clusters[j]].nzval) / twom
    end
    return tr(E) - sum(sum(E; dims=1) .^ 2)
end

"""
    silhouette_index(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})

Calculates the maximum silhouette index for a given graph clustering.

The silhouette index assesses the separation distance between the resulting clusters. A higher value indicates that nodes are well matched to their own cluster and poorly matched to neighbouring clusters. It utilizes shortest path distances between nodes.

# Arguments
- `Ad::SparseMatrixCSC`: The adjacency matrix of the graph.
- `clusters::Vector{Vector{Int}}`: A vector containing the clusters, where each cluster is a vector of vertex indices.

# Output
- `Float64`: The computed maximum silhouette index.
"""
function silhouette_index(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})
    n = size(Ad, 1)
    dist = zeros(Float64, n, n)
    for i in 1:n
        dist[:, i] = dijkstra_shortest_path(Ad, i)[1]
    end
    
    k = length(clusters)
    S_avg = zeros(Float64, k)
    
    for (ci, clust) in enumerate(clusters)
        length(clust) ≤ 1 && continue
        for v in clust
            a = sum(dist[v, clust]) / (length(clust) - 1)
            b = Inf
            for (cj, clust2) in enumerate(clusters)
                cj == ci && continue
                isempty(clust2) && continue
                b = min(b, sum(dist[v, clust2]) / length(clust2))
            end
            S_avg[ci] += (b - a) / max(a, b)
        end
        S_avg[ci] /= length(clust)
    end
    return maximum(S_avg)
end


"""
    conductance(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})

Calculates the average conductance of the identified clusters.

Conductance measures the fraction of total edge volume that points outside the cluster. A lower conductance indicates better-defined clusters with fewer inter-cluster edges relative to their internal volume.

# Arguments
- `Ad::SparseMatrixCSC`: The adjacency matrix of the graph.
- `clusters::Vector{Vector{Int}}`: A vector containing the clusters, where each cluster is a vector of vertex indices.

# Output
- `Float64`: The average conductance across all valid clusters.
"""
function conductance(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})
    n = size(Ad, 1)
    total = 0.0
    count = 0
    
    for clust in clusters
        isempty(clust) && continue
        mask = falses(n)
        mask[clust] .= true
        
        cut = sum(Ad[mask, .!mask])
        vol_S = sum(Ad[mask, :])
        vol_Sc = sum(Ad[.!mask, :])
        
        denom = min(vol_S, vol_Sc)
        denom > 0.0 && (total += cut / denom)
        count += 1
    end
    
    return count > 0 ? total / count : 0.0
end

"""
    coverage(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})

Calculates the coverage metric for a given graph clustering.

Coverage is defined as the ratio of the weight of all intra-cluster edges to the total weight of all edges in the graph. Higher coverage indicates that more edges are contained within clusters rather than between them.

# Arguments
- `Ad::SparseMatrixCSC`: The adjacency matrix of the graph.
- `clusters::Vector{Vector{Int}}`: A vector containing the clusters, where each cluster is a vector of vertex indices.

# Output
- `Float64`: The computed coverage ratio.
"""
function coverage(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})
    tw = sum(Ad)
    tw == 0.0 && return 0.0
    return sum(sum(Ad[c, c]) for c in clusters) / tw
end


"""
    performance(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})

Calculates the performance metric of a graph clustering.

Performance evaluates the fraction of correctly classified pairs of nodes (nodes in the same cluster connected by an edge, and nodes in different clusters not connected by an edge) relative to the total number of possible node pairs.

# Arguments
- `Ad::SparseMatrixCSC`: The adjacency matrix of the graph.
- `clusters::Vector{Vector{Int}}`: A vector containing the clusters, where each cluster is a vector of vertex indices.

# Output
- `Float64`: The computed performance metric.
"""
function performance(Ad::SparseMatrixCSC, clusters::Vector{Vector{Int}})
    n = size(Ad, 1)
    avg_weight = sum(Ad) / length(Ad.nzval)
    I_ad, J_ad, _ = findnz(Ad)
    
    B = sparse(I_ad, J_ad, ones(Bool, length(I_ad)), n, n)
    f = sum(sum(Ad[c, c]) / 2.0 for c in clusters)
    g = 0.0
    
    for (i, c1) in enumerate(clusters)
        for c2 in clusters[i+1:end]
            g += avg_weight * sum(.!B[c1, c2])
        end
    end
    
    return 2.0 * (f + g) / (n * (n - 1))
end


"""
    ClusterMethod

Descriptor for a clustering method.

# Examples
```julia
# 1. Standard methods (default)
methods = [
    ClusterMethod("Spectral k-means", (Ad_w, adj, lap, k, nv) -> spectral_kmeans(lap, k)),
    ClusterMethod("FluidC", (Ad_w, adj, lap, k, nv) -> FluidC_clustering(adj, k)),
    ClusterMethod("GNN", (Ad_w, adj, lap, k, nv) -> iterative_deep_clustering(adj, W -> strength_aggregation_adj(W, τ=0.005), k=k))
]

# 2. Advanced/Markov-based methods (requires MCL or Restarted variants)
advanced_methods = [
    ClusterMethod("MCL", (Ad_w, adj, lap, k, nv) -> MCL(Ad_w)[1]),
    ClusterMethod("R_MCL", (Ad_w, adj, lap, k, nv) -> R_MCL(Ad_w)[1]),
    ClusterMethod("restarted_MCL", (Ad_w, adj, lap, k, nv) -> restarted_MCL(Ad_w; max_clusters=floor(Int, 0.05*nv))),
    ClusterMethod("restarted_R_MCL", (Ad_w, adj, lap, k, nv) -> restarted_R_MCL(Ad_w; max_clusters=floor(Int, 0.05*nv)))
]

# 3. Geometric/Spectral cut methods
cut_methods = [
    ClusterMethod("Normalized Cut", (Ad_w, adj, lap, k, nv) -> normalized_cut(Ad_w, k)),
    ClusterMethod("DRA_QC_CE", (Ad_w, adj, lap, k, nv) -> DRA_QC_CE_adj(adj))
]
$(FIELDS)
"""
struct ClusterMethod
    "Name of the clustering method."
    name::String
    "Callable that takes `(adj_weighted, adj_unweighted, lap_mat, k, num_verts)` and returns `Vector{Vector{Int}}`."
    run::Function
end


"""
    archaeological_clusterings(toml_file::String; num_clusters::Int = 0, max_clusters_param::Int = 10, plot_map::Bool = false, seed::Any = nothing, incl_sources::Bool = true, double_distance::Bool = false, output_dir::String = pwd(), overwrite::Bool = false, cluster_methods::Vector{ClusterMethod} = [...], quality_measures::Vector{Tuple{String, Function}} = [...])

Determines and visually parses distinct cluster combinations structuring archaeological data maps evaluating phase transitions associated with geographical material distributions.

# Arguments
- `toml_file::String`: Source structural dataset filename containing time-series elements.

# Keywords
- `num_clusters::Int`: Defined limit overriding automated gap detection algorithms. Defaults to 0 (auto-detect).
- `max_clusters_param::Int`: Maximum limit of clusters to process for spectral gap detection. Defaults to 10.
- `plot_map::Bool`: Enables the generation of geographical map plots via GMT. Defaults to false.
- `seed::Any`: Randomization seed for procedural reproducibility. Defaults to nothing.
- `incl_sources::Bool`: Triggers source site integration into the clustered models. Defaults to true.
- `double_distance::Bool`: Flag indicating extended spatial limits evaluation. Defaults to false.
- `output_dir::String`: Target export folder. Defaults to the active user working location (`pwd()`).
- `overwrite::Bool`: Enables file deletion privileges permitting dynamic historical file replacements. Defaults to false.
- `cluster_methods::Vector{ClusterMethod}`: An array of executable clustering algorithms to evaluate. See `ClusterMethod` for possible methods.
- `quality_measures::Vector{Tuple{String, Function}}`: An array of tuples linking metric names to their respective callable scoring functions.

# Output
- `Nothing`: Executes physical system disk commits recording graphical arrays, heatmaps, and diagnostic tables.
"""
function archaeological_clusterings(toml_file::String;
                                    num_clusters::Int = 0,
                                    max_clusters_param::Int = 10,
                                    plot_map::Bool = false,
                                    seed::Any = nothing,
                                    incl_sources::Bool = true,
                                    double_distance::Bool = false,
                                    output_dir::String = pwd(),
                                    overwrite::Bool = false,
                                    cluster_methods::Vector{ClusterMethod} = [
                                        ClusterMethod("Spectral k-means", (Ad_w, adj, lap, k, nv) -> spectral_kmeans(lap, k)),
                                        ClusterMethod("FluidC", (Ad_w, adj, lap, k, nv) -> FluidC_clustering(adj, k)),
                                        ClusterMethod("GNN Clustering", (Ad_w, adj, lap, k, nv) -> iterative_deep_clustering(adj, W -> strength_aggregation_adj(W, τ=0.005), k=k))
                                    ],
                                    quality_measures::Vector{Tuple{String, Function}} = [
                                        ("modularity", modularity),
                                        ("silhouette_index", silhouette_index),
                                        ("conductance", conductance),
                                        ("coverage", coverage),
                                        ("performance", performance)
                                    ])
    data_path = isfile(toml_file) ? toml_file : joinpath(dirname(@__DIR__), "DataFolder", toml_file)
    out_dir = joinpath(output_dir, string(Dates.today()), "Clusterings")
    if isdir(out_dir)
        if overwrite == false
            error("Folder '$out_dir' already exists! To overwrite files, use the argument `overwrite=true`.")
        else
            @info "Folder already exists. Files will be overwritten (overwrite=true)."
        end
    else
        mkpath(out_dir)
    end
    out_dir = output_dir

    seed !== nothing && Random.seed!(seed)

    ArData, MaData = read_data_TOML(data_path)
    
    graphen = if incl_sources
        SoData = read_source_data(data_path, MaData)
        compare_all_connect_nodes_sources(ArData, MaData, SoData, double_distance=double_distance, overwrite=overwrite, output_dir=out_dir)
    else
        compare_all_connect_nodes(ArData, MaData)
    end

    # Copy the configuration file to the output directory and append global parameters
    target_toml = joinpath(out_dir, toml_file)
    cp(data_path, target_toml, force=true)
    open(target_toml, "a") do io
        write(io, "\nmaximum_distance = $maximum_distance\nmaterial_multiplier=$material_multiplier")
    end

    method_names  = [m.name for m in cluster_methods]
    quality_names = [q[1] for q in quality_measures]
    n_methods     = length(cluster_methods)
    n_quality     = length(quality_measures)

    cluster_dict = Dict{Tuple{String,String}, Tuple{Matrix{Int}, Vector{Vector{Vector{Int}}}}}()
    quality_dict = Dict{String, Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}}()

    for graph_entry in graphen
        adj_mats     = graph_entry[1]
        time_period  = graph_entry[2]
        distance_mat = graph_entry[5]
        num_verts    = size(adj_mats[1], 1)
        n_stages     = length(adj_mats)

        @info "Processing time period: $time_period"

        for m in cluster_methods
            cluster_dict[(time_period, m.name)] = (
                zeros(Int, n_stages, num_verts),
                [Vector{Int}[] for _ in 1:n_stages]
            )
        end
        
        quality_dict[time_period] = [
            (zeros(n_quality, n_methods), zeros(n_quality, n_methods))
            for _ in 1:n_stages
        ]

        k = if num_clusters == 0
            eigenvalue_counts = [
                identify_spectral_gap(
                    eigvals(Matrix(get_matrices(adj_mats[j])[1]))[1:min(max_clusters_param, size(adj_mats[j], 1))]
                )
                for j in eachindex(adj_mats)
            ]
            round(Int, sum(eigenvalue_counts) / n_stages)
        else
            num_clusters
        end

        max_k_seen = 0

        # Tracks the last non-skipped stage to ensure valid matching even if intermediate stages are identical
        last_valid_j = Dict(m.name => 0 for m in cluster_methods)

        for (j_mat, adj_mat) in enumerate(adj_mats)
            @info "Evaluating material stage $j_mat / $n_stages"

            # Skip duplicate stages but copy clustering references to maintain indicator matrix consistency
            if j_mat == 2 && isapprox(adj_mats[1], adj_mats[2])
                for m in cluster_methods
                    cluster_dict[(time_period, m.name)][2][2] = cluster_dict[(time_period, m.name)][2][1]
                    cluster_dict[(time_period, m.name)][1][2, :] = cluster_dict[(time_period, m.name)][1][1, :]
                end
                continue
            end

            lap_mat, _ = get_matrices(adj_mat)
            Ad_w       = adj_mat .* distance_mat
            Ad_w.nzval .= 1.0 ./ Ad_w.nzval

            # ── Clustering ──
            for m in cluster_methods
                raw = m.run(Ad_w, adj_mat, lap_mat, k, num_verts)

                prev_j = last_valid_j[m.name]
                matched = if prev_j == 0
                    raw
                else
                    prev = cluster_dict[(time_period, m.name)][2][prev_j]
                    match_clusters(prev, raw)
                end

                cluster_dict[(time_period, m.name)][2][j_mat] = matched
                last_valid_j[m.name] = j_mat
            end

            # ── Indicator Matrix ──
            for m in cluster_methods
                for (c_ind, cluster) in enumerate(cluster_dict[(time_period, m.name)][2][j_mat])
                    for v in cluster
                        cluster_dict[(time_period, m.name)][1][j_mat, v] = c_ind
                    end
                end
            end

            # ── Quality Metrics ──
            nW_quality = zeros(n_quality, n_methods)
            W_quality  = zeros(n_quality, n_methods)
            for (fi, m) in enumerate(cluster_methods)
                cl = cluster_dict[(time_period, m.name)][2][j_mat]
                for (qi, (_, qfunc)) in enumerate(quality_measures)
                    nW_quality[qi, fi] = qfunc(adj_mat, cl)
                    W_quality[qi, fi]  = qfunc(Ad_w,    cl)
                end
            end
            quality_dict[time_period][j_mat] = (nW_quality, W_quality)

            max_k_seen = max(max_k_seen,
                maximum(length(cluster_dict[(time_period, m.name)][2][j_mat])
                        for m in cluster_methods))
        end

        # ── Heatmaps ──
        ref_method_name = cluster_methods[1].name
        node_order      = reduce(vcat, cluster_dict[(time_period, ref_method_name)][2][1]; init=Int[])
        unassigned      = setdiff(1:num_verts, node_order)
        node_order      = vcat(node_order, unassigned)

        my_palette = palette(:tab20, max_k_seen)
        heat_vec = [
            heatmap(
                cluster_dict[(time_period, m.name)][1][:, node_order];
                title         = "$(time_period) – $(m.name)",
                color         = my_palette,
                colorbar      = false,
                bottom_margin = 10Plots.mm,
                left_margin   = 5Plots.mm,
                xlabel        = "Sites (sorted via initial cluster)",
                ylabel        = "Material Stages",
                yticks        = (1:4:n_stages, 0:3:n_stages-1)
            )
            for m in cluster_methods
        ]
        
        heatmap_obj = plot(heat_vec...; layout=(n_methods, 1), size=(1200, n_methods * 300))
        savefig(heatmap_obj, joinpath(out_dir, "$(time_period)_heatmaps.png"))

        # ── Quality Line Plots ──
        quality_plots = []
        for (qi, (qname, _)) in enumerate(quality_measures)
            q_data = zeros(n_stages, n_methods)
            for j in 1:n_stages
                q_data[j, :] = quality_dict[time_period][j][2][qi, :]
            end
            push!(quality_plots, plot(
                1:n_stages, q_data;
                label     = reshape(method_names, 1, :),
                title     = "$qname Timeline ($time_period)",
                xlabel    = "Material Stage",
                ylabel    = qname,
                linewidth = 2,
                legend    = :outertopright,
                margin    = 5Plots.mm
            ))
        end
        if length(quality_plots) >= 2
            qual_obj = plot(quality_plots[1], quality_plots[2]; layout=(2, 1), size=(1000, 800))
            savefig(qual_obj, joinpath(out_dir, "$(time_period)_quality_metrics_plot.png"))
        end
        
        # ── Quality Metrics Text Output ──
        @show joinpath(out_dir, "$(time_period)_quality_metrics_tables.txt")
        open(joinpath(out_dir, "$(time_period)_quality_metrics_tables.txt"), "w") do io
            for j_mat in [1, n_stages]
                nW, W = quality_dict[time_period][j_mat]
                write(io, "\n========================================\n")
                write(io, "Clustering Quality – $time_period\n")
                write(io, "Material Stage $j_mat (Unweighted)\n")
                write(io, "========================================\n")
                pretty_table(io, nW; header=method_names, row_labels=quality_names, crop=:none)
                
                write(io, "\n========================================\n")
                write(io, "Clustering Quality – $time_period\n")
                write(io, "Material Stage $j_mat (Weighted)\n")
                write(io, "========================================\n")
                pretty_table(io, W; header=method_names, row_labels=quality_names, crop=:none)
            end
        end
        @info "Successfully exported metric tables and heatmaps to $out_dir"
    end

    # ── GMT Map Plotting ──
    if plot_map
        for m in cluster_methods
            @info "Initiating GMT Map rendering for clustering method: $(m.name)"
            plot_clusters(graphen,
                [cluster_dict[(graph_entry[2], m.name)][2] for graph_entry in graphen],
                m.name, output_dir=out_dir)
        end
    end
    
    return nothing
end


"""
    plot_clusters(graphen::AbstractArray, cluster_vecs::AbstractArray, cluster_method::String; plot_paths::Bool = false, output_dir::String = pwd(), overwrite::Bool = false)

Plots clusters of archaeological sites and sources onto a geographic map using GMT.

# Arguments
- `graphen::AbstractArray`: List of graph entries containing site data and adjacency matrices.
- `cluster_vecs::AbstractArray`: List of cluster assignments corresponding to the graphs.
- `cluster_method::String`: Name of the clustering method used (utilized for titles and filenames).

# Keywords
- `plot_paths::Bool`: If `true`, lines representing the edges between nodes are drawn. Defaults to `false`.
- `output_dir::String`: Base directory to save the plots. Defaults to the current working directory.
- `overwrite::Bool`: If `false` (default), throws an error if the directory exists.

# Output
- `Nothing`: Saves the generated maps as PDF and PNG files in the specified output directory.
"""
function plot_clusters(graphen::AbstractArray, cluster_vecs::AbstractArray, cluster_method::String; plot_paths::Bool = false, output_dir::String = pwd(), overwrite::Bool = false)
    
    out_dir = output_dir
    if !isdir(out_dir)
        @error "Output directory $out_dir does not exist. Please run `archaeological_clusterings` first to generate the necessary data."
        return
    end
    for i in eachindex(graphen)
        graph_entry = graphen[i]
        time_period = graph_entry[2]

        relevant_sites = graph_entry[3]
        all_nodes = if length(graph_entry) >= 4 && !isempty(graph_entry[4]) && hasproperty(first(graph_entry[4]), :coords)
            vcat(relevant_sites, graph_entry[4])
        else
            relevant_sites
        end

        long = Float64[s.coords[1] for s in all_nodes]
        lat  = Float64[s.coords[2] for s in all_nodes]

        left, right = extrema(long)
        bottom, top = extrema(lat)
        Δlong = 2.0
        Δlat  = 1.0

        adj_mats = graph_entry[1]
        n_stages = min(length(cluster_vecs[i]), length(adj_mats))
        
        for mat_ind in 1:n_stages
            clusters_at_stage = cluster_vecs[i][mat_ind]
            n_clust = count(!isempty, clusters_at_stage)
            
        GMT.coast(
            region="$(left - Δlong)/$(bottom - Δlat)/$(right + 2.0 * Δlong)/$(top + Δlat)+r",
            proj=(name=:laea, center=[10, 50]), t=75,
            frame=:none, res=:low, area=500, shore=:thin, borders=1, land="195/209/182",
            water="deepskyblue3@50"#, land="green"
        )

        # Step 2: Overlay the opaque frame and title on top of the map
        GMT.basemap!(
            frame=:ag,
            title="$(cluster_method), $time_period, Mats: $(mat_ind - 1), #Clusters: $n_clust", 
            show=false, fmt=:pdf
        )
            
            n_clust_total = length(clusters_at_stage)
            raw_colors = distinguishable_colors(n_clust_total, [colorant"white"]; dropseed=true)
            colors = cgrad(raw_colors; categorical=true)
            
            for (ci, clust) in enumerate(clusters_at_stage)
                isempty(clust) && continue

                valid_idx = Int[v for v in clust if 1 <= v <= length(all_nodes)]
                isempty(valid_idx) && continue

                col = colors[ci]
                GMT.scatter!(long[valid_idx], lat[valid_idx];
                    marker=:circle, size=0.15,
                    markerfacecolor=[col.r col.g col.b], show=false
                )
            end

            if plot_paths
                adj_mat = adj_mats[mat_ind]
                rows = rowvals(adj_mat)
                for v in axes(adj_mat, 2)
                    v > length(all_nodes) && continue
                    for idx in nzrange(adj_mat, v)
                        w = rows[idx]
                        w > length(all_nodes) && continue

                        GMT.lines!(
                            [all_nodes[v].coords[1] all_nodes[v].coords[2];
                             all_nodes[w].coords[1] all_nodes[w].coords[2]];
                            pen=(0.2, "black@50")
                        )
                    end
                end
            end

            GMT.showfig(figname=joinpath(output_dir, "$(cluster_method)_$(time_period)_mat_$(mat_ind - 1).png"))
        end
    end
    
    @info "Cluster maps successfully generated and saved."
    return nothing
end


"""
    GraphData

Data structure for holding graph attributes.
$(FIELDS)
"""
struct GraphData
    n::Int                              # Anzahl Knoten
    edge_src::Vector{Int}               # Quellknoten jeder Kante
    edge_dst::Vector{Int}               # Zielknoten jeder Kante
    node_attr::Matrix{Float64}          # [n, n_v]   Knotenattribute
    edge_attr::Matrix{Float64}          # [E, n_e]   Kantenattribute
    global_attr::Vector{Float64}        # [n_g]      Globale Attribute
end


"""
    GraphData(Lp::SparseMatrixCSC)

Constructor for `GraphData` from a Laplacian matrix `Lp`.
"""
function GraphData(Lp::SparseMatrixCSC)
    n = size(Lp, 1)
    rows = rowvals(Lp)
    vals = nonzeros(Lp)

    src = Int[]; dst = Int[]; eattr = Float64[]
    for j in 1:n
        for idx in nzrange(Lp, j)
            i = rows[idx]
            i == j && continue
            push!(src, j)
            push!(dst, i)
            push!(eattr, vals[idx])
        end
    end

    node_attr = reshape(Float64.(diag(Lp)), n, 1)
    edge_attr = reshape(eattr, length(eattr), 1)
    global_attr = Float64[]

    return GraphData(n, src, dst, node_attr, edge_attr, global_attr)
end


"""
    GNNLayer

Encapsulates GNN feature update and aggregation functions.

# Fields
$(FIELDS)
"""
struct GNNLayer
    φ_e     # Edge Update       or nothing
    ρ_ev    # Edge→Vertex Agg.  or nothing
    φ_v     # Vertex Update     or nothing
    ρ_eg    # Edge→Global Agg.  or nothing
    ρ_vg    # Vertex→Global Agg. or nothing
    φ_g     # Global Update     or nothing
end


"""
    apply_layer(layer::GNNLayer, data::GraphData)

Applies a `GNNLayer` to the provided `GraphData` structure, updating features 
for edges, nodes, and globals.

# Arguments
- `layer::GNNLayer`: The GNN layer configuration.
- `data::GraphData`: The current graph data state.

# Output
- `GraphData`: The updated `GraphData` object.
"""
function apply_layer(layer::GNNLayer, data::GraphData)
    n = data.n
    E = size(data.edge_attr, 1)

    # Edge Update 
    if layer.φ_e !== nothing
        e_test = layer.φ_e(
            data.edge_attr[1, :],
            data.node_attr[data.edge_src[1], :],
            data.node_attr[data.edge_dst[1], :],
            data.global_attr
        )
        n_e_new = length(e_test)
        new_edge_attr = zeros(E, n_e_new)
        for k in 1:E
            new_edge_attr[k, :] = layer.φ_e(
                data.edge_attr[k, :],
                data.node_attr[data.edge_src[k], :],
                data.node_attr[data.edge_dst[k], :],
                data.global_attr
            )
        end
    else
        new_edge_attr = copy(data.edge_attr)
    end


    if layer.ρ_ev !== nothing
        incoming = [Vector{Float64}[] for _ in 1:n]
        for k in 1:E
            push!(incoming[data.edge_dst[k]], new_edge_attr[k, :])
        end

        # Determine output dimension
        first_nonempty = findfirst(!isempty, incoming)
        n_agg = first_nonempty !== nothing ? length(layer.ρ_ev(incoming[first_nonempty])) : 0
        agg_attr = zeros(n, n_agg)
        for k in 1:n
            if !isempty(incoming[k])
                agg_attr[k, :] = layer.ρ_ev(incoming[k])
            end
        end
    else
        agg_attr = zeros(n, 0)
    end

    if layer.φ_v !== nothing
        ē_dummy = size(agg_attr, 2) > 0 ? agg_attr[1, :] : Float64[]
        v_test = layer.φ_v(data.node_attr[1, :], ē_dummy, data.global_attr)
        n_v_new = length(v_test)
        new_node_attr = zeros(n, n_v_new)
        for k in 1:n
            ē_k = size(agg_attr, 2) > 0 ? agg_attr[k, :] : Float64[]
            new_node_attr[k, :] = layer.φ_v(data.node_attr[k, :], ē_k, data.global_attr)
        end
    else
        new_node_attr = copy(data.node_attr)
    end


    if layer.φ_g !== nothing
        ρ_eg_res = layer.ρ_eg !== nothing ? layer.ρ_eg(new_edge_attr) : Float64[]
        ρ_vg_res = layer.ρ_vg !== nothing ? layer.ρ_vg(new_node_attr) : Float64[]
        new_global = layer.φ_g(data.global_attr, ρ_eg_res, ρ_vg_res)
    else
        new_global = copy(data.global_attr)
    end

    return GraphData(n, data.edge_src, data.edge_dst,
                     new_node_attr, new_edge_attr, new_global)
end


"""
    strength_aggregation_adj(A::SparseMatrixCSC; τ::Float64 = 0.0, max_agg_size::Int = 0)

Aggregates the adjacency matrix `A` based on connection strengths to form clusters.

# Arguments
- `A::SparseMatrixCSC`: Adjacency matrix.

# Keywords
- `τ::Float64`: Threshold for edge inclusion.
- `max_agg_size::Int`: Maximum number of nodes per cluster.

# Output
- `Vector{Vector{Int}}`: Calculated graph partitions.
"""
function strength_aggregation_adj(A::SparseMatrixCSC; τ::Float64=0.0, max_agg_size::Int=0)
    n = size(A, 1)

    # 1. Calculate degrees (row sums of the adjacency matrix)
    degrees = vec(sum(A, dims=2))

    # 2. Build the Laplacian matrix L = D - A
    # D is a diagonal matrix containing the node degrees
    Lp = spdiagm(0 => degrees) - A

    # 3. Pass the Laplacian directly to the GNN (No Graphs.jl needed!)
    S = strength_of_connection_gnn(Lp; τ=τ)

    # 4. Standard clustering procedure based on the S matrix
    s_rows = rowvals(S); s_vals = nonzeros(S)
    strong_neighbors = [Tuple{Int,Float64}[] for _ in 1:n]
    for j in 1:n
        for idx in nzrange(S, j)
            i = s_rows[idx]
            i == j && continue
            s_vals[idx] > 0 || continue
            push!(strong_neighbors[i], (j, s_vals[idx]))
        end
    end
    
    # Sort neighbors by connection strength (descending)
    for i in 1:n
        sort!(strong_neighbors[i]; by=x -> x[2], rev=true)
    end

    # Process nodes with the lowest degree first
    order = sortperm(degrees)

    aggregated = falses(n)
    Partition  = Vector{Vector{Int}}()
    
    # Build the clusters
    for seed in order
        aggregated[seed] && continue
        agg = [seed]; aggregated[seed] = true
        for (nb, _) in strong_neighbors[seed]
            aggregated[nb] && continue
            max_agg_size > 0 && length(agg) ≥ max_agg_size && break
            push!(agg, nb); aggregated[nb] = true
        end
        push!(Partition, agg)
    end
    return Partition
end


"""
    strength_of_connection_gnn(Lp::SparseMatrixCSC; τ::Float64 = 0.0)

Computes connection strengths between nodes based on GNN-inspired edge updates 
applied to the Laplacian `Lp`.

# Arguments
- `Lp::SparseMatrixCSC`: The graph Laplacian matrix.

# Keywords
- `τ::Float64`: The strength threshold for filtering weak connections.

# Output
- `SparseMatrixCSC`: The filtered connection strength matrix.
"""
function strength_of_connection_gnn(Lp::SparseMatrixCSC; τ::Float64=0.0)
    data = GraphData(Lp)


    φ_e_strength = function(e_ij, v_src, v_dst, g)
        A_ij = e_ij[1]
        A_ii = v_dst[1]    
        A_jj = v_src[1]    
        denom = A_ii * A_jj
        s_ij = abs(denom) > 1e-15 ? A_ij^2 / denom : 0.0
        return [s_ij]
    end

    layer = GNNLayer(φ_e_strength, nothing, nothing, nothing, nothing, nothing)

    result = apply_layer(layer, data)

    E = length(data.edge_src)
    I = Int[]; J = Int[]; V = Float64[]
    for k in 1:E
        s_ij = result.edge_attr[k, 1]
        if s_ij ≥ τ
            push!(I, data.edge_dst[k])
            push!(J, data.edge_src[k])
            push!(V, s_ij)
        end
    end
    return sparse(I, J, V, data.n, data.n)
end


"""
    DRA_QC_adj(W::SparseMatrixCSC; θ::Float64 = 10.0, η::Float64 = 2.0)

Performs Density-based Recursive Aggregation (DRA) on an adjacency matrix `W` to partition the graph.

# Arguments
- `W::SparseMatrixCSC`: The adjacency matrix of the graph.

# Keywords
- `θ::Float64`: The density threshold for cluster grouping. Defaults to 10.0.
- `η::Float64`: The regularization parameter for subgroup extraction. Defaults to 2.0.

# Output
- `Vector{Vector{Int}}`: The final clustered vertex partition.
"""
function DRA_QC_adj(W::SparseMatrixCSC; θ=10.0, η=2.0)
    n = size(W, 1)
    
    deg = [length(nzrange(W, j)) for j in 1:n]
    
    d_total = vec(sum(W, dims=2))
    γ       = _compute_gamma(W, d_total, n)
    active  = copy(deg)

    Partition = Vector{Vector{Int}}()

    while maximum(active) > -1
        r = argmax(active)
        G = _dra_tentative_adj(W, r, active) 
        G = _qc_extract(G, W, γ, r, θ, η)
        active[G] .= -1
        push!(Partition, G)
    end

    return Partition
end

"""
    _dra_tentative_adj(W::SparseMatrixCSC, r::Int, active::Vector{Int})

Generates a tentative aggregate centered at root node `r` based on 1-ring and 2-ring connectivity.

# Arguments
- `W::SparseMatrixCSC`: Adjacency matrix.
- `r::Int`: The root node index.
- `active::Vector{Int}`: Indicator vector for active nodes in the pool.

# Output
- `Vector{Int}`: Tentative cluster members.
"""
function _dra_tentative_adj(W::SparseMatrixCSC, r::Int, active::Vector{Int})
    G     = Int[r]
    G_set = Set(G)
    rows  = rowvals(W)

    for idx in nzrange(W, r)
        nb = rows[idx]
        if active[nb] >= 0 && nb != r
            push!(G, nb)
            push!(G_set, nb)
        end
    end

    if length(G) <= 6
        for v in copy(G)
            for idx in nzrange(W, v)
                nb = rows[idx]
                if active[nb] >= 0 && nb ∉ G_set && nb != v
                    push!(G, nb)
                    push!(G_set, nb)
                end
            end
        end
    end
    
    return G
end


"""
    _compute_gamma(W::SparseMatrixCSC, d::Vector{Float64}, n::Int)

Computes the auxiliary values `γ` representing local connectivity potential for each node.

# Output
- `Vector{Float64}`: Gamma values for density estimation.
"""
function _compute_gamma(W, d, n)
    rows = rowvals(W)
    vals = nonzeros(W)

    # Internal Sum: s_k = Σ_{i<k, (k,i)∈E} w_ki
    s = zeros(n)
    for k in 1:n
        acc = 0.0
        for idx in nzrange(W, k)
            i = rows[idx]
            i < k && (acc += vals[idx])
        end
        s[k] = acc
    end

    # γ_j = Σ_{k>j, (j,k)∈E} (w_jk / d_k) · s_k
    γ = zeros(n)
    for j in 1:n
        for idx in nzrange(W, j)
            k = rows[idx]
            w = vals[idx]
            if k > j && d[k] != 0.0
                γ[j] += (w / d[k]) * s[k]
            end
        end
    end
    return γ
end


"""
    _qc_extract(G::Vector{Int}, W::SparseMatrixCSC, γ::Vector{Float64}, r::Int, θ::Float64, η::Float64)

Performs Quality Control (QC) and refinement on a tentative cluster aggregate `G`.

# Output
- `Vector{Int}`: Refined cluster partition.
"""
function _qc_extract(G, W, γ, r, θ, η)
    length(G) == 1 && return G

    current_η = η

    while length(G) > 1

        # 1. Bad Vertices Removal
        G = _bad_vertices_removal(G, W, γ, r, θ, current_η)
        length(G) == 1 && return G

        # 2. Cheap criterion (16)
        _criterion_16(G, W, γ, r, θ) && return G

        # 3. Large aggregates: only criterion (16) → accept
        length(G) > 1024 && return G

        # 4. Cholesky-quality test
        passed, v = _cholesky_quality_test(G, W, γ, θ)
        passed && return G

        # 5. Subgroup Extraction
        current_η += 0.5
        G_neu = _subgroup_extraction(G, W, γ, r, θ, current_η, v)

        if length(G_neu) < length(G)
            G = G_neu
        else
            G_fb = [j for j in G if j == r || W[j, r] > 0.0]
            return _bad_vertices_removal(G_fb, W, γ, r, θ, current_η)
        end
    end

    return G
end


"""
    _bad_vertices_removal(G::Vector{Int}, W::SparseMatrixCSC, γ::Vector{Float64}, r::Int, θ::Float64, η::Float64)

Iteratively removes vertices from aggregate `G` that violate the density/QC conditions.

# Output
- `Vector{Int}`: Cleaned cluster partition.
"""
function _bad_vertices_removal(G, W, γ, r, θ, η)
    G     = copy(G)
    large = length(G) > 1024

    cleaned = false
    while !cleaned
        cleaned = true
        G_set   = Set(G)            
        newG    = Int[]
        sizehint!(newG, length(G))

        for j in G
            if j == r
                push!(newG, j); continue
            end

            ext  = _ext_weight(W, j, G_set)
            w_jr = W[j, r]

            c18 = (w_jr > 0.0) && ((2ext + γ[j]) / w_jr <= θ - 1)

            c19 = false
            if !c18 && !large
                int_w = _int_weight(W, j, G_set)
                c19   = (int_w > 0.0) && ((2ext + γ[j]) / int_w <= (θ - 1) / η)
            end

            (c18 || c19) ? push!(newG, j) : (cleaned = false)
        end
        G = newG
    end
    return G
end


"""
    _ext_weight(W::SparseMatrixCSC, j::Int, G_set::Set{Int})

Calculates external edge weights of vertex `j` relative to cluster set `G_set`.

# Output
- `Float64`: Computed external weight.
"""
function _ext_weight(W, j, G_set)
    ext  = 0.0
    rows = rowvals(W)
    vals = nonzeros(W)
    for idx in nzrange(W, j)
        k = rows[idx]
        k ∉ G_set && (ext += vals[idx])
    end
    return ext
end


"""
    _int_weight(W::SparseMatrixCSC, j::Int, G_set::Set{Int})

Calculates internal edge weights of vertex `j` relative to cluster set `G_set`.

# Output
- `Float64`: Computed internal weight.
"""
function _int_weight(W, j, G_set)
    s    = 0.0
    rows = rowvals(W)
    vals = nonzeros(W)
    for idx in nzrange(W, j)
        k = rows[idx]
        (k ∈ G_set && k ≠ j) && (s += vals[idx])
    end
    return s
end


"""
    _criterion_16(G::Vector{Int}, W::SparseMatrixCSC, γ::Vector{Float64}, r::Int, θ::Float64)

Evaluates the sufficiency criterion (16) for cluster stability.

# Output
- `Bool`: Returns `true` if criterion is met.
"""
function _criterion_16(G, W, γ, r, θ)
    G_set = Set(G)
    for j in G
        j == r && continue
        w_jr = W[j, r]
        w_jr == 0.0 && return false
        ext = _ext_weight(W, j, G_set)
        (2ext + γ[j]) / w_jr > θ - 1 && return false
    end
    return true
end


"""
    _cholesky_quality_test(G::Vector{Int}, W::SparseMatrixCSC, γ::Vector{Float64}, θ::Float64)

Tests cluster `G` quality using a Cholesky decomposition of the aggregated density matrix.

# Output
- `Tuple{Bool, Union{Vector{Float64}, Nothing}}`: Success status and Fiedler-like eigenvector if failed.
"""
function _cholesky_quality_test(G, W, γ, θ)
    m = length(G)
    m == 1 && return true, nothing

    AG, XG = _build_AG_XG(G, W, γ)

    # Z_G = θ·A_G - X_G + X_G·1·(1ᵀX_G1)⁻¹·1ᵀX_G
    xg1 = XG * ones(m)
    s   = dot(ones(m), xg1)
    ZG  = θ .* AG .- XG .+ (xg1 * xg1') ./ s

    # Numerical symmetrization + test of the (m-1)×(m-1) principal submatrix
    ZG_sub = (ZG[1:m-1, 1:m-1] .+ ZG[1:m-1, 1:m-1]') ./ 2
    F      = cholesky(Symmetric(ZG_sub); check=false)
    issuccess(F) && return true, nothing

    try
        AG_s = (AG .+ AG') ./ 2
        XG_s = (XG .+ XG') ./ 2
        e    = eigen(Symmetric(AG_s), Symmetric(XG_s))
        ord  = sortperm(real.(e.values))
        v    = real.(e.vectors[:, ord[2]])
        return false, v
    catch
        return false, randn(m)
    end
end


"""
    _build_AG_XG(G::Vector{Int}, W::SparseMatrixCSC, γ::Vector{Float64})

Constructs the auxiliary Laplacian `AG` and the density matrix `XG` for spectral testing.

# Output
- `Tuple{Matrix{Float64}, Matrix{Float64}}`: Tuple `(AG, XG)`.
"""
function _build_AG_XG(G, W, γ)
    m       = length(G)
    idx_map = Dict(G[i] => i for i in 1:m)
    AG      = zeros(m, m)
    ext_vec = zeros(m)
    G_set   = Set(G)
    rows    = rowvals(W)
    vals    = nonzeros(W)

    for (li, vi) in enumerate(G)
        for idx in nzrange(W, vi)
            k = rows[idx]
            w = vals[idx]
            if k ∈ G_set && k ≠ vi
                lk = idx_map[k]
                AG[li, lk]  = -w
                AG[li, li] +=  w
            else
                ext_vec[li] += w
            end
        end
    end

    Gamma_diag = [γ[G[i]] + 2.0 * ext_vec[i] for i in 1:m]
    XG         = AG + Diagonal(Gamma_diag)
    return AG, XG
end


"""
    _subgroup_extraction(G::Vector{Int}, W::SparseMatrixCSC, γ::Vector{Float64}, r::Int, θ::Float64, η::Float64, v::Vector{Float64})

Refines cluster `G` by splitting it into subgroups based on spectral information `v`.

# Output
- `Vector{Int}`: Extracted subgroup.
"""
function _subgroup_extraction(G, W, γ, r, θ, η, v)
    r_local = findfirst(==(r), G)
    vr      = v[r_local]

    Gp_set = if vr == 0.0
        Set(G[i] for i in 1:length(G) if v[i] >= 0.0)
    else
        Set(G[i] for i in 1:length(G) if v[i] * vr >= 0.0)
    end

    new_G = Int[r]
    for j in G
        j == r && continue
        ext   = _ext_weight(W, j, Gp_set)
        w_jr  = W[j, r]
        int_w = _int_weight(W, j, Gp_set)
        c_a   = (w_jr  > 0.0) && ((2ext + γ[j]) / w_jr  <= θ - 1)
        c_b   = (int_w > 0.0) && ((2ext + γ[j]) / int_w <= (θ - 1) / η)
        (c_a || c_b) && push!(new_G, j)
    end
    return new_G
end

"""
    DRA_QC_CE_adj(W::SparseMatrixCSC; θ::Float64=10.0, η::Float64=2.0)

Public wrapper for DRA with Quality Control and Complexity Enhancement.

# Output
- `Vector{Vector{Int}}`: Clustered partition.
"""
function DRA_QC_CE_adj(W::SparseMatrixCSC; θ=10.0, η=2.0)
    Partition, _ = _dra_qc_ce_internal_adj(W; θ=θ, η=η)
    return Partition
end


"""
    DRA_QC_CE_detailed_adj(W::SparseMatrixCSC; θ::Float64=10.0, η::Float64=2.0)

Wrapper for `_dra_qc_ce_internal_adj` returning detailed partition and QC status.

# Output
- `Tuple{Vector{Vector{Int}}, Vector{Bool}}`: Detailed partitioning results.
"""
function DRA_QC_CE_detailed_adj(W::SparseMatrixCSC; θ=10.0, η=2.0)
    return _dra_qc_ce_internal_adj(W; θ=θ, η=η)
end


"""
    _dra_qc_ce_internal_adj(W::SparseMatrixCSC; θ::Float64=10.0, η::Float64=2.0)

Internal orchestrator for Density-based Recursive Aggregation (DRA) with integrated Quality Control (QC) and Complexity Enhancement (CE).

# Arguments
- `W::SparseMatrixCSC`: Adjacency matrix.

# Keywords
- `θ::Float64`: Threshold for density-based grouping.
- `η::Float64`: Regularization parameter for subgroup extraction.

# Output
- `Tuple{Vector{Vector{Int}}, Vector{Bool}}`: Tuple containing the partition and associated QC success flags.
"""
function _dra_qc_ce_internal_adj(W::SparseMatrixCSC; θ=10.0, η=2.0)
    n       = size(W, 1)
    d_total = vec(sum(W, dims=2))
    γ       = _compute_gamma(W, d_total, n)

    # Calculate node degree as the number of non-zero entries per column
    deg    = [length(nzrange(W, j)) for j in 1:n]
    active = copy(deg) 

    # Sort nodes by degree in descending order (calculated only once)
    order = sortperm(deg, rev=true)

    Partition = Vector{Vector{Int}}()
    qc_flags  = Vector{Bool}()

    # Main Loop: DRA with Quality Control
    ptr = 1 
    while true 
        # Find the next free root node
        while ptr <= n && active[order[ptr]] < 0 
            ptr += 1
        end
        ptr > n && break 
        r = order[ptr] 

        # Call the adjacency-based tentative function
        G = _dra_tentative_adj(W, r, active)
        G = _qc_extract(G, W, γ, r, θ, η)
        active[G] .= -1
        push!(Partition, G)
        push!(qc_flags, true)
    end

    # Complexity Enhancement (CE)
    n_c = length(Partition)
    if n_c > n ÷ 4
        # Identify clusters with 3 or fewer nodes
        small_idx = findall(p -> length(p) <= 3, Partition)
        
        # Safety check to ensure there are actually small clusters to process
        if !isempty(small_idx)
            # Release nodes back into the active pool
            released  = vcat(Partition[small_idx]...)
            deleteat!(Partition, sort(small_idx))
            deleteat!(qc_flags,  sort(small_idx))

            for v in released
                active[v] = deg[v]
            end

            # Second DRA pass WITHOUT Quality Control for the released nodes
            ptr = 1 
            while true 
                while ptr <= n && active[order[ptr]] < 0 
                    ptr += 1
                end
                ptr > n && break 
                r = order[ptr] 

                G = _dra_tentative_adj(W, r, active)
                active[G] .= -1
                push!(Partition, G)
                push!(qc_flags, false)
            end
        end
    end

    return Partition, qc_flags
end

"""
    build_coarse_graph(W_fine::SparseMatrixCSC{Float64, Int}, Partition::Vector{Vector{Int}})

Constructs a coarsened adjacency matrix from a fine graph and a cluster partition.

# Arguments
- `W_fine::SparseMatrixCSC{Float64, Int}`: Original fine-grained adjacency matrix.
- `Partition::Vector{Vector{Int}}`: Current clustering partition.

# Output
- `SparseMatrixCSC`: Coarsened adjacency matrix.
"""
function build_coarse_graph(W_fine::SparseMatrixCSC{Float64, Int}, Partition::Vector{Vector{Int}})
    n_fine = size(W_fine, 1)
    n_coarse = length(Partition)
    I = Int[]; J = Int[]; V = Float64[]
    for (coarse_node, cluster) in enumerate(Partition)
        for fine_node in cluster
            push!(I, fine_node); push!(J, coarse_node); push!(V, 1.0)
        end
    end
    P = sparse(I, J, V, n_fine, n_coarse)
    W_coarse = P' * W_fine * P
    W_coarse = W_coarse - spdiagm(0 => diag(W_coarse))
    dropzeros!(W_coarse)
    return W_coarse
end

"""
    iterative_deep_clustering(W_initial::SparseMatrixCSC, cluster_func::Function; min_reduction::Float64=0.85, max_iter::Int=10, k::Int=10)

Performs recursive hierarchical clustering by coarsening the graph structure.

This function iteratively applies a clustering function, coarsens the graph based on the resulting partitions, and repeats until the cluster count reduction stagnates or a minimum threshold is reached.

# Arguments
- `W_initial::SparseMatrixCSC`: Initial weighted adjacency matrix.
- `cluster_func::Function`: Clustering algorithm function to apply at each iteration.

# Keywords
- `min_reduction::Float64`: Threshold for cluster count reduction to trigger early termination. Defaults to 0.85.
- `max_iter::Int`: Maximum number of recursive levels. Defaults to 10.
- `k::Int`: Minimum target cluster count. Defaults to 10.

# Output
- `Vector{Vector{Int}}`: Final partitioned vertex assignments mapped back to the original graph.
"""
function iterative_deep_clustering(W_initial::SparseMatrixCSC, cluster_func::Function; 
                                   min_reduction::Float64=0.85, 
                                   max_iter::Int=10, k::Int=10)
    
    W_current = W_initial
    
    # First clustering on the original graph
    current_partition = cluster_func(W_current)
    global_partition = current_partition 
    
    println("Iter 1: $(size(W_initial, 1)) Startknoten -> $(length(global_partition)) Cluster")
    
    # rekursive Loop: Coarsen the graph and cluster again
    for iter in 2:max_iter
        n_old = length(global_partition)
        
        # Create a coarse graph based on the current partition
        W_current = build_coarse_graph(W_current, current_partition)
        
        # Cluster new, coarse graph
        current_partition = cluster_func(W_current)
        n_new = length(current_partition)
        
        # println("Iter $iter: $n_old Knoten -> $n_new Cluster")
        
        # Map the new coarse clusters back to the original base nodes
        new_global_partition = Vector{Vector{Int}}()
        for coarse_cluster in current_partition
            base_nodes = Int[]
            
            for c_idx in coarse_cluster
                append!(base_nodes, global_partition[c_idx])
            end
            push!(new_global_partition, base_nodes)
        end
        
        if n_new > n_old * min_reduction
            println("-> Stop: Number of clusters stagnates (Reduction <= $((1 - min_reduction)*100)%).")
            break
        end
        if n_new < k
            println("-> Stop: Very few clusters left (<= $k).")
            break
        end

        global_partition = new_global_partition
    end

    return global_partition
end