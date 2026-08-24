using SparseArrays
using Statistics



"""
    build_G(A::Dict{Int,<:AbstractMatrix})

Constructs the logical quadtree evaluation matrix `G` from coefficient matrices `A`.

# Arguments
- `A::Dict{Int, <:AbstractMatrix}`: Dictionary of coefficient matrices.

# Output
- `Dict{Int, BitMatrix}`: A dictionary of boolean matrices indicating uniform cells.
"""
function build_G(A::Dict{Int,<:AbstractMatrix})
    isempty(A) && error("A may not be empty.")
    levels = sort(collect(keys(A)))
    L = maximum(levels)
    
    G = Dict{Int, BitMatrix}()

    # Finest coefficient level
    G[L] = BitMatrix(A[L] .== 0)

    # Bottom-up evaluation
    for ℓ in (L-1):-1:0
        Gℓ = falses(2^ℓ, 2^ℓ)

        for i in 1:2^ℓ, j in 1:2^ℓ
            local_good = (A[ℓ][i,j] == 0)

            children_good =
                G[ℓ+1][2i-1, 2j-1] &&
                G[ℓ+1][2i-1, 2j]   &&
                G[ℓ+1][2i,   2j-1] &&
                G[ℓ+1][2i,   2j]

            Gℓ[i,j] = local_good && children_good
        end

        G[ℓ] = Gℓ
    end
    return G
end


"""
    adaptive_leaf_cells(G::Dict{Int,BitMatrix}, L_geom::Int)

Extracts the terminal leaf cells of the adaptive quadtree structure.

# Arguments
- `G::Dict{Int, BitMatrix}`: Dictionary of boolean quadtree matrices.
- `L_geom::Int`: The finest geometric level.

# Output
- `Vector{Tuple{Int,Int,Int}}`: Vector of leaf cells defined as `(level, row, col)`.
"""
function adaptive_leaf_cells(G::Dict{Int,BitMatrix}, L_geom::Int)
    L_coeff = maximum(keys(G))
    leaves = Tuple{Int,Int,Int}[]

    function descend(ℓ::Int, i::Int, j::Int)
        if ℓ <= L_coeff && G[ℓ][i,j]
            push!(leaves, (ℓ, i, j))
            return
        end
        if ℓ == L_geom
            push!(leaves, (ℓ, i, j))
            return
        end

        descend(ℓ+1, 2i-1, 2j-1)
        descend(ℓ+1, 2i-1, 2j)
        descend(ℓ+1, 2i,   2j-1)
        descend(ℓ+1, 2i,   2j)
    end

    descend(0, 1, 1)
    return sort!(leaves)
end

"""
    HatNode

A struct representing a node in the adaptive hat-function graph.

$(FIELDS)
"""
struct HatNode
    level::Int
    row::Float64 
    col::Float64
end

"""
    leaves_to_hat_nodes(leaves::Vector{Tuple{Int,Int,Int}}, L_geom::Int)

Converts quadtree leaf tuples into `HatNode` objects.

# Arguments
- `leaves::Vector{Tuple{Int,Int,Int}}`: List of quadtree leaf cells.
- `L_geom::Int`: The finest geometric grid level.

# Output
- `Vector{HatNode}`: List of initialized hat nodes.
"""
function leaves_to_hat_nodes(leaves::Vector{Tuple{Int,Int,Int}}, L_geom::Int)
    hat_nodes = HatNode[]
    for (ℓ, i, j) in leaves
        s = 2^(L_geom - ℓ)
        # Exact center of the support on the finest grid
        r_center = (i - 1) * s + (s / 2.0)
        c_center = (j - 1) * s + (s / 2.0)
        push!(hat_nodes, HatNode(ℓ, r_center, c_center))
    end
    return hat_nodes
end


"""
    hat_adjacency_matrix(hat_nodes::Vector{HatNode}, L_geom::Int)

Creates a directed, unweighted adjacency matrix for hat-function nodes using 8-connectivity.

# Arguments
- `hat_nodes::Vector{HatNode}`: List of nodes in the graph.
- `L_geom::Int`: Finest geometric resolution.

# Output
- `SparseMatrixCSC{Int, Int}`: Adjacency matrix of connected nodes.
"""
function hat_adjacency_matrix(hat_nodes::Vector{HatNode}, L_geom::Int)
    n = length(hat_nodes)
    N_fine = 2^L_geom
    
    @debug "Creating lookup table for adjacency..."
    pixel_to_node = zeros(Int, N_fine, N_fine)
    
    for idx in 1:n
        node = hat_nodes[idx]
        s = 2^(L_geom - node.level) 
        
        r_start = Int(node.row - s/2) + 1
        c_start = Int(node.col - s/2) + 1
        
        pixel_to_node[r_start:(r_start+s-1), c_start:(c_start+s-1)] .= idx
    end
    
    @debug "Scanning pixel grid for 8-connectivity..."
    I = Int[]; J = Int[]; V = Int[]
    
    for r in 1:N_fine
        for c in 1:N_fine
            u = pixel_to_node[r, c]
            u == 0 && continue
            
            # 1. Right neighbor
            if c < N_fine
                v = pixel_to_node[r, c+1]
                if u != v && v != 0
                    push!(I, u); push!(J, v); push!(I, v); push!(J, u)
                end
            end
            
            # 2. Bottom neighbor
            if r < N_fine
                v = pixel_to_node[r+1, c]
                if u != v && v != 0
                    push!(I, u); push!(J, v); push!(I, v); push!(J, u)
                end
            end
            
            # 3. Diagonal Bottom-Right
            if r < N_fine && c < N_fine
                v = pixel_to_node[r+1, c+1]
                if u != v && v != 0
                    push!(I, u); push!(J, v); push!(I, v); push!(J, u)
                end
            end
            
            # 4. Diagonal Bottom-Left
            if r < N_fine && c > 1
                v = pixel_to_node[r+1, c-1]
                if u != v && v != 0
                    push!(I, u); push!(J, v); push!(I, v); push!(J, u)
                end
            end
        end
    end
    
    @debug "Creating final sparse matrix..."
    Ad = sparse(I, J, 1, n, n)
    Ad.nzval .= 1
    
    return Ad
end


"""
    get_weighted_hat_matrix(Ad_hat01, hat_nodes, TopData, L_geom; Cost_func, raster)

Applies travel costs to the unweighted hat adjacency matrix using terrain slope.

# Arguments
- `Ad_hat01::SparseMatrixCSC{Int,Int}`: Unweighted adjacency matrix.
- `hat_nodes::Vector{HatNode}`: List of graph nodes.
- `TopData::TopographicalData`: Topographical boundaries and data.
- `L_geom::Int`: Finest geometric level.

# Keywords
- `Cost_func::Function`: Function determining the travel cost (default: `adjusted_Hiking`).
- `raster::AbstractMatrix`: Elevation raster (default: `TopData.CostMap`).

# Output
- `SparseMatrixCSC{Float64, Int}`: Weighted adjacency matrix representing travel costs.
"""
function get_weighted_hat_matrix(Ad_hat01::SparseMatrixCSC{Int,Int},
                                 hat_nodes::Vector{HatNode}, 
                                 TopData::TopographicalData,
                                 L_geom::Int;
                                 Cost_func::Function = adjusted_Hiking,
                                 raster::AbstractMatrix = TopData.CostMap)
    
    # Extract height at the center of the hat function
    function hat_center_height(node::HatNode, raster::AbstractMatrix)
        r = clamp(round(Int, node.row), 1, size(raster, 1))
        c = clamp(round(Int, node.col), 1, size(raster, 2))
        return raster[r, c]
    end

    # Penalty factor for level jumps
    αfun = (node_u, node_v) -> begin
        α0 = 0.1; β = 0.5
        α0 * 2.0^(β * (L_geom - min(node_u.level, node_v.level)))
    end

    n = length(hat_nodes)
    
    # Real geographic coordinates of the centers
    centers = [ind2coords((clamp(round(Int, node.row), 1, size(raster, 1)), 
                           clamp(round(Int, node.col), 1, size(raster, 2))), TopData) 
               for node in hat_nodes]
               
    heights = [hat_center_height(node, raster) for node in hat_nodes]

    rows, cols, _ = findnz(Ad_hat01)
    vals = Vector{Float64}(undef, length(rows))

    for t in eachindex(rows)
        v = rows[t]
        u = cols[t]

        dist_uv = distance(centers[u], centers[v])

        # Prevent self-loops or zero distances
        if dist_uv <= 0
            vals[t] = Inf
            continue
        end

        αuv = αfun(hat_nodes[u], hat_nodes[v])
        S_uv = αuv * (heights[v] - heights[u]) / dist_uv
        
        # Final assignment of costs (e.g., hiking time in hours)
        vals[t] = dist_uv * Cost_func(S_uv)
    end

    return sparse(rows, cols, vals, n, n)
end

"""
    find_nearest_hat_node(target_coord::Tuple{Float64, Float64}, hat_nodes::Vector{HatNode}, TopData::TopographicalData)

Finds the index of the node in the adaptive graph that is closest to a geographic coordinate.

# Arguments
- `target_coord::Tuple{Float64, Float64}`: Target `(latitude, longitude)`.
- `hat_nodes::Vector{HatNode}`: List of available graph nodes.
- `TopData::TopographicalData`: Map boundaries and data.

# Output
- `Int`: Index of the nearest node.
"""
function find_nearest_hat_node(target_coord::Tuple{Float64, Float64}, hat_nodes::Vector{HatNode}, TopData::TopographicalData)
    best_idx = 1
    min_dist = Inf
    
    for (i, node) in enumerate(hat_nodes)
        # Calculate real coordinates of the node
        r = clamp(round(Int, node.row), 1, size(TopData.CostMap, 1))
        c = clamp(round(Int, node.col), 1, size(TopData.CostMap, 2))
        node_coord = ind2coords((r, c), TopData)
        
        # Distance using Haversine
        d = distance(target_coord, node_coord)
        
        if d < min_dist
            min_dist = d
            best_idx = i
        end
    end
    
    return best_idx
end


@inline function wavelet_offset(N_m)
    if N_m == N_2
        return -1
    elseif N_m == N_2_boundary
        return 1
    end
    return 0
end

@inline side_from_level(L::Int, offset::Int) = 2^L + offset

function max_level_for_offset(min_side::Int, offset::Int)
    if offset == -1
        return floor(Int, log2(min_side + 1))
    elseif offset == 1
        min_side > 1 || error("For N_2_boundary, min(size(CostMap)) > 1 must hold.")
        return floor(Int, log2(min_side - 1))
    end
    return floor(Int, log2(min_side))
end

"""
    extract_site_centered_quadratic_map(TopoData, start_site, target_site, N_m; margin_ratio=0.5)

Extracts a square sub-region around two target sites.

The side length is dynamically scaled to `2^L + offset` (depending on the wavelet type)
to ensure compatibility with subsequent wavelet pipelines.

# Arguments
- `TopoData::TopographicalData`: Full topographical map.
- `start_site::Tuple`: Start coordinate.
- `target_site::Tuple`: Target coordinate.
- `N_m::N_m_and_Wavelet`: Wavelet configuration.

# Keywords
- `margin_ratio::Float64`: Extra padding around the sites (default: `0.5`).

# Output
- `Tuple{TopographicalData, Int}`: Extracted map data and the determined level `L`.
"""
function extract_site_centered_quadratic_map(TopoData::TopographicalData,
                                             start_site::Tuple{Float64,Float64},
                                             target_site::Tuple{Float64,Float64},
                                             N_m;
                                             margin_ratio::Float64=0.5)
    margin_ratio >= 0 || error("margin_ratio muss >= 0 sein.")

    m, n = size(TopoData.CostMap)
    min_side = min(m, n)
    offset = wavelet_offset(N_m)

    L_max = max_level_for_offset(min_side, offset)

    start_lexi = coords2lexi(start_site, TopoData)
    target_lexi = coords2lexi(target_site, TopoData)

    start_ij = lexi2ind(start_lexi, TopoData)
    target_ij = lexi2ind(target_lexi, TopoData)

    r1, c1 = start_ij
    r2, c2 = target_ij

    row_min, row_max = minmax(r1, r2)
    col_min, col_max = minmax(c1, c2)

    needed = max(row_max - row_min + 1, col_max - col_min + 1)
    target_span = max(needed, ceil(Int, needed * (1 + margin_ratio)))

    L = 0
    while L <= L_max && side_from_level(L, offset) < target_span
        L += 1
    end

    L <= L_max || error("Sites plus margin do not fit into a valid square cutout.")

    side = side_from_level(L, offset)

    row_lo = max(1, row_max - side + 1)
    row_hi = min(row_min, m - side + 1)
    col_lo = max(1, col_max - side + 1)
    col_hi = min(col_min, n - side + 1)

    row_lo <= row_hi || error("No valid row range found for site-centered cutout.")
    col_lo <= col_hi || error("No valid column range found for site-centered cutout.")

    center_r = (r1 + r2) / 2
    center_c = (c1 + c2) / 2

    r_start_candidate = round(Int, center_r - (side - 1) / 2)
    c_start_candidate = round(Int, center_c - (side - 1) / 2)

    r_start = clamp(r_start_candidate, row_lo, row_hi)
    c_start = clamp(c_start_candidate, col_lo, col_hi)

    r_end = r_start + side - 1
    c_end = c_start + side - 1

    cut = TopoData.CostMap[r_start:r_end, c_start:c_end]

    lBot_new = ind2coords((r_end, c_start), TopoData)
    rTop_new = ind2coords((r_start, c_end), TopoData)

    return TopographicalData(cut, size(cut), lBot_new, rTop_new, (1, 1)), L
end


"""
    plot_site_centered_cutout(TopData, TopoData_cut, start_site, target_site; show_plot, backend, save_prefix)

Generates a diagnostic visual with two panels (or GMT images):
1) Full map with marked cutout bounds.
2) Extracted cutout.

# Arguments
- `TopData::TopographicalData`: Original map.
- `TopoData_cut::TopographicalData`: Extracted sub-map.
- `start_site::Tuple{Float64,Float64}`: Start coordinates.
- `target_site::Tuple{Float64,Float64}`: Target coordinates.
"""
function plot_site_centered_cutout(TopData::TopographicalData,
                                   TopoData_cut::TopographicalData,
                                   start_site::Tuple{Float64,Float64},
                                   target_site::Tuple{Float64,Float64};
                                   show_plot::Bool=true,
                                   backend::Symbol=:auto,
                                   save_prefix::Union{Nothing,String}=nothing)
    use_gmt = backend == :gmt || (backend == :auto && @isdefined(GMT))

    if use_gmt
        reg_full = "$(TopData.lBot[1])/$(TopData.lBot[2])/$(TopData.rTop[1])/$(TopData.rTop[2])+r"
        reg_cut = "$(TopoData_cut.lBot[1])/$(TopoData_cut.lBot[2])/$(TopoData_cut.rTop[1])/$(TopoData_cut.rTop[2])+r"

        x1, y1 = TopoData_cut.lBot
        x2, y2 = TopoData_cut.rTop
        rect_xy = [
            x1 y1;
            x2 y1;
            x2 y2;
            x1 y2;
            x1 y1
        ]
        @info "Creating GMT visualization..."
        GMT.grdimage("@earth_relief_01m", region=reg_full, proj=:M15c, cmap=:geo, frame=:ag)
        GMT.lines!(rect_xy, pen=(2.0, :red))
        GMT.scatter!([start_site[1], target_site[1]], [start_site[2], target_site[2]], marker=:circle, size=0.25, markerfacecolor=:white, markeredgecolor=:black, show=true)

        if save_prefix !== nothing
            GMT.savefig("$(save_prefix)_full.png")
        end

        GMT.grdimage("@earth_relief_01m", region=reg_cut, proj=:M15c, cmap=:geo, frame=:ag)
        GMT.scatter!([start_site[1], target_site[1]], [start_site[2], target_site[2]], marker=:circle, size=0.25, markerfacecolor=:white, markeredgecolor=:black, show=true)

        if save_prefix !== nothing
            GMT.savefig("$(save_prefix)_cutout.png")
        end
        return nothing
    end

    m_full, n_full = size(TopData.CostMap)
    x_full = range(TopData.lBot[1], TopData.rTop[1], length=n_full)
    y_full = range(TopData.rTop[2], TopData.lBot[2], length=m_full)

    p_full = heatmap(
        x_full, y_full, TopData.CostMap;
        aspect_ratio=:equal,
        legend=false,
        colorbar=false,
        title="Full Map + Cutout",
        xlabel="Longitude",
        ylabel="Latitude"
    )

    x1, y1 = TopoData_cut.lBot
    x2, y2 = TopoData_cut.rTop
    plot!(p_full, [x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], color=:red, lw=2, label=false)
    scatter!(p_full, [start_site[1], target_site[1]], [start_site[2], target_site[2]], color=:white, ms=4, label=false)

    m_cut, n_cut = size(TopoData_cut.CostMap)
    x_cut = range(TopoData_cut.lBot[1], TopoData_cut.rTop[1], length=n_cut)
    y_cut = range(TopoData_cut.rTop[2], TopoData_cut.lBot[2], length=m_cut)

    p_cut = heatmap(
        x_cut, y_cut, TopoData_cut.CostMap;
        aspect_ratio=:equal,
        legend=false,
        colorbar=false,
        title="Extracted Cutout",
        xlabel="Longitude",
        ylabel="Latitude"
    )
    scatter!(p_cut, [start_site[1], target_site[1]], [start_site[2], target_site[2]], color=:white, ms=4, label=false)

    p = plot(p_full, p_cut, layout=(1, 2), size=(1200, 450))
    show_plot && display(p)
    return p
end


"""
    run_adaptive_hat_pipeline(tiff_file, N_m; r, N, start_site, target_site, margin_ratio)

Executes the complete pipeline from loading a TIFF file to creating a weighted adaptive adjacency matrix.

Optionally extracts a site-centered cutout if start and target sites are provided.
"""
function run_adaptive_hat_pipeline(tiff_file::String, N_m;
                                   r::Float64=-1.0,
                                   N::Int=-1,
                                   start_site::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                                   target_site::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                                   margin_ratio::Float64=0.5)
    A1 = tiff_to_matrix(tiff_file)
    TopData = TopographicalData(A1, size(A1), (-13,35), (46,70), (1,1))

    if start_site !== nothing && target_site !== nothing
        TopoData_quad, maxL = extract_site_centered_quadratic_map(TopData, start_site, target_site, N_m, margin_ratio=margin_ratio)
    else
        TopoData_quad, maxL = extract_quadratic_map(TopData, wavelet_offset(N_m))
    end

    L_geom = maxL
    
    # Pass keyword arguments to get_adaptive_matrices
    A, comp_rate = get_adaptive_matrices(TopoData_quad.CostMap, N_m, r=r, N=N)
    G = build_G(A)
    leaves = adaptive_leaf_cells(G, L_geom)
    
    hat_nodes = leaves_to_hat_nodes(leaves, L_geom)
    Ad_hat01 = hat_adjacency_matrix(hat_nodes, L_geom)
    Ad_weighted = get_weighted_hat_matrix(Ad_hat01, hat_nodes, TopoData_quad, L_geom)
    
    @debug "Weighted graph generated " 
    
    return TopoData_quad, Ad_weighted, hat_nodes, leaves, L_geom, comp_rate
end


"""
    run_uniform_reference_pipeline(tiff_file, L_ref, N_m; start_site, target_site, margin_ratio)

Generates a uniform reference graph rapidly by circumventing the adaptive evaluation bottlenecks.
"""
function run_uniform_reference_pipeline(tiff_file::String, L_ref::Int, N_m;
                                        start_site::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                                        target_site::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                                        margin_ratio::Float64=0.5,
                                        plot_cutout::Bool=false,
                                        cutout_plot_backend::Symbol=:auto)
   
    A1 = tiff_to_matrix(tiff_file)
    TopData = TopographicalData(A1, size(A1), (-13,35), (46,70), (1,1))

    if start_site !== nothing && target_site !== nothing
        TopoData_quad, maxL = extract_site_centered_quadratic_map(TopData, start_site, target_site, N_m, margin_ratio=margin_ratio)
        if plot_cutout
            plot_site_centered_cutout(TopData, TopoData_quad, start_site, target_site, backend=cutout_plot_backend)
        end
    else
        TopoData_quad, maxL = extract_quadratic_map(TopData, wavelet_offset(N_m))
    end

    L_geom = maxL

    if L_ref > maxL
        @warn "L_ref ($L_ref) is larger than maxL ($maxL). Enforcing L_ref = maxL."
        L_ref = maxL
    end

    @info "Extracting uniform reference grid at level $L_ref..."
    
    if L_ref < maxL
        coarse_raster = compression_2D(TopoData_quad.CostMap, L_ref, N_m=N_m, threshold=false, type=:iso)
    else
        coarse_raster = TopoData_quad.CostMap
    end

    n_cells = size(coarse_raster, 1)
    s = 2^(maxL - L_ref) 

    # Generate HatNodes
    hat_nodes = HatNode[]
    for j in 1:n_cells
        for i in 1:n_cells
            r_center = (i - 1) * s + (s / 2.0)
            c_center = (j - 1) * s + (s / 2.0)
            push!(hat_nodes, HatNode(L_ref, r_center, c_center))
        end
    end

    I = Int[]; J = Int[]; V = Float64[]
    
    α0 = 0.1; β = 0.5
    α_const = α0 * 2.0^(β * (maxL - L_ref))

    for j in 1:n_cells
        for i in 1:n_cells
            u_idx = (j - 1) * n_cells + i
            center_u = ind2coords((clamp(round(Int, hat_nodes[u_idx].row), 1, size(TopoData_quad.CostMap, 1)), 
                                   clamp(round(Int, hat_nodes[u_idx].col), 1, size(TopoData_quad.CostMap, 2))), TopoData_quad)
            h_u = coarse_raster[i, j]

            for dj in -1:1, di in -1:1
                if di == 0 && dj == 0 continue end
                
                ni = i + di
                nj = j + dj
                
                if 1 <= ni <= n_cells && 1 <= nj <= n_cells
                    v_idx = (nj - 1) * n_cells + ni
                    center_v = ind2coords((clamp(round(Int, hat_nodes[v_idx].row), 1, size(TopoData_quad.CostMap, 1)), 
                                           clamp(round(Int, hat_nodes[v_idx].col), 1, size(TopoData_quad.CostMap, 2))), TopoData_quad)
                    h_v = coarse_raster[ni, nj]

                    dist_uv = distance(center_u, center_v)
                    
                    if dist_uv > 0
                        S_uv = α_const * (h_v - h_u) / dist_uv
                        cost = dist_uv * adjusted_Hiking(S_uv)
                        
                        push!(I, u_idx)
                        push!(J, v_idx)
                        push!(V, cost)
                    end
                end
            end
        end
    end

    Ad_weighted = sparse(I, J, V, length(hat_nodes), length(hat_nodes))
    
    leaves = [(L_ref, Int(ceil(node.row/s)), Int(ceil(node.col/s))) for node in hat_nodes]
    
    return TopoData_quad, Ad_weighted, hat_nodes, leaves, L_geom, 0.0
end


### Flagged for CleanUp
"""
    example_adaptive(wavelet_type, r_tol; filepath, start_site, target_site, output_dir)

Example application demonstrating adaptive least-cost path evaluation and plotting.
"""
function example_adaptive(wavelet_type, r_tol;
                          filepath::String = joinpath(@__DIR__, "src", "DataFolder", "-13_35_46_70.tiff"),
                          start_site::Tuple{Float64, Float64} = (-5.4, 38.6),
                          target_site::Tuple{Float64, Float64} = (10.0, 45.0),
                          output_dir::String = pwd(),
                          overwrite::Bool = false)

    TopoData_quad, Ad_weighted, hat_nodes, leaves, L_geom, _ = run_adaptive_hat_pipeline(
        filepath,
        wavelet_type,
        r=r_tol,
        start_site=start_site,
        target_site=target_site
    )

    start_idx = find_nearest_hat_node(start_site, hat_nodes, TopoData_quad)
    target_idx = find_nearest_hat_node(target_site, hat_nodes, TopoData_quad)

    @info "Removing water connections..."
    keepNodes = adaptive_WaterLandBool(leaves, TopoData_quad, L_geom, [start_idx, target_idx])
    Ad_weighted = delConnections(Ad_weighted, keepNodes)
    
    labels, sizes = components(Ad_weighted)

    start_komponente = labels[start_idx]
    ziel_komponente = labels[target_idx]

    @debug "Start node is in component: $start_komponente"
    @debug "Target node is in component: $ziel_komponente"

    if start_komponente == ziel_komponente
        @info "Path is clear! Start and target nodes are connected over land."
    else
        @warn "ATTENTION: Start and target nodes are separated by water! Path computation might fail."
    end

    @info "Calculating shortest path from node $start_idx to $target_idx..."

    allocated_list, predecessor = dijkstra_shortest_path(Ad_weighted, start_idx)
    path_indices = get_path(target_idx, start_idx, predecessor)
    
    if isempty(path_indices)
        @warn "No path could be found!"
        return
    end

    reverse!(path_indices)

    path_coords = Tuple{Float64, Float64}[]
    for idx in path_indices
        node = hat_nodes[idx]
        r = clamp(round(Int, node.row), 1, size(TopoData_quad.CostMap, 1))
        c = clamp(round(Int, node.col), 1, size(TopoData_quad.CostMap, 2))
        push!(path_coords, ind2coords((r, c), TopoData_quad))
    end

    @info "LCP successfully calculated! Path length (nodes): $(length(path_coords))"
    @info "Total travel cost: $(allocated_list[target_idx])"
    
    @info "Rendering GMT map..."

    left_border = min(start_site[1], target_site[1])
    right_border = max(start_site[1], target_site[1])
    bottom_border = min(start_site[2], target_site[2])
    top_border = max(start_site[2], target_site[2])

    Δlong = 3; Δlat = 3
    reg = "$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+Δlong)/$(top_border+Δlat)+r"

    GMT.grdimage(filepath, region=reg, proj=:M15c, cmap=:geo)
    GMT.coast!(shorelines=:thin, borders=1, frame=:ag, area=500)

    for k in 1:length(path_coords)-1
        GMT.lines!([path_coords[k][1] path_coords[k][2]; path_coords[k+1][1] path_coords[k+1][2]], pen=(2.0,:white))
    end
    
    if output_dir != pwd()
        if isdir(output_dir)
            if overwrite == false
                error("Folder '$output_dir' already exists! To overwrite files, use the argument `overwrite=true`.")
            else
                @info "Folder already exists. Files will be overwritten (overwrite=true)."
            end
        else
            mkpath(output_dir)
        end
    end

    GMT.scatter!(
        [start_site[1], target_site[1]], 
        [start_site[2], target_site[2]],
        fmt=:png, 
        marker=:circle,
        markeredgecolor=0, 
        size=0.3, 
        markerfacecolor=:blue, 
        show=true,
        figname=joinpath(output_dir, "Single_Adaptive_LCP_GMT")
    )

    @info "GMT plot successfully created and saved at: $(joinpath(output_dir, "Single_Adaptive_LCP_GMT.png"))"
end

"""
    plot_lcp_gmt(filepath, path_coords, start_site, target_site, filename; output_dir)

Generates a GMT map plot for a specific generated path.
"""
function plot_lcp_gmt(filepath::String, path_coords::Vector{}, 
                      start_site::Tuple{Number, Number}, target_site::Tuple{Number, Number}, 
                      filename::String;
                      output_dir::String = pwd(),
                      overwrite::Bool = false)
    
    left_border = min(start_site[1], target_site[1])
    right_border = max(start_site[1], target_site[1])
    bottom_border = min(start_site[2], target_site[2])
    top_border = max(start_site[2], target_site[2])

    Δlong = 3; Δlat = 3
    reg = "$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+Δlong)/$(top_border+Δlat)+r"

    GMT.grdimage(filepath, region=reg, proj=:M15c, cmap=:geo, t=75)
    GMT.coast!(shorelines=:thin, borders=1, frame=:ag, area=500)

    for k in 1:length(path_coords)-1
        GMT.lines!([path_coords[k][1] path_coords[k][2]; path_coords[k+1][1] path_coords[k+1][2]], pen=(2.0,:black))
    end

    if output_dir != pwd()
        if isdir(output_dir)
            if overwrite == false
                error("Folder '$output_dir' already exists! To overwrite files, use the argument `overwrite=true`.")
            else
                @info "Folder already exists. Files will be overwritten (overwrite=true)."
            end
        else
            mkpath(output_dir)
        end
    end

    GMT.scatter!(
        [start_site[1], target_site[1]], 
        [start_site[2], target_site[2]],
        fmt=:png, 
        marker=:circle,
        markeredgecolor=0, 
        size=0.3, 
        markerfacecolor=:blue, 
        show=false, 
        figname=joinpath(output_dir, filename)
    )
end


"""
    plot_lcp_gmt_multi(filepath, path_specs, start_site, target_site, filename; output_dir)

Generates a comparative GMT plot containing multiple paths defined in `path_specs`.
"""
function plot_lcp_gmt_multi(filepath::String,
                            path_specs::Vector{},
                            start_site::Tuple{Number, Number}, target_site::Tuple{Number, Number},
                            filename::String;
                            output_dir::String = pwd(),
                            overwrite::Bool = false)

    left_border = min(start_site[1], target_site[1])
    right_border = max(start_site[1], target_site[1])
    bottom_border = min(start_site[2], target_site[2])
    top_border = max(start_site[2], target_site[2])

    Δlong = 3; Δlat = 3
    reg = "$(left_border-Δlong)/$(bottom_border-Δlat)/$(right_border+Δlong)/$(top_border+Δlat)+r"

    GMT.grdimage(filepath, region=reg, proj=:M15c, cmap=:geo, t=75)
    GMT.coast!(shorelines=:thin, borders=1, frame=:ag, area=500)

    for (path_coords, linecolor) in path_specs
        for k in 1:length(path_coords)-1
            GMT.lines!([path_coords[k][1] path_coords[k][2]; path_coords[k+1][1] path_coords[k+1][2]], pen=(2.0, linecolor))
        end
    end

    if output_dir != pwd()
        if isdir(output_dir)
            if overwrite == false
                error("Folder '$output_dir' already exists! To overwrite files, use the argument `overwrite=true`.")
            else
                @info "Folder already exists. Files will be overwritten (overwrite=true)."
            end
        else
            mkpath(output_dir)
        end
    end

    GMT.scatter!(
        [start_site[1], target_site[1]],
        [start_site[2], target_site[2]],
        fmt=:png,
        marker=:circle,
        markeredgecolor=0,
        size=0.3,
        markerfacecolor=:white,
        show=false,
        figname=joinpath(output_dir, filename)
    )
end

"""
    get_component_stats(Ad::SparseMatrixCSC, start_idx::Int)

Evaluates graph component properties.

# Output
- `Tuple{Int, Int}`: Node and edge counts of the connected component.
"""
function get_component_stats(Ad::SparseMatrixCSC, start_idx::Int)
    labels, _ = components(Ad)
    
    start_comp_id = labels[start_idx]
    vertices_in_comp = findall(==(start_comp_id), labels)
    
    V_comp = length(vertices_in_comp)
    E_comp = sum(Ad.colptr[v+1] - Ad.colptr[v] for v in vertices_in_comp) ÷ 2
    
    return V_comp, E_comp
end


"""
    benchmark_adaptive_lcpa(filepath, start_site, target_site; margin_ratio, L_ref, plot_cutout, cutout_plot_backend, wavelet_choices, r_values, N_values)

Runs a comprehensive performance and accuracy benchmark between different compression thresholds.

# Output
- `Vector{Tuple}`: Array of tracking results across different runs. Generates LaTeX table to stdout.

# Keywords
- `wavelet_choices::Vector{Symbol}`: Choose from `[:N_1]`, `[:N_2]`, or `[:N_1, :N_2]`.
"""
function benchmark_adaptive_lcpa(filepath::String, start_site::Tuple{Float64, Float64}, 
                                 target_site::Tuple{Float64, Float64}; 
                                 margin_ratio::Float64=0.5,
                                 L_ref::Int=10,
                                 plot_cutout::Bool=true,
                                 cutout_plot_backend::Symbol=:auto,
                                 wavelet_choices::Vector{Symbol}=[:N_1, :N_2],
                                 r_values::Vector{Float64}=Float64[], 
                                 N_values::Vector{Int}=Int[])
    
    # Types: (Method, Threshold_Type, Threshold_Value, Comp_Rate, Vertices, Edges, Cost, Distance, Rel_Error)
    results = Tuple{String, String, Float64, Float64, Int, Int, Float64, Float64, Float64}[]

    # Resolve and validate explicit wavelet selection.
    valid_choices = Symbol[:N_1, :N_2]
    if isempty(wavelet_choices)
        error("wavelet_choices must not be empty. Choose from :N_1 and :N_2.")
    end
    invalid_choices = [w for w in unique(wavelet_choices) if !(w in valid_choices)]
    if !isempty(invalid_choices)
        error("Unsupported wavelet choices $(invalid_choices). Allowed: :N_1 and :N_2.")
    end

    wavelet_types = [(w == :N_1 ? N_1 : N_2, String(w)) for w in unique(wavelet_choices)]
    reference_wavelet = wavelet_types[1][1]

    # ---------------------------------------------------------
    # 1. Calculate Reference Solution (uniform grid on local cutout)
    # ---------------------------------------------------------
    @info "=== Berechne Referenzlösung (uniformes lokales Subset) ==="
    TopoData_ref, Ad_ref, hat_nodes_ref, leaves_ref, L_geom_ref, _ = run_uniform_reference_pipeline(
        filepath, L_ref, reference_wavelet,
        start_site=start_site,
        target_site=target_site,
        margin_ratio=margin_ratio,
        plot_cutout=plot_cutout,
        cutout_plot_backend=cutout_plot_backend
    )

    start_idx_ref = find_nearest_hat_node(start_site, hat_nodes_ref, TopoData_ref)
    target_idx_ref = find_nearest_hat_node(target_site, hat_nodes_ref, TopoData_ref)
    
    keepNodes_ref = adaptive_WaterLandBool(leaves_ref, TopoData_ref, L_geom_ref, [start_idx_ref, target_idx_ref])
    Ad_ref = delConnections(Ad_ref, keepNodes_ref)
    
    V_ref_comp, E_ref_comp = get_component_stats(Ad_ref, start_idx_ref)
    
    cost_list_ref, pred_ref = dijkstra_shortest_path(Ad_ref, start_idx_ref)
    ref_cost = cost_list_ref[target_idx_ref]
    println("Referenzkosten: $ref_cost")
    
    path_indices_ref = get_path(target_idx_ref, start_idx_ref, pred_ref)
    reverse!(path_indices_ref)
    
    path_coords_ref = [ind2coords((clamp(round(Int, hat_nodes_ref[idx].row), 1, size(TopoData_ref.CostMap, 1)), 
                                   clamp(round(Int, hat_nodes_ref[idx].col), 1, size(TopoData_ref.CostMap, 2))), TopoData_ref) 
                       for idx in path_indices_ref]

    ref_dist = 0.0
    for k in 1:(length(path_coords_ref)-1)
        ref_dist += distance(path_coords_ref[k], path_coords_ref[k+1])
    end
    
    push!(results, ("Reference", "-", 0.0, 0.0, V_ref_comp, E_ref_comp, ref_cost, ref_dist, 0.0))
    plot_lcp_gmt(filepath, path_coords_ref, start_site, target_site, "LCP_Reference_FullGrid")

    for N in N_values
        @info "=== Berechne Adaptive LCPA für N = $N ==="

        comparison_paths = Dict{String, Vector{Tuple{Float64,Float64}}}()

        for (wt, wt_name) in wavelet_types
            @info "--- Wavelet $wt_name ---"

            TopoData_ad, Ad_adapt, hat_nodes_ad, leaves_ad, L_geom_ad, comp_rate = run_adaptive_hat_pipeline(
                filepath, wt,
                N=N,
                start_site=start_site,
                target_site=target_site,
                margin_ratio=margin_ratio
            )

            start_idx_ad = find_nearest_hat_node(start_site, hat_nodes_ad, TopoData_ad)
            target_idx_ad = find_nearest_hat_node(target_site, hat_nodes_ad, TopoData_ad)

            keepNodes_ad = adaptive_WaterLandBool(leaves_ad, TopoData_ad, L_geom_ad, [start_idx_ad, target_idx_ad])
            Ad_adapt = delConnections(Ad_adapt, keepNodes_ad)

            V_ad_comp, E_ad_comp = get_component_stats(Ad_adapt, start_idx_ad)

            cost_list_ad, pred_ad = dijkstra_shortest_path(Ad_adapt, start_idx_ad)
            ad_cost = cost_list_ad[target_idx_ad]

            path_indices_ad = get_path(target_idx_ad, start_idx_ad, pred_ad)
            reverse!(path_indices_ad)

            path_coords_ad = [ind2coords((clamp(round(Int, hat_nodes_ad[idx].row), 1, size(TopoData_ad.CostMap, 1)),
                                          clamp(round(Int, hat_nodes_ad[idx].col), 1, size(TopoData_ad.CostMap, 2))), TopoData_ad)
                              for idx in path_indices_ad]

            ad_dist = 0.0
            for k in 1:(length(path_coords_ad)-1)
                ad_dist += distance(path_coords_ad[k], path_coords_ad[k+1])
            end

            rel_error = abs(ad_cost - ref_cost) / ref_cost * 100.0

            push!(results, (wt_name, "N", Float64(N), comp_rate * 100.0, V_ad_comp, E_ad_comp, ad_cost, ad_dist, rel_error))

            filename = "LCP_Adaptive_$(wt_name)_N_$(N)"
            plot_lcp_gmt(filepath, path_coords_ad, start_site, target_site, filename)
            comparison_paths[wt_name] = path_coords_ad
        end

        if haskey(comparison_paths, "N_1") && haskey(comparison_paths, "N_2")
            plot_lcp_gmt_multi(
                filepath,
                [(path_coords_ref, :black), (comparison_paths["N_1"], :red), (comparison_paths["N_2"], :blue)],
                start_site,
                target_site,
                "LCP_Comparison_N_$(N)"
            )
        end
    end

    for r in r_values
        @info "=== Berechne Adaptive LCPA für r = $r ==="

        comparison_paths = Dict{String, Vector{Tuple{Float64,Float64}}}()

        for (wt, wt_name) in wavelet_types
            @info "--- Wavelet $wt_name ---"

            TopoData_ad, Ad_adapt, hat_nodes_ad, leaves_ad, L_geom_ad, comp_rate = run_adaptive_hat_pipeline(
                filepath, wt,
                r=r,
                start_site=start_site,
                target_site=target_site,
                margin_ratio=margin_ratio
            )

            start_idx_ad = find_nearest_hat_node(start_site, hat_nodes_ad, TopoData_ad)
            target_idx_ad = find_nearest_hat_node(target_site, hat_nodes_ad, TopoData_ad)

            keepNodes_ad = adaptive_WaterLandBool(leaves_ad, TopoData_ad, L_geom_ad, [start_idx_ad, target_idx_ad])
            Ad_adapt = delConnections(Ad_adapt, keepNodes_ad)

            V_ad_comp, E_ad_comp = get_component_stats(Ad_adapt, start_idx_ad)

            cost_list_ad, pred_ad = dijkstra_shortest_path(Ad_adapt, start_idx_ad)
            ad_cost = cost_list_ad[target_idx_ad]

            path_indices_ad = get_path(target_idx_ad, start_idx_ad, pred_ad)
            reverse!(path_indices_ad)

            path_coords_ad = [ind2coords((clamp(round(Int, hat_nodes_ad[idx].row), 1, size(TopoData_ad.CostMap, 1)),
                                          clamp(round(Int, hat_nodes_ad[idx].col), 1, size(TopoData_ad.CostMap, 2))), TopoData_ad)
                              for idx in path_indices_ad]

            ad_dist = 0.0
            for k in 1:(length(path_coords_ad)-1)
                ad_dist += distance(path_coords_ad[k], path_coords_ad[k+1])
            end

            rel_error = abs(ad_cost - ref_cost) / ref_cost * 100.0

            push!(results, (wt_name, "r", r, comp_rate * 100.0, V_ad_comp, E_ad_comp, ad_cost, ad_dist, rel_error))

            filename = "LCP_Adaptive_$(wt_name)_r_$(replace(string(r), "." => "_"))"
            plot_lcp_gmt(filepath, path_coords_ad, start_site, target_site, filename)
            comparison_paths[wt_name] = path_coords_ad
        end

        if haskey(comparison_paths, "N_1") && haskey(comparison_paths, "N_2")
            plot_lcp_gmt_multi(
                filepath,
                [(path_coords_ref, :black), (comparison_paths["N_1"], :red), (comparison_paths["N_2"], :blue)],
                start_site,
                target_site,
                "LCP_Comparison_r_$(replace(string(r), "." => "_"))"
            )
        end
    end

    printstyled("\n=== Benchmark Abgeschlossen ===\n", color=:green)
    println("\n% --- LaTeX Table Output ---")
    println("\\begin{table}[htbp]")
    println("    \\centering")
    println("    \\begin{tabular}{llrrrrrr}")
    println("        \\toprule")
    println(raw"        \textbf{Method} & \textbf{Type} & \textbf{Value} & \textbf{Comp. (\%)} & \textbf{$|V_{\text{comp}}|$} & \textbf{$|E_{\text{comp}}|$} & \textbf{Total Cost} & \textbf{Distance} & \textbf{Error (\%)} \\")
    println("        \\midrule")
    
    for (method, t_type, t_val, comp_rate, v, e, cost, dist, err) in results
        val_str = t_type == "N" ? string(Int(t_val)) : (t_type == "r" ? @sprintf("%.3f", t_val) : "-")
        @printf("        %s & %s & %s & %.2f & %d & %d & %.2f & %.2f & %.3f \\\\\n", method, t_type, val_str, comp_rate, v, e, cost, dist, err)
    end
    
    println("        \\bottomrule")
    println("    \\end{tabular}")
    println(raw"    \caption{Comparison of adaptive LCPA complexity (measured on the fully connected component), path distance, and cost accuracy.}")
    println("    \\label{tab:lcpa_benchmark}")
    println("\\end{table}")
    println("% --------------------------\n")

    return results
end

### Flagged for CleanUp

# filepath = "src\\DataFolder\\-13_35_46_70.tiff"
# start_site = (-5.4, 38.6)
# target_site = (10.0, 45.0)
# r_values = [0.01, 0.025, 0.05, 0.1, 0.2]
# N_values = [5000, 20000, 50000, 75000, 150000]
# benchmark_adaptive_lcpa(filepath, start_site, target_site, r_values=r_values, N_values=N_values,margin_ratio=0.25)


# # Alternative start- end target sites für weitere Tests
# filepath = "src\\DataFolder\\-13_35_46_70.tiff"
# start_site = (10.0,45.0)
# target_site = (14.0,50.0)
# # r_values = [0.0005, 0.01, 0.025, 0.05, 0.1, 0.2]
# # N_values = [5000, 20000, 50000, 75000, 150000, 350000]
# r_values = [0.0005]
# N_values = [5000, 150000]
# benchmark_adaptive_lcpa(filepath, start_site, target_site, r_values=r_values, N_values=N_values,margin_ratio=0.72)
