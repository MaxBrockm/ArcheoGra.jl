@inline G_to_int(G::BitMatrix) = Int.(G)

"""
    get_A_j(j::Int, Img_decomp::Matrix{Float64}, N_m)

Extracts the binary coefficient matrix for level `j` based on the wavelet decomposition.

# Arguments
- `j::Int`: The decomposition level.
- `Img_decomp::Matrix{Float64}`: The fully decomposed wavelet coefficient matrix.
- `N_m::N_m_and_Wavelet`: The chosen wavelet configuration.

# Output
- `Matrix{Bool}`: A boolean matrix indicating where wavelet coefficients are non-zero.
"""
function get_A_j(j::Int, Img_decomp::Matrix{Float64}, N_m)
    # Extract exact boundaries based on the wavelet type
    Δ_j = N_m.Delta_l(j)
    Δ_j1 = N_m.Delta_l(j+1)
    
    # The quadtree structure always expects a strictly 2^j x 2^j grid
    limit = min(Δ_j, 2^j)
    
    H = zeros(Float64, 2^j, 2^j)
    V = zeros(Float64, 2^j, 2^j)
    
    # Horizontal details
    H[1:limit, :] .= abs.(Img_decomp[1:limit, Δ_j+1:Δ_j1])
    if Δ_j > 2^j
        for extra in (2^j+1):Δ_j
            H[end, :] .+= abs.(Img_decomp[extra, Δ_j+1:Δ_j1])
        end
    end
    
    # Vertical details
    V[:, 1:limit] .= abs.(Img_decomp[Δ_j+1:Δ_j1, 1:limit])
    if Δ_j > 2^j
        for extra in (2^j+1):Δ_j
            V[:, end] .+= abs.(Img_decomp[Δ_j+1:Δ_j1, extra])
        end
    end
    
    D = abs.(Img_decomp[Δ_j+1:Δ_j1, Δ_j+1:Δ_j1])
    
    Img_coeffs = H .+ V .+ D
    A_j = Img_coeffs .!= 0
    
    return A_j
end

"""
    get_adaptive_matrices(Img::AbstractMatrix, N_m; r=-1.0, N=-1)

Generates a dictionary of adaptive coefficient matrices across all levels based on thresholding.

# Arguments
- `Img::AbstractMatrix`: The original raster image.
- `N_m::N_m_and_Wavelet`: The wavelet configuration.

# Keywords
- `r::Float64`: Relative error threshold (default: `-1.0`).
- `N::Int`: Number of largest coefficients to retain (default: `-1`).

# Output
- `Tuple{Dict{Int, Matrix{Int}}, Float64}`: Dictionary of coefficient matrices and the achieved compression rate.
"""
function get_adaptive_matrices(Img, N_m; r::Float64=-1.0, N::Int=-1)
    # Calculate exact L depending on wavelet structure
    L = 0
    if N_m == N_2
        L = Int(log2(size(Img,1) + 1))
    elseif N_m == N_2_boundary
        L = Int(log2(size(Img,1) - 1))
    else
        L = Int(log2(size(Img,1)))
    end
    
    A = Dict{Int,Matrix{Int}}()
    
    # Pass both r and N down to compression_2D
    _, Img_coeffs = compression_2D(Img, 2, N_m=N_m, threshold=true, r=r, N=N, type=:iso, adapt_ref=true)
    
    comp_rate = 1.0 - (count(!iszero, Img_coeffs) / length(Img_coeffs))

    A[L] = zeros(Int, 2^L, 2^L)
    for j in L-1:-1:0
        # Pass N_m down to properly resolve the indices
        A[j] = get_A_j(j, Img_coeffs, N_m)
    end
    return A, comp_rate
end

"""
    leaf_edge_masks(max_blocks::Vector{Tuple{Int,Int,Int}}, Lpix::Int)

Generates edge boundary masks on the finest pixel grid for a given set of quadtree blocks.

# Arguments
- `max_blocks::Vector{Tuple{Int,Int,Int}}`: Vector of block tuples `(ℓ, i, j)`.
- `Lpix::Int`: The finest pixel grid level (generates `N = 2^Lpix` pixels per dimension).

# Output
- `Tuple{BitMatrix, BitMatrix}`: A tuple `(hmask, vmask)` where:
  - `hmask`: Horizontal segments mask of size `(N+1, N)`.
  - `vmask`: Vertical segments mask of size `(N, N+1)`.
"""
function leaf_edge_masks(max_blocks::Vector{Tuple{Int,Int,Int}}, Lpix::Int)
    N = 2^Lpix

    hmask = falses(N + 1, N)  # horizontal unit-segments
    vmask = falses(N, N + 1)   # vertical unit-segments

    for (ℓ, i, j) in max_blocks
        0 <= ℓ <= Lpix || error("Block-Level ℓ=$ℓ ist nicht kompatibel mit Lpix=$Lpix")

        s = 2^(Lpix - ℓ)   # Block size in pixels
        row0 = (i - 1) * s
        row1 = i * s
        col0 = (j - 1) * s
        col1 = j * s

        # Top and bottom edges
        hmask[row0 + 1, col0 + 1 : col1] .= true
        hmask[row1 + 1, col0 + 1 : col1] .= true

        # Left and right edges
        vmask[row0 + 1 : row1, col0 + 1] .= true
        vmask[row0 + 1 : row1, col1 + 1] .= true
    end

    return hmask, vmask
end


"""
    cut_leaves(Lpix::Int, leaves::Vector{Tuple{Int,Int,Int}})

Extract a subset of leaf nodes up to a specific pixel index.

This function searches the provided vector of `leaves` for the first tuple where the first element matches `Lpix`. It then returns a new vector containing all elements prior to this match.

# Arguments
- `Lpix::Int`: The target pixel index used as the stopping criterion.
- `leaves::Vector{Tuple{Int,Int,Int}}`: A vector of tuples representing the leaf nodes, where the first integer of each tuple corresponds to a pixel index.

# Output
- `Vector{Tuple{Int,Int,Int}}`: A truncated vector containing all leaf tuples that appear strictly before the matching `Lpix`.
"""
@inline cut_leaves(Lpix::Int, leaves::Vector{Tuple{Int,Int,Int}}) = leaves[1:(findfirst(t -> t[1] == Lpix, leaves)-1)]


"""
    plot_quadtree_leaves(max_blocks::Vector{Tuple{Int,Int,Int}}; Lpix, img=nothing, lw=2, linecolor=:black)

Plots the boundaries of the quadtree leaf blocks on the finest pixel grid.

# Arguments
- `max_blocks::Vector{Tuple{Int,Int,Int}}`: Vector of quadtree blocks `(ℓ, i, j)`.

# Keywords
- `Lpix::Int`: The finest pixel grid level.
- `img::Union{Nothing, AbstractMatrix}`: Optional background image matrix of size `N×N` to render as a heatmap.
- `lw::Real`: Line width (default: 2).
- `linecolor::Symbol`: Line color (default: `:black`).

# Output
- `Plots.Plot`: The generated plot object.
"""
function plot_quadtree_leaves(max_blocks::Vector{Tuple{Int,Int,Int}};
                              Lpix::Int,
                              img=nothing,
                              lw=2,
                              linecolor=:black)

    N = 2^Lpix
    hmask, vmask = leaf_edge_masks(max_blocks, Lpix)
    
    @info "Creating plot with N=$N (Horizontal segments: $(count(x->x, hmask)), Vertical segments: $(count(x->x, vmask)))"
    
    if img === nothing
        p = plot(
            xlim=(0, N), ylim=(0, N),
            aspect_ratio=1, legend=false, grid=false, ticks=0:N, framestyle=:box, yflip=true
        )
    else
        size(img) == (N, N) || error("img must have size $(N)×$(N)")
        p = heatmap(
            0.5:1:(N-0.5), 0.5:1:(N-0.5), img;
            yflip=true, colorbar=false, legend=false, aspect_ratio=1, xlims=(0, N), ylims=(0, N), framestyle=:box
        )
    end

    # Horizontal segments
    for y in 0:N
        y % 1000 == 0 && @debug "Processing row $y"
        for x in 1:N
            if hmask[y + 1, x]
                plot!(p, [x - 1, x], [y, y], color=linecolor, lw=lw, label=false)
            end
        end
    end

    # Vertical segments
    for y in 1:N
        y % 1000 == 0 && @debug "Processing column $y"
        for x in 0:N
            if vmask[y, x + 1]
                plot!(p, [x, x], [y - 1, y], color=linecolor, lw=lw, label=false)
            end
        end
    end

    return p
end

"""
    leaf_cells_to_level_matrix(leaves::Vector{Tuple{Int,Int,Int}}, ℓ::Int)

Fills a representation matrix up to a specific level with the quadtree leaf data.

For a leaf `(ℓ′, i, j)`:
- If `ℓ′ <= ℓ`, the corresponding block is filled with the level `ℓ′`.
- If `ℓ′ > ℓ`, the block remains `0`.

# Arguments
- `leaves::Vector{Tuple{Int,Int,Int}}`: Vector of leaf cells.
- `ℓ::Int`: Target representation level.

# Output
- `Matrix{Int}`: Matrix of size `2^ℓ × 2^ℓ` containing the mapping.
"""
function leaf_cells_to_level_matrix(leaves::Vector{Tuple{Int,Int,Int}}, ℓ::Int)
    N = 2^ℓ
    M = zeros(Int, N, N)

    for (ℓ′, i, j) in leaves
        if ℓ′ > ℓ
            continue
        end

        s = 2^(ℓ - ℓ′)   # Block size on the target representation level ℓ

        r1 = (i - 1) * s + 1
        r2 = i * s
        c1 = (j - 1) * s + 1
        c2 = j * s

        M[r1:r2, c1:c2] .= ℓ′
    end

    return M
end


"""
    group_leaves_by_level(leaves::Vector{<:Tuple})

Groups a sorted list of leaves by their decomposition level.

# Output
- `Tuple{Dict{Int,Int}, Dict{Int,Vector{Tuple{Int,Int}}}}`: A tuple containing a mapping of levels to their starting indices, and a dictionary of grouped elements.
"""
function group_leaves_by_level(leaves::Vector{<:Tuple})
    level_start = Dict{Int,Int}()
    by_level = Dict{Int,Vector{Tuple{Int,Int}}}()

    current_level = nothing

    for (k, leaf) in pairs(leaves)
        ℓ, i, j = leaf

        if current_level !== ℓ
            current_level = ℓ
            level_start[ℓ] = k
            by_level[ℓ] = Tuple{Int,Int}[]
        end

        push!(by_level[ℓ], (i, j))
    end

    return level_start, by_level
end


"""
    leaf_adjacency_matrix_directed(leaves::Vector{Tuple{Int,Int,Int}}, ℓ_0::Int)

Constructs a directed unweighted adjacency matrix mapping the connectivity of adaptive leaf cells.

Nodes are connected if they share an edge or corner (King's move logic). Transitions 
between different levels (fine to coarse and vice versa) are appropriately mapped.

# Arguments
- `leaves::Vector{Tuple{Int,Int,Int}}`: Sorted list of leaf nodes.
- `ℓ_0::Int`: Coarsest level limit.

# Output
- `SparseMatrixCSC{Int, Int}`: A directed sparse adjacency matrix.
"""
function leaf_adjacency_matrix_directed(leaves::Vector{Tuple{Int,Int,Int}}, ℓ_0::Int)
    isempty(leaves) && return spzeros(Int, 0, 0)

    level_start, by_level = group_leaves_by_level(leaves)
    L = maximum(first(t) for t in leaves)

    # pos_in_level[ℓ][(i,j)] = position of (i,j) within by_level[ℓ]
    pos_in_level = Dict{Int, Dict{Tuple{Int,Int}, Int}}()
    for (ℓ, vec) in by_level
        d = Dict{Tuple{Int,Int}, Int}()
        for (p, ij) in pairs(vec)
            d[ij] = p
        end
        pos_in_level[ℓ] = d
    end

    global_index(ℓ, i, j) = level_start[ℓ] + pos_in_level[ℓ][(i,j)] - 1

    king_moves = (
        (-1, -1), (-1, 0), (-1, 1),
        ( 0, -1),           ( 0, 1),
        ( 1, -1), ( 1, 0), ( 1, 1)
    )

    n = length(leaves)
    I = Int[]
    J = Int[]
    V = Int[]

    # Traverse from fine to coarse levels
    for ℓ in L:-1:ℓ_0
        haskey(by_level, ℓ) || continue

        M = leaf_cells_to_level_matrix(leaves, ℓ)
        N = size(M, 1)

        for (i, j) in by_level[ℓ]
            idx1 = global_index(ℓ, i, j)

            for (di, dj) in king_moves
                r = i + di
                c = j + dj

                if !(1 <= r <= N && 1 <= c <= N)
                    continue
                end

                ℓnbr = M[r, c]

                # 0 => Only finer structures exist here on this representation level
                if ℓnbr == 0
                    continue
                end

                if ℓnbr == ℓ
                    idx2 = global_index(ℓ, r, c)
                else
                    # gröberes Level
                    s = 2^(ℓ - ℓnbr)
                    i2 = fld(r - 1, s) + 1
                    j2 = fld(c - 1, s) + 1
                    idx2 = global_index(ℓnbr, i2, j2)
                end

                if idx1 != idx2
                    push!(I, idx1); push!(J, idx2); push!(V, 1)
                    push!(I, idx2); push!(J, idx1); push!(V, 1)
                end
            end
        end
    end

    A = sparse(I, J, V, n, n)

    A.nzval .= 1 # Ensure identical directed edges appear as 1

    return A
end


"""
    get_leaves_adjacency_matrix(twoD_Data::AbstractMatrix, ℓ_0::Int, r::Real)

Builds the unweighted adjacency matrix directly from the underlying data and a compression threshold.

# Arguments
- `twoD_Data::AbstractMatrix`: The 2D matrix representing the grid or map.
- `ℓ_0::Int`: Coarsest limit level.
- `r::Real`: Relative error threshold for the adaptive generation.

# Output
- `Tuple{SparseMatrixCSC, Vector}`: The unweighted adjacency matrix and the list of generated leaves.
"""
function get_leaves_adjacency_matrix(twoD_Data, ℓ_0, r)
    A, comp_rate = get_adaptive_matrices(twoD_Data, N_1, r=r)
    G = build_G(A)
    leaves = adaptive_leaf_cells(G, floor(Int, log2(size(twoD_Data, 1))))
    A01 = leaf_adjacency_matrix_directed(leaves, ℓ_0)
    return A01, leaves
end


"""
    get_weighted_adaptive_matrix(A01::SparseMatrixCSC, sorted_leaves::Vector, TopData::TopographicalData, L_geom::Int; Cost_func, raster)

Applies terrain-based travel costs to the unweighted directed adjacency matrix of the adaptive graph.

# Arguments
- `A01::SparseMatrixCSC{Int,Int}`: Unweighted adjacency matrix.
- `sorted_leaves::Vector{Tuple{Int,Int,Int}}`: Vector of adaptive leaves.
- `TopData::TopographicalData`: Topographical data structure.
- `L_geom::Int`: Finest geometric level.

# Keywords
- `Cost_func::Function`: Function computing relative travel cost from slope (default: `adjusted_Hiking`).
- `raster::AbstractMatrix`: Elevation raster (default: `TopData.CostMap`).

# Output
- `SparseMatrixCSC{Float64, Int}`: A weighted sparse adjacency matrix.
"""
function get_weighted_adaptive_matrix(A01::SparseMatrixCSC{Int,Int},
    sorted_leaves::Vector{Tuple{Int,Int,Int}}, TopData::TopographicalData,
    L_geom::Int;
    Cost_func::Function = adjusted_Hiking,
    # αfun::Function = (leaf_u, leaf_v) -> 0.1,
    raster::AbstractMatrix = TopData.CostMap)
    function leaf_block_bounds(leaf::Tuple{Int,Int,Int}, L_geom::Int)
        ℓ, i, j = leaf
        s = 2^(L_geom - ℓ)
        r1 = (i - 1) * s + 1; r2 = i * s; c1 = (j - 1) * s + 1; c2 = j * s
        return r1, r2, c1, c2
    end

    function leaf_mean_height(leaf::Tuple{Int,Int,Int}, raster::AbstractMatrix, L_geom::Int)
        r1, r2, c1, c2 = leaf_block_bounds(leaf, L_geom)
        block = @view raster[r1:r2, c1:c2]
        return mean(block)
    end

    αfun = (leaf_u, leaf_v) -> begin
        ℓu, _, _ = leaf_u
        ℓv, _, _ = leaf_v
        α0 = 0.1; β = 0.5
        α0 * 2.0^(β * (L_geom - min(ℓu, ℓv)))
    end

    n = length(sorted_leaves)
    centers = [ind2coords(sorted_leaves[k], TopData, L_geom) for k in 1:n]
    heights = [leaf_mean_height(sorted_leaves[k], raster, L_geom) for k in 1:n]

    rows, cols, _ = findnz(A01)
    vals = Vector{Float64}(undef, length(rows))

    for t in eachindex(rows)
        v = rows[t]
        u = cols[t]

        # Note: Ad[v,u] represents the cost to travel from u to v
        dist_uv = distance(centers[u], centers[v])

        if dist_uv <= 0
            vals[t] = Inf
            continue
        end

        αuv = αfun(sorted_leaves[u], sorted_leaves[v])
        S_uv = αuv * (heights[v] - heights[u]) / dist_uv
        vals[t] = dist_uv * Cost_func(S_uv)
    end

    return sparse(rows, cols, vals, n, n)

end


##### Flagged for cleanUp
"""
    example_weighted_adaptive_matrix(tiff_file::String, ℓ_0::Int, r::Real)

Example routine demonstrating the creation of a weighted adaptive matrix from a TIFF file.

# Arguments
- `tiff_file::String`: Path to the topographical image file.
- `ℓ_0::Int`: Coarsest level limit.
- `r::Real`: Threshold parameter.

# Output
- `Tuple`: Returns the quadratic `TopoData`, the weighted adjacency matrix, the leaves array, and the max level.
"""
function example_weighted_adaptive_matrix(tiff_file, ℓ_0, r)
    A1 = tiff_to_matrix(tiff_file) 
    TopData = TopographicalData(A1, size(A1), (-13, 35), (46, 70), (1, 1))
    TopoData_quad, maxL = extract_quadratic_map(TopData, 0)
    
    @info "Right-Top Bounds: $(TopoData_quad.rTop), Left-Bottom Bounds: $(TopoData_quad.lBot)"
    
    Ad_leaves01, leaves = get_leaves_adjacency_matrix(TopoData_quad.CostMap, ℓ_0, r)
    @info "Number of generated leaves: $(length(leaves))"
    
    Ad_weighted = get_weighted_adaptive_matrix(Ad_leaves01, leaves, TopData, maxL)
    @info "Matrix size: $(size(Ad_weighted)), Edges: $(nnz(Ad_weighted)/2)"

    return TopoData_quad, Ad_weighted, leaves, floor(Int, log2(size(TopoData_quad.CostMap, 1)))
end