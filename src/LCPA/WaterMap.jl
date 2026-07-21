"""
    WaterLandBool(TopData::TopographicalData)

Creates a sparse boolean matrix indicating whether each grid cell is water (0) or land (1).

The function iterates over the entire cost map and evaluates the geographic 
coordinates against a global land-sea mask.

# Arguments
- `TopData::TopographicalData`: The topographical data structure containing the grid size and boundaries.

# Output
- `SparseMatrixCSC{Bool, Int}`: A sparse matrix of the same dimensions as `TopData.CostMap`, 
  containing `true` (1) for land and `false` (0) for water.
"""
function WaterLandBool(TopData::TopographicalData)
    mask = GeoDatasets.LandSeaMask(;resolution='f', grid=1.25)
    WaterLandBool = zeros(Bool, size(TopData.CostMap))
    
    @showprogress for i in 1:size(TopData.CostMap)[1]
        for j in 1:size(TopData.CostMap)[2]
            Coords = ind2coords((i, j), TopData)
            WaterLandBool[i,j] = GeoDatasets.is_land(mask, Coords[2], Coords[1])
        end
    end
    @info "Creating WaterLandBool is done"
    return sparse(WaterLandBool)
end


"""
    adaptive_WaterLandBool(leaves::Vector{Tuple{Int,Int,Int}}, TopData::TopographicalData, L_geom::Int, existence_list_leaves::AbstractVector{Int})

Determines which adaptive leaf cells contain land elements to filter out pure water cells.

The function tests five points per leaf cell: the center and the four corners. 
If any of these points lie on land, the entire cell is kept. Additionally, 
cells containing target sites (specified in `existence_list_leaves`) are unconditionally preserved.

# Arguments
- `leaves::Vector{Tuple{Int,Int,Int}}`: A vector of adaptive quadtree leaf cells `(ℓ, i, j)`.
- `TopData::TopographicalData`: The topographical data structure.
- `L_geom::Int`: The finest geometric level of the quadtree decomposition.
- `existence_list_leaves::AbstractVector{Int}`: Indices of leaves containing the start/target sites (these are never deleted).

# Output
- `Vector{Bool}`: A boolean array of the same length as `leaves`, containing `true` if the leaf should be kept.
"""
function adaptive_WaterLandBool(leaves::Vector{Tuple{Int,Int,Int}}, TopData::TopographicalData, L_geom::Int, existence_list_leaves::AbstractVector{Int})
    mask = GeoDatasets.LandSeaMask(;resolution='f', grid=1.25)
    N_leaves = length(leaves)
    keepNodes = zeros(Bool, N_leaves)
    
    # Grid dimensions for corner point calculations
    N_grid = 2^L_geom
    Δx = (TopData.rTop[1] - TopData.lBot[1]) / N_grid
    Δy = (TopData.rTop[2] - TopData.lBot[2]) / N_grid
    
    @showprogress for k in 1:N_leaves
        if k in existence_list_leaves
            keepNodes[k] = true
            continue
        end
        
        ℓ, i, j = leaves[k]
        s = 2^(L_geom - ℓ)
        
       # Calculate coordinates of the block boundaries
        x_left   = TopData.lBot[1] + (j - 1) * s * Δx
        x_right  = TopData.lBot[1] + j * s * Δx
        y_top    = TopData.rTop[2] - (i - 1) * s * Δy
        y_bottom = TopData.rTop[2] - i * s * Δy
        
       # 5 test points: Center + 4 corners
        pts = [
            ((x_left + x_right) / 2, (y_top + y_bottom) / 2), # Center
            (x_left, y_bottom),                               # Bottom Left
            (x_right, y_bottom),                              # Bottom Right
            (x_left, y_top),                                  # Top Left
            (x_right, y_top)                                  # Top Right
        ]
        
        # If ANY point is on land, keep the cell (true)
        keepNodes[k] = any(pt -> GeoDatasets.is_land(mask, pt[2], pt[1]), pts)
    end
    
    @info "Creating defensive adaptive WaterLandBool is done"
    return keepNodes
end
