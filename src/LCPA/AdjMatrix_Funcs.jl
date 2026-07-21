"""
    weighted_adjacency_matrix(cost_raster::AbstractMatrix{<:Number})

Computes a symmetric adjacency matrix based on a cost raster grid.

The connectivity is determined using 8-neighbourhood relations.
For horizontal and vertical neighbors, edge weights are calculated as the average cost 
of the two adjacent cells. For diagonal neighbors, edge weights are calculated as 
the average cost divided by `√2`.

# Arguments
- `cost_raster::AbstractMatrix{<:Number}`: A 2D matrix containing terrain costs.

# Output
- `SparseMatrixCSC{Float64, Int}`: A symmetric sparse adjacency matrix of size `(m * n) × (m * n)`, 
  where `m` and `n` are the dimensions of the input cost raster.
"""
function weighted_adjacency_matrix(cost_raster::AbstractMatrix{<:Number})
    (m, n) = size(cost_raster)

    # Vertical connections
    valvec = Vector{Float64}(undef, n*(m-1))
    rowvec = collect((n+1):(m*n))
    colvec = collect(1:(m-1)*n)
    for j=1:m-1
        valvec[(j-1)*n+1:j*n] = 0.5*(cost_raster[m-j+1,:] + cost_raster[m-j,:])
    end
    Ad1 = sparse(rowvec, colvec, valvec, m*n, m*n)

    # Horizontal connections
    valvec = Vector{Float64}(undef, m*(n-1))
    colvec = zeros(m*(n-1))
    rowvec = zeros(m*(n-1))
    for j=1:n-1
        valvec[(j-1)*m+1:j*m] = 0.5*(cost_raster[:, j+1] + cost_raster[:, j])
        colvec[(j-1)*m+1:j*m] = collect(((m*n)+j+1-n:-n:(j+1)))
        rowvec[(j-1)*m+1:j*m] = collect((m*n)+j-n:-n:j)
    end
    Ad2 = sparse(rowvec, colvec, valvec, m*n, m*n)

    # Diagonal connections (top right to left bottom)
    valvec = Vector{Float64}(undef, (m-1)*(n-1))
    rowvec = zeros((m-1)*(n-1))
    colvec = zeros((m-1)*(n-1))
    for j=1:m-1
        valvec[(j-1)*(n-1)+1:j*(n-1)] = (cost_raster[m-j+1, 1:n-1] + cost_raster[m-j, 2:n])/sqrt(2)
        rowvec[(j-1)*(n-1)+1:j*(n-1)] = collect((j*n+2):(j+1)*n)
        colvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j-1)*n+1):(j*n-1))
    end
    Ad3 = sparse(rowvec, colvec, valvec, m*n, m*n)

    # Diagonal connections (top left to right bottom)
    valvec = Vector{Float64}(undef, (m-1)*(n-1))
    rowvec = zeros((m-1)*(n-1))
    colvec = zeros((m-1)*(n-1))
    for j=1:m-1
        valvec[(j-1)*(n-1)+1:j*(n-1)] = (cost_raster[m-j, 1:n-1]+cost_raster[m-j+1, 2:n])/sqrt(2)
        colvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j-1)*n+2):j*n)
        rowvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j*n+1):((j+1)*n-1)))
    end
    Ad4 = sparse(rowvec, colvec, valvec, m*n, m*n)
    Ad = Ad1 + Ad2 + Ad3 + Ad4
    Ad = Ad'+ Ad

    @info "Adjacency matrix is done"
    return Ad
end

"""
    weighted_adjacency_matrix_haversine(TopData::AbstractLCPAData, Cost_func::Function)

Constructs an asymmetric weighted adjacency matrix using Haversine distances and elevation differences.

Slope values between neighbouring cells are calculated using the geographic distance 
and elevation differences. These slope values are then passed through the provided 
cost function to determine the edge weights.

# Arguments
- `TopData::AbstractLCPAData`: Topographical data structure containing the cost map and boundaries.
- `Cost_func::Function`: A function mapping slope (rise over run) to travel costs.

# Output
- `SparseMatrixCSC{Float64, Int}`: An asymmetric sparse adjacency matrix of size `(m * n) × (m * n)`.
"""
function weighted_adjacency_matrix_haversine(TopData::AbstractLCPAData, Cost_func::Function)
    elevation_raster = TopData.CostMap
    typeof(Float64.(elevation_raster)) == Vector{Float64} ? (m,n) = (length(elevation_raster), 1) : (m, n) = size(elevation_raster)

    # Vertical connections
    valvec = Vector{Float64}(undef, n*(m-1))
    colvec = collect((n+1):(m*n))
    rowvec = collect(1:(m-1)*n)
    for j=1:m-1
        dist = distance(ind2coords((m-j, 1), TopData), ind2coords((m-j+1, 1), TopData))
        valvec[(j-1)*n+1:j*n] = 0.1*(elevation_raster[m-j+1,:] - elevation_raster[m-j,:])/dist
    end
    
    AdI1 = sparse(rowvec, colvec, Cost_func.(valvec), m*n, m*n)
    AdII1 = sparse(rowvec, colvec, Cost_func.(-valvec), m*n, m*n)

    # Horizontal connections
    valvec = Vector{Float64}(undef, m*(n-1))
    colvec = zeros(m*(n-1))
    rowvec = zeros(m*(n-1))
    
    disthori = [distance(ind2coords((j,1), TopData), ind2coords((j,2), TopData)) for j=1:m]
    for j=1:n-1
        valvec[(j-1)*m+1:j*m] = 0.1*(elevation_raster[:, j+1] - elevation_raster[:, j])./disthori
        colvec[(j-1)*m+1:j*m] = collect(((m*n)+j+1-n:-n:(j+1)))
        rowvec[(j-1)*m+1:j*m] = collect((m*n)+j-n:-n:j)
    end
    AdI2 = sparse(rowvec, colvec, Cost_func.(-valvec), m*n, m*n)
    AdII2 = sparse(rowvec, colvec, Cost_func.(valvec), m*n, m*n)

    # Diagonal connections (top right to left bottom)
    valvec = Vector{Float64}(undef, (m-1)*(n-1))
    colvec = zeros((m-1)*(n-1))
    rowvec = zeros((m-1)*(n-1))
    for j=1:m-1
        dist = distance(ind2coords((m-j, 1), TopData), ind2coords((m-j+1, 2), TopData))
        
        valvec[(j-1)*(n-1)+1:j*(n-1)] = 0.1*(elevation_raster[m-j+1, 1:n-1] - elevation_raster[m-j, 2:n])/dist
        colvec[(j-1)*(n-1)+1:j*(n-1)] = collect((j*n+2):(j+1)*n)
        rowvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j-1)*n+1):(j*n-1))
    end
    AdI3 = sparse(rowvec, colvec, Cost_func.(valvec), m*n, m*n)
    AdII3 = sparse(rowvec, colvec, Cost_func.(-valvec), m*n, m*n)
    
    # Diagonal connections (top left to right bottom)
    valvec = Vector{Float64}(undef, (m-1)*(n-1))
    colvec = zeros((m-1)*(n-1))
    rowvec = zeros((m-1)*(n-1))
    for j=1:m-1
        dist = distance(ind2coords((m-j, 2), TopData), ind2coords((m-j+1, 1), TopData))
        
        valvec[(j-1)*(n-1)+1:j*(n-1)] = 0.1*(elevation_raster[m-j, 1:n-1] - elevation_raster[m-j+1, 2:n])/dist
        rowvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j-1)*n+2):j*n)
        colvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j*n+1):((j+1)*n-1)))
    end
    AdI4 = sparse(rowvec, colvec, Cost_func.(-valvec), m*n, m*n)
    AdII4 = sparse(rowvec, colvec, Cost_func.(valvec), m*n, m*n)

    AdI = AdI1 + AdI2 + AdI3 + AdI4
    AdII = AdII1 + AdII2 + AdII3 + AdII4
    Ad =  AdI + AdII'
    @info "Haversine Cost Function Adjacency matrix is done"
    return Ad
end

"""
    WeightConfig

A configuration container defining directional weighting and distance functions.

This structure specifies how vertical, horizontal, and diagonal relationships 
between grid cells are combined and normalized for graph construction.

# Fields
- `vert_func::Function`: Combination rule for vertical neighbors.
- `hori_func::Function`: Combination rule for horizontal neighbors.
- `diag1_func::Function`: Combination rule for the first diagonal neighbors (top-right to bottom-left).
- `diag2_func::Function`: Combination rule for the second diagonal neighbors (top-left to bottom-right).
- `vert_dist::Function`: Distance or normalization function for the vertical direction.
- `hori_dist::Function`: Distance or normalization function for the horizontal direction.
- `diag1_dist::Function`: Distance or normalization function for the first diagonal direction.
- `diag2_dist::Function`: Distance or normalization function for the second diagonal direction.
- `symmetric::Bool`: Set to `true` if weights are symmetric (`w(u, v) == w(v, u)`), `false` otherwise.

# Example
```julia
# Simple additive configuration
cfg = WeightConfig(
    (a,b)->a+b, (a,b)->a+b, (a,b)->a+b, (a,b)->a+b,
    (a,b)->1.0, (a,b)->1.0, (a,b)->√2, (a,b)->√2,
    true
)

cfg.vert_func(1, 2)     # → 3
cfg.diag1_dist(1, 1.5)  # → 1.414...
```
"""
struct WeightConfig
    vert_func::Function
    hori_func::Function
    diag1_func::Function
    diag2_func::Function
    vert_dist::Function
    hori_dist::Function
    diag1_dist::Function
    diag2_dist::Function
    symmetric::Bool
end

"""
    weighted_adjacency_matrix_combi(LCPAData::AbstractLCPAData, WConfig::WeightConfig; Cost_func::Function = x -> x)

Constructs an adjacency matrix combining a customizable weight configuration and a cost function.

Computes transition costs between adjacent cells by dividing the combined cell values 
by the directional distance, as defined in `WConfig`. The result is then evaluated 
using the provided `Cost_func`.

# Arguments
- `LCPAData::AbstractLCPAData`: Data structure containing the grid (e.g., `CostMap`).
- `WConfig::WeightConfig`: Configuration defining interpolation and distance functions.

# Keywords
- `Cost_func::Function`: A function applied to the computed slopes (default: identity `x -> x`).

# Output
- `SparseMatrixCSC{Float64, Int}`: Sparse adjacency matrix of size `(m * n) × (m * n)`.
"""
function weighted_adjacency_matrix_combi(LCPAData::AbstractLCPAData, WConfig::WeightConfig; Cost_func::Function = x -> x)
    if typeof(LCPAData) == TopographicalData
        cost_raster = LCPAData.CostMap ### CostMap must be in km
    end
    typeof(Float64.(cost_raster)) == Vector{Float64} ? (m,n) = (length(cost_raster), 1) : (m, n) = size(cost_raster)

    # Vertical connections
    valvec = zeros(n*(m-1)) 
    colvec = collect((n+1):(m*n))
    rowvec = collect(1:(m-1)*n)
    for j=1:m-1
        dist = WConfig.vert_dist(ind2coords((m-j, 1), LCPAData), ind2coords((m-j+1, 1), LCPAData))
        valvec[(j-1)*n+1:j*n] = WConfig.vert_func(cost_raster[m-j+1,:], cost_raster[m-j,:])/dist
    end
    AdI1 = sparse(rowvec, colvec, Cost_func.(valvec), m*n, m*n)
    AdII1 = sparse(rowvec, colvec, Cost_func.(-valvec), m*n, m*n)

    # Horizontal connections
    valvec = zeros(m*(n-1)) 
    colvec = zeros(m*(n-1))
    rowvec = zeros(m*(n-1))
    
    disthori = [WConfig.hori_dist(ind2coords((j,1), LCPAData), ind2coords((j,2), LCPAData)) for j=1:m]
    for j=1:n-1
        valvec[(j-1)*m+1:j*m] = WConfig.hori_func(cost_raster[:, j+1], cost_raster[:, j])./disthori
        colvec[(j-1)*m+1:j*m] = collect(((m*n)+j+1-n:-n:(j+1)))
        rowvec[(j-1)*m+1:j*m] = collect((m*n)+j-n:-n:j)
    end
    AdI2 = sparse(rowvec, colvec, Cost_func.(-valvec), m*n, m*n)
    AdII2 = sparse(rowvec, colvec, Cost_func.(valvec), m*n, m*n)

    # Diagonal connections (top right to left bottom)
    valvec = zeros((m-1)*(n-1)) 
    colvec = zeros((m-1)*(n-1))
    rowvec = zeros((m-1)*(n-1))
    for j=1:m-1
        dist = WConfig.diag1_dist(ind2coords((m-j, 1), LCPAData), ind2coords((m-j+1, 2), LCPAData))

        valvec[(j-1)*(n-1)+1:j*(n-1)] = WConfig.diag1_func(cost_raster[m-j+1, 1:n-1], cost_raster[m-j, 2:n])/dist
        colvec[(j-1)*(n-1)+1:j*(n-1)] = collect((j*n+2):(j+1)*n)
        rowvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j-1)*n+1):(j*n-1))
    end
    AdI3 = sparse(rowvec, colvec, Cost_func.(valvec), m*n, m*n)
    AdII3 = sparse(rowvec, colvec, Cost_func.(-valvec), m*n, m*n)
    
    # Diagonal connections (top left to right bottom)
    valvec = zeros((m-1)*(n-1))
    colvec = zeros((m-1)*(n-1))
    rowvec = zeros((m-1)*(n-1))
    for j=1:m-1
        dist = WConfig.diag2_dist(ind2coords((m-j, 2), LCPAData), ind2coords((m-j+1, 1), LCPAData))
        
        valvec[(j-1)*(n-1)+1:j*(n-1)] = WConfig.diag2_func(cost_raster[m-j, 1:n-1], cost_raster[m-j+1, 2:n])/dist
        rowvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j-1)*n+2):j*n)
        colvec[(j-1)*(n-1)+1:j*(n-1)] = collect(((j*n+1):((j+1)*n-1)))
    end
    AdI4 = sparse(rowvec, colvec, Cost_func.(-valvec), m*n, m*n)
    AdII4 = sparse(rowvec, colvec, Cost_func.(valvec), m*n, m*n)

    if WConfig.symmetric == false
        AdI = AdI1 + AdI2 + AdI3 + AdI4
        AdII = AdII1 + AdII2 + AdII3 + AdII4
        @debug "Cost Function Adjacency matrix is done"
        return AdI + AdII'
    elseif WConfig.symmetric == true
        @debug "Cost Function Adjacency matrix is done"
        Ad = AdI1 + AdII2 + AdI3 + AdII4
        return Ad + Ad'
    end
end

"""
    Kondo_Average

A pre-defined symmetric `WeightConfig` based on the Kondo averaging scheme.

This configuration defines symmetric weighting where simple addition is used to combine 
neighboring values, and constant factors (`2` and `√2`) are used for spatial normalization.
The final argument is set to `true` to declare that weights are symmetric.

# Example
```julia
Kondo_Average.combine_x(1.0, 2.0)  # → 3.0
```
"""
const Kondo_Average = WeightConfig(
    (a,b) -> a+b,
    (a,b) -> a+b,
    (a,b) -> a+b,
    (a,b) -> a+b,
    (a,b) -> 2,
    (a,b) -> 2,
    (a,b) -> sqrt(2),
    (a,b) -> sqrt(2),
    true
)

"""
    Haversine_Config

A pre-defined asymmetric `WeightConfig` using the Haversine distance metric.

This configuration measures directed differences (`a - b`) and uses the real geographic 
distance between cell coordinates to compute slope-like relationships.
The final argument is set to `false` to indicate non-symmetric weighting.

# Example
```julia
Haversine_Config.dist_x(lat1, lat2)
```
"""
const Haversine_Config = WeightConfig(
    (a,b) -> a-b,
    (a,b) -> a-b,
    (a,b) -> a-b,
    (a,b) -> a-b,
    (a,b) -> distance(a,b),
    (a,b) -> distance(a,b),
    (a,b) -> distance(a,b),
    (a,b) -> distance(a,b),
    false
)