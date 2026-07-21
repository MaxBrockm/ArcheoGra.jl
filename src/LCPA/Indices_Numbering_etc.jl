"""
    characteristic(x::Number, a::Number, b::Number, rBorder::Number)

Calculates the characteristic (indicator) function for a given interval.

Evaluates the interval as `[a, b)` if `b < rBorder`, and as a closed interval 
`[a, b]` if `b == rBorder`.

# Arguments
- `x::Number`: The input value to evaluate.
- `a::Number`: Lower boundary of the interval.
- `b::Number`: Upper boundary of the interval.
- `rBorder::Number`: The absolute right limit of the grid domain.

# Output
- `Int`: Returns `1` if `x` is inside the interval, and `0` otherwise.
"""
@inline function characteristic(x::Number, a::Number, b::Number, rBorder::Number)
    return b==rBorder ? (a<=x<=b ? 1 : 0) : (a<=x<b ? 1 : 0)
end


"""
    coords2lexi(site_coords::Tuple, TopData::AbstractLCPAData)

Computes the linear index in lexicographical order for a given coordinate on a grid.

# Arguments
- `site_coords::Tuple`: A tuple containing `(latitude, longitude)` of the target site.
- `TopData::AbstractLCPAData`: The topographical data structure containing the map grid and boundaries.

# Output
- `Int`: The 1D lexicographical node index corresponding to the coordinate.
"""
@inline function coords2lexi(site_coords::Tuple, TopData::AbstractLCPAData)
    a, b = TopData.lBot
    c, d = TopData.rTop
    m, n = size(TopData.CostMap)

    (lat, long) = site_coords # lat = N,S ; long = E, W

    # Calculate cell dimensions
    x = abs(c-a)/n
    y = abs(d-b)/m

    # Direct calculation instead of loops
    i = clamp(floor(Int, (lat - a) / x) + 1, 1, n)
    j = clamp(floor(Int, (long - b) / y), 0, m - 1)

    # Handle right/top border (closed interval at boundary)
    if lat >= c; i = n; end
    if long >= d; j = m - 1; end

    return j*n+i # Cell of the site in a m×n raster
end

"""
    coords2leaf_index(site_coords::Tuple, leaves::Vector{Tuple{Int,Int,Int}}, TopData::AbstractLCPAData, L_geom::Int)

Finds the linear graph index (from 1 to N) of the adaptive leaf cell containing the specified coordinate.

# Arguments
- `site_coords::Tuple`: A tuple containing `(latitude, longitude)` of the site.
- `leaves::Vector{Tuple{Int,Int,Int}}`: A vector of leaf cells represented as `(ℓ, i, j)`.
- `TopData::AbstractLCPAData`: The typographical data structure containing map parameters.
- `L_geom::Int`: The finest geometric level of the quadtree decomposition.

# Output
- `Int`: The 1D linear index of the leaf within the `leaves` vector.
"""
function coords2leaf_index(site_coords::Tuple, leaves::Vector{Tuple{Int,Int,Int}}, TopData::AbstractLCPAData, L_geom::Int)
    # 1. Find the index on the finest grid
    lexi = coords2lexi(site_coords, TopData)
    fine_i, fine_j = lexi2ind(lexi, TopData)

    # 2. Locate the element in the quadtree
    for (idx, leaf) in enumerate(leaves)
        ℓ, i, j = leaf
        s = 2^(L_geom - ℓ)
        
        r1 = (i - 1) * s + 1
        r2 = i * s
        c1 = (j - 1) * s + 1
        c2 = j * s

        # Check if the fine pixel lies within the block boundary of this leaf
        if r1 <= fine_i <= r2 && c1 <= fine_j <= c2
            return idx
        end
    end
    error("Coordinate $site_coords could not be mapped to any leaf.")
end


"""
    get_xy(Topdata::AbstractLCPAData)

Calculates the cell dimensions of a given (possibly compressed) raster grid.

# Arguments
- `Topdata::AbstractLCPAData`: Structure containing spatial limits, base size, and current resolution parameters.

# Output
- `Tuple{Float64, Float64}`: A tuple `(x, y)` representing the cell width and height, respectively.
"""
@inline function get_xy(Topdata::AbstractLCPAData)
    return (Topdata.rTop.-Topdata.lBot).*reverse((Topdata.kl ./Topdata.size))
end


"""
    ind2coords(ij::Tuple{Int,Int}, Topdata::AbstractLCPAData)

Calculates the geographic center coordinates for a specific 2D matrix entry in a compressed grid.

# Arguments
- `ij::Tuple{Int,Int}`: Matrix indices `(row, column)`.
- `Topdata::AbstractLCPAData`: Structure containing grid resolution and territory boundaries.

# Output
- `Tuple{Float64, Float64}`: A tuple `(cx, cy)` representing the longitude and latitude coordinates of the cell center.
"""
@inline function ind2coords(ij::Tuple{Int,Int}, Topdata::AbstractLCPAData)
    i, j = ij
    m, n = size(Topdata.CostMap)
    
    # 1. Bounds checking
    if i < 1 || i > m || j < 1 || j > n
        throw(BoundsError(Topdata.CostMap, (i, j)))
    end
    
    x, y = get_xy(Topdata)
    a = Topdata.size[1] ÷ Topdata.kl[1]
    b = Topdata.size[2] ÷ Topdata.kl[2]
    
    lBot_x, lBot_y = Topdata.lBot
    rTop_x, rTop_y = Topdata.rTop
    
    # Calculate X-coordinate (depends only on j)
    if j <= b
        cx = lBot_x + (j - 0.5) * x
    else
        distx = rTop_x - lBot_x
        cx = rTop_x - 0.5 * (distx - b * x)
    end
    
    # Calculate Y-coordinate (depends only on i)
    if i <= a
        cy = rTop_y - (i - 0.5) * y
    else
        disty = rTop_y - lBot_y
        cy = lBot_y + 0.5 * (disty - a * y)
    end
    
    return (cx, cy)
end


"""
    ind2coords(leaf::Tuple{Int,Int,Int}, TopData::AbstractLCPAData, L_geom::Int)

Calculates the geographic center coordinates of an adaptive quadtree leaf cell.

Assumes that `TopData.CostMap` is square and sized in powers of two.

# Arguments
- `leaf::Tuple{Int,Int,Int}`: A tuple `(ℓ, i, j)` specifying the level, row, and column of the leaf.
- `TopData::AbstractLCPAData`: Topographical data structure containing boundary definitions.
- `L_geom::Int`: The finest geometric level of the quadtree decomposition.

# Output
- `Tuple{Float64, Float64}`: A tuple `(cx, cy)` representing the center coordinates of the adaptive leaf cell.
"""
@inline function ind2coords(leaf::Tuple{Int,Int,Int}, TopData::AbstractLCPAData, L_geom::Int)
    ℓ, i, j = leaf

    0 <= ℓ <= L_geom || error("Leaf level ℓ=$ℓ is not compatible with L_geom=$L_geom")
    1 <= i <= 2^ℓ || error("Row index i must be in 1:$(2^ℓ)")
    1 <= j <= 2^ℓ || error("Column index j must be in 1:$(2^ℓ)")

    N = 2^L_geom
    Δx = (TopData.rTop[1] - TopData.lBot[1]) / N
    Δy = (TopData.rTop[2] - TopData.lBot[2]) / N

    s = 2^(L_geom - ℓ)

    x_left   = TopData.lBot[1] + (j - 1) * s * Δx
    x_right  = TopData.lBot[1] + j * s * Δx
    y_top    = TopData.rTop[2] - (i - 1) * s * Δy
    y_bottom = TopData.rTop[2] - i * s * Δy

    return ((x_left + x_right) / 2, (y_bottom + y_top) / 2) # Coordinate of middlepoint
end


"""
    lexiMat2Vec(Mat::SparseMatrixCSC{Bool, Int64})

Converts a 2D sparse matrix into a 1D boolean vector using a custom lexicographical conversion.

# Arguments
- `Mat::SparseMatrixCSC{Bool, Int64}`: The sparse input matrix.

# Output
- `Vector{Bool}`: A flattened 1D array mapped in lexicographical order.
"""
function lexiMat2Vec(Mat::SparseMatrixCSC{Bool, Int64})
    m, n = size(Mat)
    I_idx, J_idx, _ = findnz(Mat)
    Boolvec = zeros(Bool, m*n)
    for k in eachindex(I_idx)
        Boolvec[(m-I_idx[k])*n+J_idx[k]] = 1
    end
    return Boolvec
end

"""
    lexi2ind(lexi::Int, TopData::TopographicalData)

Converts a 1D lexicographical grid index back to 2D matrix indices `(row, column)`.

# Arguments
- `lexi::Int`: The linear lexicographical cell index.
- `TopData::TopographicalData`: The underlying topographic data structure defining the grid dimensions.

# Output
- `Tuple{Int, Int}`: A tuple `(i, j)` containing the matrix row and column indices.
"""
@inline function lexi2ind(lexi::Int, TopData::TopographicalData)
    m, n = size(TopData.CostMap) 
    i = m-ceil(lexi/n)+1
    j = mod(lexi,n)==0 ? n : mod(lexi, n)
    return (Int(i), Int(j))
end