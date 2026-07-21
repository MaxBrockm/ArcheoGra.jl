"""
    coeffi_mat_Mℓ_0(a_mask::AbstractVector, ℓ::Int, Delta::Function)

Constructs the scaling coefficient matrix `Mℓ_0` for a given level `ℓ` using the scaling mask `a_mask`.  

This matrix represents the refinement relation between scaling functions across 
successive levels (low-pass filtering).

# Arguments
- `a_mask::AbstractVector`: Scaling mask (low-pass filter) coefficients.
- `ℓ::Int`: Current scale (level) in the multiresolution decomposition.
- `Delta::Function`: Function returning the number of scaling coefficients at a given level.

# Output
- `SparseMatrixCSC{Float64, Int}`: Sparse matrix defining the scaling (approximation) transformation between levels `ℓ` and `ℓ + 1`.
"""
function coeffi_mat_Mℓ_0(a_mask::AbstractVector, ℓ, Delta::Function)
    n = length(a_mask)-1
    Mℓ_0 = spzeros(Delta(ℓ+1), Delta(ℓ))
    
    for index in 0:Delta(ℓ)-1
        Mℓ_0[(1+2*index):(2*index+n+1), index+1] = a_mask
    end

    return Mℓ_0
end

"""
    coeffi_mat_Mℓ_1(b_mask::AbstractVector, ℓ::Int, Delta::Function, Nabla::Function)

Constructs the wavelet coefficient matrix `Mℓ_1` for a given level `ℓ` using the wavelet mask `b_mask`.

This matrix encodes the relationship between scaling and wavelet coefficients across levels 
(high-pass filtering).

# Arguments
- `b_mask::AbstractVector`: Wavelet mask (high-pass filter) coefficients.
- `ℓ::Int`: Current scale (level) in the multiresolution decomposition.
- `Delta::Function`: Function returning the number of scaling coefficients at a given level.
- `Nabla::Function`: Function returning the number of wavelet coefficients at a given level.

# Output
- `SparseMatrixCSC{Float64, Int}`: Sparse matrix defining the wavelet (detail) transformation between levels `ℓ` and `ℓ + 1`.
"""
function coeffi_mat_Mℓ_1(b_mask::AbstractVector, ℓ, Delta::Function, Nabla::Function)
    n = length(b_mask)-1
    Mℓ_1 = spzeros(Delta(ℓ+1), Nabla(ℓ))

    for index in 0:Delta(ℓ)-1
        Mℓ_1[(1+2*index):(2*index+n+1), index+1] = b_mask
    end

    return Mℓ_1
end

"""
    Mℓ_and_Gℓ(ℓ::Int, N_m::N_m_and_Wavelet)

Constructs the transformation matrices `Mℓ_0`, `Mℓ_1`, `Gℓ_0`, and `Gℓ_1` used in multiscale wavelet analysis.

These matrices depend on the chosen wavelet type `N_m` (e.g., `N_1`, `N_2`, `N_2_boundary`) 
and the scale level `ℓ`. They define the transformations between different resolution levels.

# Arguments
- `ℓ::Int`: The current scale (level) of the multiresolution analysis.
- `N_m::N_m_and_Wavelet`: A structure containing the filter masks and the level sizing functions.

# Output
- `Tuple`: Returns `(Mℓ_0, Mℓ_1, Gℓ_0, Gℓ_1)`:
  - `Mℓ_0`: Analysis scaling matrix.
  - `Mℓ_1`: Analysis wavelet matrix.
  - `Gℓ_0`: Reconstruction scaling matrix.
  - `Gℓ_1`: Reconstruction wavelet matrix.

# Example
```julia
Mℓ_0, Mℓ_1, Gℓ_0, Gℓ_1 = Mℓ_and_Gℓ(3, N_2)
```
"""
function Mℓ_and_Gℓ(ℓ, N_m::N_m_and_Wavelet) # Δ_ℓ = 2^ℓ for N_1 bzw 2^ℓ-1 for N_2 bzw 2^ℓ+1 for N_2^⋆
    Mℓ_0 = spzeros(N_m.Delta_l(ℓ+1), N_m.Delta_l(ℓ))
    Mℓ_1 = spzeros(N_m.Delta_l(ℓ+1), N_m.Nabla_l(ℓ))
    if N_m == N_1 || N_m == N_2
        Mℓ_0 = (1/sqrt(2))*coeffi_mat_Mℓ_0(N_m.Mask, ℓ, N_m.Delta_l)
        Mℓ_1 = (1/sqrt(2))*coeffi_mat_Mℓ_1(N_m.Wavelet_Mask, ℓ, N_m.Delta_l, N_m.Nabla_l)
    elseif N_m == N_2_boundary
        n = length(N_m.Mask)-1
        Mℓ_0[1:2,1] = [1, 1/2]
        for index in 1:N_m.Delta_l(ℓ)-2
            Mℓ_0[(2*index):(2*index+n), index+1] = N_m.Mask
        end
        Mℓ_0[end-1:end, end] = [1/2, 1]
        Mℓ_0 /= sqrt(2)

        n = length(N_m.Wavelet_Mask)-1
        Mℓ_1[1:6,1] = [3/4, -7/16, 3/8, -11/16, 1/4, 1/8]
        for index in 1:N_m.Delta_l(ℓ)-3
            Mℓ_1[(2*index):(2*index+n), index+1] = N_m.Wavelet_Mask
        end
        Mℓ_1[end-5:end, N_m.Delta_l(ℓ)-1] = [1/8, 1/4, -11/16, 3/8, -7/16, 3/4]
        Mℓ_1 /= sqrt(2)
    end
    Gℓ_0 = spzeros(N_m.Delta_l(ℓ+1), N_m.Delta_l(ℓ))
    Gℓ_1 = spzeros(N_m.Delta_l(ℓ+1), N_m.Nabla_l(ℓ))
    if N_m == N_1
        Gℓ_0 = transpose(Mℓ_0)
        Gℓ_1 = transpose(Mℓ_1)
    elseif N_m == N_2
        Mℓ_1[end,end] = (1/sqrt(2))
        Gℓ_0 = sqrt(2)*coeffi_mat_Mℓ_0(N_m.inv_Mask, ℓ, N_m.Delta_l)'
        Gℓ_1_temp = spzeros(N_m.Delta_l(ℓ+1), N_m.Nabla_l(ℓ))
        Gℓ_1_temp[1:2, 1] = [1, -1/2]; Gℓ_1_temp[end-1:end, end] = [-1/2, 1]
        n = length(N_m.inv_Wavelet_Mask)-1
        for index in 1:N_m.Delta_l(ℓ)-1
            Gℓ_1_temp[(2*index):(2*index+n), index+1] = N_m.inv_Wavelet_Mask
        end
        Gℓ_1 = sqrt(2)*Gℓ_1_temp'
    elseif N_m == N_2_boundary
        n = length(N_m.Wavelet_Mask)-1
        Gℓ_0[1:3,1] = [5/4, 3/2, -3/4]
        Gℓ_0[1:5,2] = [-1/8, 1/4, 13/8, 1/2, -1/4]
        for index in 2:N_m.Delta_l(ℓ)-2
            Gℓ_0[(2*index-1):(2*index+n-1), index+1] = N_m.inv_Mask
        end
        Gℓ_0[end-4:end, end-1] = [-1/4, 1/2, 13/8, 1/4, -1/8]
        Gℓ_0[end-2:end, end] = [-3/4, 3/2, 5/4]

        n = length(N_m.inv_Wavelet_Mask)-1
        Gℓ_1[1:3,1] = [1, -2, 1]
        Gℓ_1[1:5,2] = [1/2, -1, 0, 1, -1/2]
        for index in 2:2^ℓ-3
            Gℓ_1[(2*index+1):(2*index+n+1), index+1] = N_m.inv_Wavelet_Mask
        end
        Gℓ_1[end-4:end, end-1] = [-1/2, 1, 0, -1, 1/2]
        Gℓ_1[end-2:end, end] = [1, -2, 1]
        Gℓ_0 = (1/sqrt(2))*Gℓ_0'
        Gℓ_1 = (1/sqrt(2))*Gℓ_1'
    end
    return Mℓ_0, Mℓ_1, Gℓ_0, Gℓ_1
end

"""
    create_Tℓ_inv(Gℓ_0::SparseMatrixCSC, Gℓ_1::SparseMatrixCSC, ℓ::Int, L::Int, Delta_l::Function)

Constructs the inverse transformation matrix `Tℓ_inv` at level `ℓ` for an overall level `L`.

This matrix combines the reconstruction operators `Gℓ_0` and `Gℓ_1` and extends them to match  
the full resolution dimension defined by `Delta_l(L)`.

# Arguments
- `Gℓ_0::SparseMatrixCSC`: Reconstruction matrix for the scaling space at level `ℓ`.
- `Gℓ_1::SparseMatrixCSC`: Reconstruction matrix for the wavelet space at level `ℓ`.
- `ℓ::Int`: Current level of decomposition.
- `L::Int`: Final (maximum) level of decomposition.
- `Delta_l::Function`: Function returning the number of coefficients at a given level.

# Output
- `SparseMatrixCSC`: Sparse inverse transformation matrix for level `ℓ`.
"""
function create_Tℓ_inv(Gℓ_0, Gℓ_1, ℓ::Int, L::Int, Delta_l::Function)
    Gℓ = vcat(Gℓ_0, Gℓ_1)
    if L == ℓ+1 
        return Gℓ
    else
        Tℓ_inv = [Gℓ spzeros(size(Gℓ,1), Delta_l(L)-Delta_l(ℓ+1)); spzeros(Delta_l(L)-Delta_l(ℓ+1), size(Gℓ,2)) sparse(I, Delta_l(L)-Delta_l(ℓ+1), Delta_l(L)-Delta_l(ℓ+1))]#blockdiag(vcat(Gℓ_0,Gℓ_1), sparse(I, 2^L-2^(ℓ+1), 2^L-2^(ℓ+1)))
        return Tℓ_inv
    end
end

"""
    create_Tℓ(Mℓ_0::SparseMatrixCSC, Mℓ_1::SparseMatrixCSC, ℓ::Int, L::Int, Delta_l::Function)

Constructs the forward transformation matrix `Tℓ` at level `ℓ` for an overall level `L`.

This matrix combines the analysis operators `Mℓ_0` and `Mℓ_1` and extends them  
to the full resolution defined by `Delta_l(L)`.

# Arguments
- `Mℓ_0::SparseMatrixCSC`: Analysis matrix for scaling coefficients.
- `Mℓ_1::SparseMatrixCSC`: Analysis matrix for wavelet coefficients.
- `ℓ::Int`: Current level of decomposition.
- `L::Int`: Final (maximum) level of decomposition.
- `Delta_l::Function`: Function returning the number of coefficients at a given level.

# Output
- `SparseMatrixCSC`: Sparse forward transformation matrix for level `ℓ`.
"""
function create_Tℓ(Mℓ_0, Mℓ_1, ℓ::Int, L::Int, Delta_l::Function)
    Mℓ = hcat(Mℓ_0,Mℓ_1)
    if L == ℓ+1 
        return Mℓ
    else
        Tℓ = [Mℓ spzeros(size(Mℓ,1), Delta_l(L)-Delta_l(ℓ+1)); spzeros(Delta_l(L)-Delta_l(ℓ+1), size(Mℓ,2)) sparse(I, Delta_l(L)-Delta_l(ℓ+1), Delta_l(L)-Delta_l(ℓ+1))]#blockdiag(hcat(Mℓ_0,Mℓ_1), sparse(I, 2^L-2^(ℓ+1), 2^L-2^(ℓ+1)))
        return Tℓ
    end
end

"""
    compression_1D(data_input::AbstractVector, ℓ_0::Int; N_m::N_m_and_Wavelet=N_1, N::Int=-1, r::Real=-1)

Performs one-dimensional wavelet compression on an input signal.

Decomposes, thresholds, and reconstructs the data according to the chosen  
wavelet model `N_m` and compression strategy. Only one thresholding method (`N` or `r`) can be used at a time.

# Arguments
- `data_input::AbstractVector`: Input signal to be compressed.
- `ℓ_0::Int`: Lowest decomposition level (coarsest scale).

# Keywords
- `N_m::N_m_and_Wavelet`: Wavelet model defining the masks and level functions (default: `N_1`).
- `N::Int`: Hard threshold - retain only the `N` largest coefficients (default: `-1` / disabled).
- `r::Real`: Relative error threshold, where `r ∈ (0, 1)` (default: `-1` / disabled).

# Output
- `Vector{Float64}`: The reconstructed (compressed) signal.
"""
function compression_1D(data_input, ℓ_0; N=-1, r=-1, N_m::N_m_and_Wavelet=N_1)
    Delta(x) = N_m.Delta_l(x)
    c = copy(data_input) # c_L
    L = 0
    
    if N_m == N_2
        try 
            L = Int(log2(size(data_input,1)+1))
        catch 
            error("Input size must be 2^L-1 for N_2, L Integer")
        end
    elseif N_m == N_2_boundary
        try 
            L = Int(log2(size(data_input,1)-1))
        catch 
            error("Input size must be 2^L+1 for N_2 with boundary adaptations, L Integer")
        end
    else
        try 
            L = Int(log2(size(data_input,1)))
        catch 
            error("Input size must be 2^L for N_1, L Integer")
        end
    end

    p1 = plot(c, label="input", lw=2)

    # Inverse Transform
    dL = zeros(Delta(L))
    for ℓ in (L-1):-1:ℓ_0
        _,_, Gℓ_0, Gℓ_1 = Mℓ_and_Gℓ(ℓ, N_m)
        dL[Delta(ℓ)+1:Delta(ℓ+1)] = Gℓ_1*c
        c = Gℓ_0*c
    end

    # Thresholding
    if N != -1 && r != -1
        error("Choose a way of thresholding")
    elseif N != -1
        dL = largestNTerm(dL, N)
    elseif r != -1 
        dL = relErrorTerm(dL, r)
    end
    printstyled("Compression rate=$(1 - (count(!iszero, dL))/(Delta(L)))\n", color=:blue)

    dL[1:Delta(ℓ_0)] = c

    # Forward Transform (Reconstruction)
    for ℓ in ℓ_0:1:L-1
        Mℓ_0, Mℓ_1, _, _ = Mℓ_and_Gℓ(ℓ, N_m)
        c = Mℓ_0*c + Mℓ_1*dL[Delta(ℓ)+1:Delta(ℓ+1)]
    end
    return c
end

"""
    largestNTerm(d::AbstractVector, N::Int)

Applies hard thresholding by retaining only the `N` largest coefficients in absolute value.

All other entries are set to zero.

# Arguments
- `d::AbstractVector`: Input coefficient vector.
- `N::Int`: Number of coefficients to retain (`1 ≤ N ≤ length(d)`).

# Output
- `Vector{Float64}`: Thresholded vector containing only the `N` largest terms.
"""
function largestNTerm(d::AbstractVector, N::Int)
    if N >= length(d) || N <= 0
        error("N must be in the range [1, length(d)].")
    end
    dc = copy(d)

    # Sort descending by absolute value
    d_sorted = sort(abs.(d), rev=true)

    # Identify the threshold value to keep exactly N elements
    smallest_kept_value = d_sorted[N] 
    dc[abs.(d[:]).<smallest_kept_value] .= 0.0
    return dc
end

"""
    relErrorTerm(d::AbstractVector, r::Real)

Applies relative error thresholding to a coefficient vector.

Zeroes out the smallest coefficients until the removed energy (in the ℓ₂-norm sense) 
reaches a relative value `r`.

# Arguments
- `d::AbstractVector`: Input coefficient vector.
- `r::Real`: Relative error tolerance (`0 < r < 1`).

# Output
- `Vector{Float64}`: Thresholded vector preserving the specified relative energy.
"""
function relErrorTerm(d::AbstractVector, r::Real)
    if r <=0 || r >= 1
        @warn "r should strictly be in the range (0, 1)."
    end
    max_error = r*norm(d,2)
    d_sorted_indices = sortperm(abs.(d))
    dc = copy(d)
    error_sqrt = 0.0

    for index in d_sorted_indices
        v = d[index] 
        error_sqrt += v^2 # Accumulate the energy of the "deleted" entries

        if sqrt(error_sqrt) <= max_error 
            dc[index] = 0.0
        else
            break
        end
    end
    return dc
end

"""
    compression_2D(Img_input::Matrix, ℓ_0::Int; N_m::N_m_and_Wavelet=N_1, threshold::Bool=false, r::Real=-1, N::Int=-1, type::Symbol=:iso, adapt_ref::Bool=false)

Performs two-dimensional wavelet compression on an input image.

Supports isotropic (`:iso`) and anisotropic (`:ani`) wavelet transforms and optional 
thresholding based on absolute or relative criteria.

# Arguments
- `Img_input::Matrix`: Input image (2D matrix) to compress.
- `ℓ_0::Int`: Lowest decomposition level.

# Keywords
- `N_m::N_m_and_Wavelet`: Wavelet model defining masks and level functions (default: `N_1`).
- `threshold::Bool`: Whether to perform coefficient thresholding (default: `false`).
- `r::Real`: Relative error tolerance (used if `threshold=true`).
- `N::Int`: Number of largest coefficients to retain (used if `threshold=true`).
- `type::Symbol`: Transform type; `:iso` for isotropic or `:ani` for anisotropic.
- `adapt_ref::Bool`: If `true`, returns both the reconstructed image and the thresholded coefficients.

# Output
- `Matrix` or `Tuple{Matrix, Matrix}`: The compressed image (and the coefficients if `adapt_ref=true`).

# Notes
- Prints the achieved compression ratio when thresholding is active.
- Returns the reconstructed image if `threshold=true`, otherwise the  
  low–resolution approximation at level `ℓ_0`.
"""
function compression_2D(Img_input::Matrix, ℓ_0::Int; N_m::N_m_and_Wavelet=N_1, threshold::Bool=false, r::Real=-1, N::Int=-1, type::Symbol=:iso, adapt_ref::Bool=false)
    Delta(x) = N_m.Delta_l(x)
    Nabla(x) = N_m.Nabla_l(x)
    Img = copy(Img_input)
    L = 0
    if N_m == N_2
        try 
            L = Int(log2(size(Img,1)+1))
        catch 
            error("Image size must be (2^L-1)×(2^L-1) for N_2, L Integer")
        end
    elseif N_m == N_2_boundary
        try 
            L = Int(log2(size(Img,1)-1))
        catch 
            error("Image size must be (2^L+1)×(2^L+1) for N_2 with boundary adaptations, L Integer")
        end
    else
        try 
            L = Int(log2(size(Img,1)))
        catch 
            error("Image size must be (2^L)×(2^L) for N_1, L Integer")
        end
    end
    
    # Inverse Transform (Decomposition)
    for ℓ in L-1:-1:ℓ_0
        _, _, Gℓ_0, Gℓ_1 = Mℓ_and_Gℓ(ℓ, N_m)
        if type == :iso
            a = Gℓ_0 * Img[1:Delta(ℓ+1), 1:Delta(ℓ+1)] * transpose(Gℓ_0) # Approximation
            h = Gℓ_0 * Img[1:Delta(ℓ+1), 1:Delta(ℓ+1)] * transpose(Gℓ_1) # Horizontal
            v = Gℓ_1 * Img[1:Delta(ℓ+1), 1:Delta(ℓ+1)] * transpose(Gℓ_0) # Vertical
            d = Gℓ_1 * Img[1:Delta(ℓ+1), 1:Delta(ℓ+1)] * transpose(Gℓ_1) # Diagonal
            Img_iso_temp = [a h; v d]
            Img[1:Delta(ℓ+1), 1:Delta(ℓ+1)] = Img_iso_temp
        elseif type == :ani
            Tℓ_inv = create_Tℓ_inv(Gℓ_0, Gℓ_1, ℓ, L, Delta)
            Img = Tℓ_inv*Img*transpose(Tℓ_inv)
        end
    end
 
    if threshold == true
        Img_vec = reshape(Img, :)
        if N != -1 && r != -1
            error("Choose only one method of thresholding (N or r).")
        elseif N != -1
            Img_vec = largestNTerm(Img_vec, N)
        elseif r != -1
            Img_vec = relErrorTerm(Img_vec, r)
        end
        Img = reshape(Img_vec, Delta(L), Delta(L))
        Img_coeefs = copy(Img)

        printstyled("Compression Ratio = $(1 - (count(!iszero, Img_vec))/(Delta(L)*Delta(L)))\n", color=:blue)

        # Forward Transform (Reconstruction)
        for ℓ in ℓ_0:1:(L-1)
            Mℓ_0, Mℓ_1, _, _ = Mℓ_and_Gℓ(ℓ, N_m)
            if type == :iso
                a = Mℓ_0 * Img[1:Delta(ℓ), 1:Delta(ℓ)] * transpose(Mℓ_0) # Approximation
                h = Mℓ_0 * Img[1:Delta(ℓ), Delta(ℓ)+1:Delta(ℓ+1)] * transpose(Mℓ_1) # Horizontal
                v = Mℓ_1 * Img[Delta(ℓ)+1:Delta(ℓ+1), 1:Delta(ℓ)] * transpose(Mℓ_0) # Vertical
                d = Mℓ_1 * Img[Delta(ℓ)+1:Delta(ℓ+1), Delta(ℓ)+1:Delta(ℓ+1)] * transpose(Mℓ_1) # Diagonal

                Img_iso_temp = a + h + v + d
                Img[1:Delta(ℓ+1), 1:Delta(ℓ+1)] = Img_iso_temp
            elseif type == :ani
                Tℓ = create_Tℓ(Mℓ_0, Mℓ_1, ℓ, L, Delta)
                Img = Tℓ*Img*transpose(Tℓ)  
            end
        end
        if adapt_ref == true
            return Img, Img_coeefs
        else
            return Img
        end
    else 
        return Img[1:Delta(ℓ_0), 1:Delta(ℓ_0)]./(2^(L-ℓ_0))
    end
end


"""
    heatmap_custom(Img::AbstractMatrix)

Generates a custom heatmap for matrix `Img` with a symmetric colour scale around zero.

Negative values are shown in blue tones, while positive values transition through a 
green-yellow gradient.

# Arguments
- `Img::AbstractMatrix`: Input matrix to visualize.

# Output
- Returns a `Plots.Plot` object with a customized colour scheme.
"""
function heatmap_custom(Img)
    dark_blue = RGB{Float64}(0, 0, 0.4)
    light_blue = RGB{Float64}(0.6, 0.8, 1.0)
    green = RGB{Float64}(0, 0.4, 0)
    yellow = RGB{Float64}(1, 1, 0)

    # Create gradient from green to yellow
    n_steps = 128
    green_yellow = range(green, stop=yellow, length=n_steps)
    blue_lightblue = range(dark_blue, stop=light_blue, length=n_steps)

    # Combine blue for <= 0 with gradient for > 0
    colors_combined = [collect(blue_lightblue); collect(green_yellow)]

    # Create colormap
    cmap = cgrad(colors_combined, [0.0; range(0.01, 1.0, length=n_steps)...])
   
    # Symmetric value range around 0, so white is exactly in the middle
    zmax = maximum(abs, Img)
    
    # Build plot
    p = heatmap(
        reverse(Img, dims=1), 
        xticks=(2 .^(7:11), string.(2 .^(7:11))), 
        yticks=(2 .^(7:11), string.(2 .^(7:11))), 
        size=(650, 600), 
        colorbar=true, 
        color=cmap, 
        clims=(-zmax, zmax), 
        right_margin=15Plots.mm
    )
    
    display(p)
    return p
end

### Flagged for CleanUp
"""
    example_compression(ℓ_0::Int; filepath::String)

Example routine demonstrating 2D wavelet compression on a sample TIFF image.

# Arguments
- `ℓ_0::Int`: The decomposition level for the compression algorithm.

# Keywords
- `filepath::String`: Path to the input `.tiff` image.
"""
function example_compression(ℓ_0; 
                             filepath::String = joinpath(@__DIR__, "src", "DataFolder", "-13_35_46_70.tiff"))
    Img2 = tiff_to_matrix(filepath) 
    start_x = 50; start_y= 2^10; LL=10
    Img2 = Img2[start_x:start_x+2^LL-1,end-2^LL+1-400:end-400]

    Img4 = compression_2D(Img2, ℓ_0, type=:ani)
    display(heatmap_custom(Img4))
end