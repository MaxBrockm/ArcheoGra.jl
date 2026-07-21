# Wavelet-based Data Compression

To handle massive topographical datasets spanning large geographic areas, `ArcheoGra` utilizes wavelet transformations to compress spatial data while preserving relevant topological features.

## Transformation Matrices and Coefficients

```@docs
coeffi_mat_Mℓ_0
coeffi_mat_Mℓ_1
N_m_and_Wavelet
Mℓ_and_Gℓ
N_1
N_2
N_2_boundary
create_Tℓ_inv
create_Tℓ
```

## Compression and Error Analysis

```@docs
compression_1D
compression_2D
largestNTerm
relErrorTerm
```

## Visualization and Examples

```@docs
heatmap_custom
example_compression
```