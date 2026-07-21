# Data Management, Maps, and I/O

Topographical data is usually provided as raster graphics. This section documents the functions used to read these `.tiff` files, handle water bodies, and convert between spatial coordinates and internal matrix indices.

## Input and Output

```@docs
tiff_to_matrix
write_Output
```

## Water Masks and Boundaries

```@docs
WaterLandBool
adaptive_WaterLandBool
```

## Indices and Coordinate Numbering

Efficient matrix operations require converting 2D spatial coordinates into 1D lexicographical indices and vice versa.

```@docs
characteristic
coords2lexi
coords2leaf_index
ind2coords(ij::Tuple{Int,Int}, Topdata::AbstractLCPAData)
ind2coords(leaf::Tuple{Int,Int,Int}, TopData::AbstractLCPAData, L_geom::Int)
lexiMat2Vec
get_xy
lexi2ind
```