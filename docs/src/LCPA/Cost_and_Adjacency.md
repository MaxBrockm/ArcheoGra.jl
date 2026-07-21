# Cost Functions and Adjacency Matrices

Before a path can be calculated, the topographical map must be converted into a weighted graph. The weights (costs) are determined by mathematical models of human movement.

## Cost Models

```@docs
Tobler_Hiking
KS_Hiking
adjusted_Hiking
f_LloberaSluckin
distance(u::Tuple, v::Tuple)
```

## Configuration Constants

```@docs
WeightConfig
Kondo_Average
Haversine_Config
```

## Adjacency Matrix Generation

```@docs
weighted_adjacency_matrix
weighted_adjacency_matrix_haversine
weighted_adjacency_matrix_combi
```