# Advanced Clustering and Graph Neural Networks

This section documents advanced analytical methodologies applied to the networks, including various community detection algorithms, Markov Clustering (MCL), and integrated Graph Neural Network layers. The code for GNN-Clustering and DRA-clustering was kindly provided by Lukas Schmitz and based on his master's thesis: Schmitz, L. (2026). Multigrid Methoden für Graph-Laplace Matrizen. Master's thesis, University of Cologne. The DRA algorithm is based on (Napov, A and  Notay, Y, *An Efficient Multigrid Method for Graph Laplacian Systems II: Robust Aggregation*, SIAM Journal on Scientific Computing, 2017) and the GNN clustering on (Moore et. al., *Graph Neural Networks and Applied Linear Algebra*, SIAM review, 2025).

## Clustering Algorithms

```@docs
extract_clusters
identify_spectral_gap
spectral_kmeans
FluidC_clustering
FluidC_clustering_modified
MCL
R_MCL
restarted_MCL
restarted_R_MCL
normalized_cut
ClusterMethod
archaeological_clusterings
```

## Cluster Matching and Operations

```@docs
mcl_inflate_prune!
column_normalise!
build_coarse_adjacency
jaccard_distance_matrix
hungarian_match
match_clusters
```

## Graph Neural Network Components

```@docs
GraphData
GNNLayer
apply_layer
strength_aggregation_adj
strength_of_connection_gnn
```

## DRA Clustering approach

```@docs
DRA_QC_adj
DRA_QC_CE_adj
DRA_QC_CE_detailed_adj
```