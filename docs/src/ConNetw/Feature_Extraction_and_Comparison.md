# Feature Extraction and Graph Comparison

To analyse the stochastically generated networks robustly, multiple metrics and feature extraction routines are provided. These tools allow for rigorous comparison between different network iterations or distinct time periods. The feature comparison approach is based on The NetSimile approach from [Berlingerio, M et. al., *NetSimile: A Scalable Approach to Size-Independent Network Similarity*, 2012](https://arxiv.org/pdf/1209.2684). 

## Network Feature Extraction

```@docs
feature_extraction
feature_aggregation
network_features
extract_degree_distribution
extract_clustering_coefficient
extract_egonet_features
extract_path_features
extract_density
```

## Distance and Statistical Comparison

```@docs
matrix_spectral_distance
effective_resistance
canberra_distance
SRGG
mahalanobis_distance
mahalanobis_pinv
standardized_distance_squared
empirical_quantile
empirical_cdf_interpolated
is_average
graph_sensitivity
create_graph_statistic_Haversine
```

## High-Level Comparison Pipelines

```@docs
produce_results
compare_all_is_connected
compare_all_connect_nodes
compare_all_is_connected_sources
compare_all_connect_nodes_sources
get_matrices
```