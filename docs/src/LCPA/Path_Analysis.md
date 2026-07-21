# Path Analysis and Routing

This section contains the core algorithmic implementations to compute the shortest paths across the weighted adjacency matrices, including Dijkstra's algorithm and large-scale parallelized distance computations.

## Shortest Path Algorithms

```@docs
neighbours
delConnections
compress
get_path
dijkstra_shortest_path
```

## Multi-Pair LCPA Pipelines

```@docs
all_pairs_LCPA
all_pairs_LCPA_GMT
all_pairs_LCPA_adaptive
all_pairs_LCPA_adaptive_GMT
extract_quadratic_map
```