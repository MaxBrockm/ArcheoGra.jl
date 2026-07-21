# Adaptive Grid Refinement

For optimized memory usage and highly performant pathfinding, topographical data can be represented as an adaptive quadtree (leaves and hat nodes). Areas with high variance (e.g., mountains) are resolved finely, while uniform areas (e.g., plains) are grouped.

## Quadtree Generation

```@docs
build_G
adaptive_leaf_cells
get_A_j
get_adaptive_matrices
leaf_edge_masks
cut_leaves
leaf_cells_to_level_matrix
group_leaves_by_level
```

## Adaptive Adjacency Structures

```@docs
leaf_adjacency_matrix_directed
get_leaves_adjacency_matrix
get_weighted_adaptive_matrix
leaves_to_hat_nodes
hat_adjacency_matrix
get_weighted_hat_matrix
```

## Pipelines, Examples, and Plots

```@docs
plot_quadtree_leaves
run_adaptive_hat_pipeline
benchmark_adaptive_lcpa
example_adaptive
example_weighted_adaptive_matrix
```