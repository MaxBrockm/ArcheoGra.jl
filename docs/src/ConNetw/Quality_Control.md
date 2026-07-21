# Quality Control and Validation

Functions to evaluate the structural integrity and quality of the generated clusterings and networks.

## General Quality Metrics

```@docs
modularity
silhouette_index
conductance
coverage
performance
```


## Internal Quality Functions

```@docs
_ncut_bisect
_dra_tentative_adj
_compute_gamma
_qc_extract
_bad_vertices_removal
_ext_weight
_int_weight
_criterion_16
_cholesky_quality_test
_build_AG_XG
_subgroup_extraction
_dra_qc_ce_internal_adj
build_coarse_graph
iterative_deep_clustering
```