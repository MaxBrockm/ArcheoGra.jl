# Connect Network Documentation

This submodule provides the foundational tools to construct networks from archaeological datasets. It handles the initial setup, global parameters, and fundamental data preparation required to model prehistoric mobility and social connectivity.

## Global Variables and Setup

The following global constants define the standard parameters for distance and raw material connectivity.

```@docs
maximum_distance
material_multiplier
setup_global_vars
```

## Project Initialization and Data Preparation

```@docs
prepare_network_data
get_sites_for_time
write_results
_output_path
```