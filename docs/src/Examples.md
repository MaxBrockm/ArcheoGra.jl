# Examples & Tutorials

This section provides practical examples of how to use the package for archaeological network generation, spectral analysis, spatial clustering, and Least Cost Path Analysis (LCPA). 

You can run these examples using the built-in sample datasets (`Setup_FullLBK.toml` and `Setup_Magdalenien_new.toml`) or adapt them for your own custom data and configuration files.

---

## 1. Archaeological Networks (Haversine & Material Flows)

These functions construct spatial and material-weighted networks based on Haversine distance calculations between archaeological sites.

### Using Built-in Datasets
To generate networks and analyze features for the package-provided LBK and Magdalenian datasets, use `example_networks_given_data_fullLBK_magda_new`.

```julia
using ArcheoGra 

# Run with default settings (plots generated, no optional feature analysis printed)
results = example_networks_given_data_fullLBK_magda_new()

# Run with feature analysis enabled and custom output directory
results = example_networks_given_data_fullLBK_magda_new(
    output_dir = "./my_outputs",
    analyze_features = true,
    overwrite = true
)

```

```@docs
example_networks_given_data_fullLBK_magda_new
```
**Key Parameters:**
* `output_dir::String` (default: `pwd()`): Directory where output subfolders and files will be saved.
* `overwrite::Bool` (default: `false`): Overwrite existing output directories if they exist.
* `plot_compare::Bool` (default: `true`): Whether to generate comparison plots.
* `plot_results::Bool` (default: `true`): Whether to plot final network results.
* `analyze_features::Bool` (default: `false`): If `true`, prints network feature summaries and runs `is_average`.


### Using Your Own Data

To run the network pipeline on your own custom dataset, create a TOML configuration file (analogous to the ones in src/DataFolder) pointing to your site and material tables, then call `example_networks_from_own_data`.

```julia
# Run with a custom TOML configuration file
graphs = example_networks_from_own_data(
    "path/to/my_custom_config.toml";
    incl_sources = true,
    double_distance = false,
    output_dir = "./my_outputs",
    analyze_features = true
)
```
```@docs
example_networks_from_own_data
```

**Setup Checklist for Custom Data:**
1. Create a TOML configuration file defining your archaeological tables and sources.
2. Ensure file paths inside the TOML correctly point to your local data files.
3. Pass the path of your TOML file to the function.

### Detailed TOML Guide (Field-by-Field)

This section explains exactly how to build TOML setup files for ConNetw-based examples.

#### 1. General Structure

Your TOML can contain one or multiple blocks. Each block represents one archaeological input table (for one time period).

```toml
[file_Magda]
filePath = "src\\DataFolder\\Magda_GridPts_clean_comma.csv"
file_type = "csv"
delim = ","
material_cols_start = 4
material_cols_end = 68
material_age = "Magdanew"
material_distance = 200
grid_id = 1
long_ind = 2
lat_ind = 3
source_data_path = "src\\DataFolder\\Magda_sources_to_grdpts_corr_comma.csv"
source_data_delim = ","
source_data_name_col = 1
source_data_long_col = 3
source_data_lat_col = 4
source_data_mult = 1
source_grid_id = 2
```

You can add additional blocks, for example:

```toml
[file_oldLBK]
filePath = "src\\DataFolder\\oLBK_ds02_belegt.txt"
file_type = "txt"
delim = "\t"
material_cols_start = 6
material_cols_end = 9
material_age = "oLBK"
material_distance = 350
grid_id = 1
long_ind = 2
lat_ind = 3
```

#### 2. Path Resolution (Important)

Relative paths in TOML are resolved relative to the folder of the TOML file itself.

That means:
1. If your TOML is in a project folder, keep paths relative to that TOML.
2. You do not need to depend on the current working directory.
3. Absolute paths are also allowed.

#### 3. Required Archaeological Fields

For each block, these fields are required for archaeological data loading:

* `filePath`
* `file_type`
* `material_cols_start`
* `material_cols_end`
* `material_age`
* `material_distance`
* `grid_id`
* `long_ind`
* `lat_ind`

#### 4. Optional Archaeological Fields

These are optional and only needed if present in your table:

* `sheetname` (for xlsx)
* `delim` (for txt or csv)
* `name_ind`
* `class_ind`

#### 5. Allowed file_type Values

For archaeological site tables:
1. csv
2. txt
3. xlsx

For source tables:
1. txt
2. xlsx

Note: `source_file_type` defaults to txt when omitted.

#### 6. Delimiters: What Is Allowed and How It Is Interpreted

In the parser, `delim` and `source_data_delim` are read as strings, and only the first character is used.

Valid and common values:
1. "," for comma-separated files
2. "\t" for tab-separated files
3. ";" for semicolon-separated files
4. "|" for pipe-separated files

Recommendations:
* Always pass delimiter values as a string of length 1, or as "\t" for tabs.
* Avoid multi-character delimiters because only the first character is used.

#### 7. material_distance Rules

You can define material_distance in two ways:

Single number:
* The same value is used for all material columns in the range `material_cols_start:material_cols_end`.

Vector/array of numbers:
* Must match the number of material columns exactly.

Example for vector form:

```toml
material_distance = [200, 220, 180, 250]
```

#### 8. Source-Data Fields (When incl_sources=true)

If you want source-based models, include these fields in each block:

* `source_data_path`
* `source_data_name_col`
* `source_data_long_col`
* `source_data_lat_col`
* `source_grid_id`

Optional source fields:
* `source_data_delim` (default behavior depends on `source_file_type`)
* `source_file_type` (txt or xlsx)
* `source_sheetname` (required when `source_file_type` is xlsx)
* `source_data_mult`

Behavior notes:
* `source_data_mult` defaults to 1.0 if omitted.
* Source material names are internally combined with `material_age` to match material IDs.

#### 9. Column Index Convention

All indices in TOML are 1-based (Julia convention).

Examples:
* `long_ind = 2` means longitude is in the second column.
* `material_cols_start = 4` and `material_cols_end = 68` means material values are read from columns 4 through 68.

#### 10. Typical Pitfalls and How to Avoid Them

1 Wrong delimiter
* Symptom: columns are shifted or parsing fails.
* Fix: verify `delim`/`source_data_delim` and match the actual file format.

2 Wrong column indices
* Symptom: wrong coordinates or empty material sets.
* Fix: test with a small file first and verify `long_ind`, `lat_ind`, and material column range.

3 Mismatch in material_distance vector length
* Symptom: runtime conversion or indexing errors.
* Fix: ensure vector length equals number of material columns.

4 Missing source fields while `incl_sources=true`
* Symptom: source graph creation fails.
* Fix: add required `source_data_*` fields or run with `incl_sources=false`.

5 xlsx source without `source_sheetname`
* Symptom: sheet lookup fails.
* Fix: add `source_sheetname` explicitly.

#### 11. Minimal Valid Templates

Minimal site-only template:

```toml
[file_example]
filePath = "data\\sites.csv"
file_type = "csv"
delim = ","
material_cols_start = 5
material_cols_end = 12
material_age = "ExamplePeriod"
material_distance = 300
grid_id = 1
long_ind = 2
lat_ind = 3
```

Minimal site-plus-source template:

```toml
[file_example]
filePath = "data\\sites.csv"
file_type = "csv"
delim = ","
material_cols_start = 5
material_cols_end = 12
material_age = "ExamplePeriod"
material_distance = 300
grid_id = 1
long_ind = 2
lat_ind = 3
source_data_path = "data\\sources.csv"
source_data_delim = ","
source_data_name_col = 1
source_data_long_col = 2
source_data_lat_col = 3
source_grid_id = 2
```

## 2. Spectral Analysis

Spectral analysis functions generate spectral plots and evaluate graph structures across multiple runs or random seeds.

### Using Built-in Datasets
```julia
# Run spectral plots for package datasets
out_dir = example_spectral_plots_given_data(
    num_graphs = 100,
    seed = 42,
    output_dir = "./my_outputs"
)
```
```@docs
example_spectral_plots_given_data
```
### Using Your Own Data
```julia
# Run spectral plots with your custom TOML configuration
out_dir = example_spectral_plots_from_own_data(
    "path/to/my_custom_config.toml";
    num_graphs = 50,
    seed = 123,
    incl_sources = true,
    double_distance = false
)
```
```@docs
example_spectral_plots_from_own_data
```
## 3. Spatial Clustering

Clustering functions help identify regional groupings and structural clusters within your archaeological network data.

### Using Built-in Datasets
```julia
# Run clustering on package-provided datasets
cluster_dirs = example_cluster_results_given_data(
    num_clusters = 0,        # 0 means automatic determination
    max_clusters_param = 10,
    plot_map = true,
    seed = 42
)
```
```@docs
example_cluster_results_given_data
```
### Using Your Own Data
```julia
# Run clustering using your custom TOML setup
cluster_dir = example_cluster_results_from_own_data(
    "path/to/my_custom_config.toml";
    incl_sources = true,
    double_distance = false,
    num_clusters = 3,
    plot_map = true
)
```
```@docs
example_cluster_results_from_own_data
```

## 4. Least Cost Path Analysis

Least Cost Path Analysis computes optimal movement pathways across topographic cost surfaces (provided as georeferenced TIFF files). You can choose between **Uniform-grid** pipelines and **Adaptive** wavelet-compressed pipelines.

### Uniform-Grid LCPA
* **Using built-in TIFF file:**
  ```julia
  result = example_site_to_site_lcpa_uniform_given_file(
      start_site = (-5.4, 38.6),
      target_site = (10.0, 45.0),
      L_ref = 10,
      margin_ratio = 0.5
  )
  ```
* **Using your own TIFF file:**
  ```julia
  result = example_site_to_site_lcpa_uniform_own_file(
      "path/to/elevation_data.tiff";
      start_site = (x_start, y_start),
      target_site = (x_target, y_target),
      L_ref = 10
  )
  ```

```@docs
example_site_to_site_lcpa_uniform_given_file
example_site_to_site_lcpa_uniform_own_file
```

  ### Adaptive LCPA (Wavelet-Based)
Adaptive LCPA uses wavelet compression (`:N_1` or `:N_2`) to speed up computations on large landscape grids while controlling error tolerances (`r_tol`).

* **Using built-in TIFF file:**
  ```julia
  result = example_site_to_site_lcpa_adaptive_given_file(
      wavelet_choice = :N_1,    # Choose explicitly: :N_1 or :N_2
      r_tol = 0.05,
      start_site = (-5.4, 38.6),
      target_site = (10.0, 45.0)
  )
  ```
* **Using your own TIFF file:**
  ```julia
  result = example_site_to_site_lcpa_adaptive_own_file(
      "path/to/elevation_data.tiff";
      start_site = (x_start, y_start),
      target_site = (x_target, y_target),
      wavelet_choice = :N_2,
      r_tol = 0.05
  )
  ```
```@docs
example_site_to_site_lcpa_adaptive_given_file
example_site_to_site_lcpa_adaptive_own_file
```

  ### Benchmarking Adaptive LCPA
If you want to evaluate different compression tolerances (`r_values`), vertex budget sizes (`N_values`), and wavelet types (`:N_1`, `:N_2`) against each other, use the benchmark wrapper:

```julia
benchmark_results = example_benchmark_adaptive_lcpa(
    "path/to/elevation_data.tiff";
    start_site = (-5.4, 38.6),
    target_site = (10.0, 45.0),
    wavelet_choices = [:N_1, :N_2],
    r_values = [0.01, 0.05, 0.1],
    N_values = [10000, 50000]
)
```
*This automatically executes the benchmark suite, prints/returns the results, and writes a summary CSV file to the output directory.*
```@docs
example_benchmark_adaptive_lcpa
```