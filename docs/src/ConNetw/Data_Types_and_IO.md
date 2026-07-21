# Data Structures and Input/Output

This section details the robust data structures used to store archaeological and source data, alongside the functions required to parse datasets from standard `.toml`, `.csv`, and `.xlsx` formats.

## Core Data Structures

```@docs
ArchDataEntry
SourceDataEntry
Source_DataSet
Arch_DataSet
MaterialEntry
MaterialDataList
file_props
PreparedTimePeriod
```

## I/O and Data Management

```@docs
read_data_TOML
read_data
read_source_data
add_data!
add_material!
num_entries(dataset::Arch_DataSet)
num_entries(dataset::Arch_DataSet, property::Symbol, keyvalue::Any)
get_coordinates
check_presence
get_material_id
```