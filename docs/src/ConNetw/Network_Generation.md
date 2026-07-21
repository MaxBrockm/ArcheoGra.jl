# Network Generation and Connectivity

These functions handle the core logic of translating archaeological sites into a connected graph. This includes applying distance measures (like Haversine), calculating connection probabilities via Rayleigh functions, and ensuring the final stochastic graph forms a single connected component.

## Distance and Basic Connectivity

```@docs
distance(u::ArchDataEntry, v::ArchDataEntry)
distance_allocation(data::AbstractArray{ArchDataEntry})
distance_allocation(data::AbstractArray{ArchDataEntry}, source_data::Source_DataSet)
distance_allocation(source_data::Source_DataSet)
adjacency_data(Material_data::MaterialDataList, data::Arch_DataSet, time_period::String)
adjacency_data(Material_data::MaterialDataList, prep::PreparedTimePeriod)
adjacency_sources(data::Arch_DataSet, SourceData::Source_DataSet, MaterialData::MaterialDataList, time_period::String)
adjacency_sources(prep::PreparedTimePeriod, SourceData::Source_DataSet, MaterialData::MaterialDataList)
RayleighSISO
connect_nodes(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, material_weight_matrices::Dict, relevant_sites::Vector{ArchDataEntry}; symmetric::Bool = true)
connect_nodes(start_node::Int, end_node::Int, distance_matrix::AbstractArray{Float64}, symmetric::Bool = true; distance_matrix_back::AbstractArray{Float64} = distance_matrix)
```

## Haversine Specific Networks

```@docs
Haversine_connect_nodes
Haversine_connect_nodes_sources
Haversine_is_connected
Haversine_is_connected_sources
Haversine_is_connected_source2source
```

## Connected Components Assurance

```@docs
components
connect_components!
custom_argmin
```