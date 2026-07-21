"""
    components(adj_mat::AbstractArray{<:Number})

Function to determine all components of given graph. Function is an adjusted version of the function "connect_components" from the Graphs.jl package.
Changes are made, such that the function can determine the components just from the adjacency matrix, instaed of creating the corresponding graph.
The output of this function is a vector of labels, where the i-th entry is the label of the component of the i-th vertex. Second output is a list of 
vectors, where each vector consists of the vertices of a corresponding component of the graph.  

# Arguments:
- `adj_mat::AbstractArray{<:Number}`: Adjacency matrix
"""
function components(adj_mat::AbstractArray{<:Number})
    labels = zeros(Int, size(adj_mat,1))
    component_list = []
    for vertex in axes(adj_mat,1) ## Find clusters by going through all vertices
        labels[vertex] != 0 && continue
        labels[vertex] = vertex
        Q = [vertex]
        component = [vertex]
        while !isempty(Q)
            src = popfirst!(Q) # find connected vertices of first vertex in list 
            for neighbour in findall(x->x!=0,adj_mat[:,src])
                if labels[neighbour] == 0 ## if vertex hasnt already been assigned to cluster, put in list 
                    labels[neighbour] = vertex
                    push!(Q, neighbour); push!(component, neighbour)
                end
            end
        end
        push!(component_list, component)
    end
    return labels, component_list
end

"""
    connect_components!(adjacency_matrix_daten::AbstractArray{<:Number}, labels::AbstractVector{Int}, component_list::Vector{Any}, dist_matrix::AbstractArray{<:Number})

Function to connect all components from a graph. Function uses only the adjacency matrix as a graph description. The component distinction comes form [`components(adj_mat::AbstractArray{<:Number})`](@ref).
Components are connect, such that a edge is drawn between vertices of two different components with the shortest distance between them. Distances are taken from distance calculations in [`distance_allocation(data::AbstractArray{ArchDataEntry})`](@ref).

# Arguements:
- `adjacency_matrix_daten::AbstractArray{<:Number}`: Adjacency matrix of graph
- `labels::AbstractVector{Int}`: Vector whose i-th entry is the component label of the i-th vertex.
- `component_list::Vector{Any}`: List of vectors, where each vector contains the vertices that are contained in the corresponding component. 
- `dist_matrix::AbstractArray{<:Number}`: Matrix with distances between vertices
"""
function connect_components!(adjacency_matrix_daten::AbstractArray{<:Number}, labels::AbstractVector{Int}, component_list::Vector{Any}, dist_matrix::AbstractArray{<:Number})
    total_nodes = size(adjacency_matrix_daten,1)
    ## start with first cluster and add every other cluster to it 
    first_comp_pointer = zeros(Bool, total_nodes); first_comp_pointer[component_list[1]].=1
    while labels != ones(Int, total_nodes)
        ## find shortest connection between connected sites (cluster 1) and all other vertices
        tobeconnected = custom_argmin(dist_matrix.* (.!first_comp_pointer .* first_comp_pointer'))
        adjacency_matrix_daten[tobeconnected...] = 1; adjacency_matrix_daten[reverse(tobeconnected)...] = 1
        connected_comp = findfirst(x->tobeconnected[1]in x, component_list) ## which cluster has been connected
        labels[labels.==labels[tobeconnected[1]]] .= 1
        first_comp_pointer[component_list[connected_comp]] .= 1 
    end
    return 
end

"""
    custom_argmin(A::AbstractArray{<:Number})

Function to determine the smallest nonzero value of an array

# Arguments:
- `A::AbstractArray{<:Number}`: Array from which the smallest nonzero value is to be determined. 
"""
function custom_argmin(A::AbstractArray{<:Number})
    si,sj = 1,1; min_val = Inf
    for i in axes(A,1), j in axes(A,2) 
        (0 < A[i,j] < min_val) ? (si,sj, min_val) = (i,j, A[i,j]) : 0 
    end
    return [si,sj]
end