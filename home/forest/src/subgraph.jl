# struct BaseGraph <: AbstractGraph
#     num_nodes::Int
#     num_edges::Int
#     total_pop::Int
#     #
#     pop_col::String
#     edge_weights::String
#     area_col::Union{String,Nothing}
#     node_border_col::Union{String,Nothing}
#     edge_perimeter_col::Union{String,Nothing}
#     oriented_nbrs_col::Union{String,Nothing}
#     #
#     populations::Array{Int,1}             # of length(num_nodes); seems redundant/not needed
#     adj_matrix::SparseMatrixCSC{Int,Int}  # not redundant, but I'm not sure if needed; also not quite an adjacency matrix
#     edge_src::Array{Int,1}                # of length(num_edges); seems redundant/not needed
#     edge_dst::Array{Int,1}                # of length(num_edges); seems redundant/not needed
#     neighbors::Array{Array{Int64,1},1}    # seems redundant/not needed
#     simple_graph::SimpleWeightedGraph     # the base SimpleWeightedGraph
#     node_attributes::Array{Dict{String,Any}}
#     edge_attributes::Dict{Set{Int},Dict{String,Any}}
# end

""""""
function cluster_base_graph(
    base_graph::BaseGraph, 
    clusters::Dict{String,T};
    node_key::String="nodes",
    cluster_key::String="clusters",
    copy_huh::Bool=true
) where T <: Any
    if copy_huh
        node_attributes = deepcopy(base_graph.node_attributes)
    else
        node_attributes = base_graph.node_attributes
    end
    sub_edges = Vector{Edge}(undef, 0)
    edges = SimpleWeightedGraphs.edges(base_graph.simple_graph)
    node_border_col = base_graph.node_border_col
    edge_perimeter_col = base_graph.edge_perimeter_col

    for edge in keys(base_graph.edge_attributes)
        n1, n2 = collect(edge)
        cluster_nodes = get_cluster(clusters[cluster_key], node_attributes[n1], 
                                    node_key)
        if in_cluster_huh(cluster_nodes, node_attributes[n2])
            e1, e2 = Edge(n1, n2), Edge(n2, n1) 
            if e1 in edges
                push!(sub_edges, e1)
            elseif e2 in edges
                push!(sub_edges, e2)
            end
        else
            edge_atr = base_graph.edge_attributes[Set([n1, n2])]
            edge_perim = edge_atr[edge_perimeter_col]
            node_attributes[n1][node_border_col] += edge_perim
            node_attributes[n2][node_border_col] += edge_perim
        end
    end
    simple_graph, vmap = induced_subgraph(base_graph.simple_graph, sub_edges)

    pop_col = base_graph.pop_col
    populations = [base_graph.node_attributes[ii][pop_col] for ii in vmap]
    total_pop = sum(populations)

    node_attributes = [node_attributes[ii] for ii in vmap]
    edge_attributes = Dict{Set{Int},Dict{String,Any}}()

    for n1 = 1:nv(simple_graph)
        for n2 in neighbors(simple_graph, n1)
            new_key = Set([n1, n2])
            old_key = Set([vmap[n1], vmap[n2]])
            edge_attributes[new_key] = base_graph.edge_attributes[old_key]
            weight = base_graph.simple_graph.weights[vmap[n1], vmap[n2]]
            simple_graph.weights[n1, n2] = weight
            simple_graph.weights[n2, n1] = weight
        end
    end

    return BaseGraph(
        nv(simple_graph),
        ne(simple_graph),
        total_pop,
        base_graph.pop_col,
        base_graph.edge_weights,
        base_graph.area_col,
        base_graph.node_border_col,
        base_graph.edge_perimeter_col,
        base_graph.oriented_nbrs_col,
        base_graph.mcd_col,
        simple_graph,
        node_attributes,
        edge_attributes
    )
end

""""""
function get_cluster(
    clusters::Vector,
    node_attributes::Dict{String,Any},
    node_key::String
)
    to_return = nothing
    for cluster in clusters
        nodes = cluster[node_key]
        for node in nodes
            in_node_huh = in_node(node, node_attributes)
            if in_node_huh
                if to_return != nothing
                    throw(
                        ArgumentError(
                            "Node with attributes ", node_attributes, " is in ",
                            "multiple clusters defined by ", clusters,
                        )
                    )
                else
                    to_return = nodes
                end
            end
        end
    end
    return to_return
end

function in_cluster_huh(
    cluster_nodes::Union{Nothing,Vector}, 
    node_attributes::Dict{String,Any}
)  
    if cluster_nodes == nothing
        return false
    end
    for node in cluster_nodes
        in_node_huh = in_node(node, node_attributes)
        if in_node_huh
            return true
        end
    end
    return false
end

function in_node(
    node::Dict{String, Any},
    node_attributes::Dict{String, Any}
)
    match_features = [node_attributes[key] == value for (key, value) in node]
    return !(false in match_features)
end