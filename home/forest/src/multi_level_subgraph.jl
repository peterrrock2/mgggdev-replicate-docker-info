struct MultiLevelSubGraph #<: AbstractGraph
    node_set::Dict{Tuple{Vararg{String}}, Any}
    graphs_by_level::Vector{Dict{Tuple{Vararg{String}}, SimpleWeightedGraph}}
    vmaps::Vector{Dict{Tuple{Vararg{String}}, Vector{Int}}}
    modified_node_populations::Vector{Dict{Tuple{Vararg{String}}, Int}}
    modified_nbr_weights::Dict{Set{Tuple{Vararg{String}}}, Any}
    parent::MultiLevelGraph
end

""""""
function get_total_graph_population(
    subgraph::MultiLevelSubGraph
)::Int
    return subgraph.modified_node_populations[1][()] # by construction
end

""""""
function get_node_population(
    subgraph::MultiLevelSubGraph,
    node::Tuple
)::Int
    level = length(node)
    modified_node_populations = subgraph.modified_node_populations
    num_levels = subgraph.parent.num_levels
    if level < num_levels && haskey(modified_node_populations[level+1], node)
        pop = modified_node_populations[level+1][node]
    else
        graph = subgraph.parent.graphs_by_level[level]
        id = subgraph.parent.partition_to_ids[level][node]
        pop_col = graph.pop_col
        pop = graph.node_attributes[id][pop_col]
    end
    return pop
end

""""""
function sum_node_data(
    graph::MultiLevelGraph,
    node_set::Dict{Tuple{Vararg{String}}, Any},
    column::String,
    level::Int=1
)
    value = 0
    node_attributes = graph.graphs_by_level[level].node_attributes
    for node in keys(node_set)
        if intact_huh(node_set[node])
            node_id = graph.partition_to_ids[level][node]
            value += node_attributes[node_id][column]
        else
            value += sum_node_data(graph, node_set[node], column, level+1)
        end
    end
    return value
end

""""""
function sum_node_data(
    subgraph::MultiLevelSubGraph,
    column::String,
    level::Int=1
)
    return sum_node_data(subgraph.parent, subgraph.node_set, column)
end

""""""
@inline function intersects_huh(subgraph::MultiLevelSubGraph, node)
    return intersects_huh(subgraph.node_set, node)
end

""""""
@inline function strictly_contained_huh(subgraph::MultiLevelSubGraph, node)
    return strictly_contained_huh(subgraph.node_set, node)
end

""""""
@inline function get_node_set(subgraph::MultiLevelSubGraph, node)
    return get_node_set(subgraph.node_set, node)
end

""""""
function get_node_population(
    subgraph::MultiLevelSubGraph,
    level::Int,
    id::Int
)::Int
    modified_node_populations = subgraph.modified_node_populations
    num_levels = subgraph.parent.num_levels
    node = subgraph.parent.id_to_partitions[level][id]
    if haskey(modified_node_populations[level+1], node) && level < num_levels
        pop = modified_node_populations[level+1][node]
    else
        graph = subgraph.parent.graphs_by_level[level]
        pop_col = graph.pop_col
        pop = graph.node_attributes[id][pop_col]
    end
    return pop
end

""""""
function get_simple_subgraph(
    subgraph::MultiLevelSubGraph,
    level::Int,
    key::Tuple
)
    if haskey(subgraph.graphs_by_level[level], key)
        simple_sub_graph = subgraph.graphs_by_level[level][key]
        vmap = subgraph.vmaps[level][key]
        return simple_sub_graph, vmap
    else
        graph = subgraph.parent
        level = length(key)
        index = graph.partition_to_ids[level][key]
        fine_graph = graph.coarse_to_fine_graphs[level][index]
        return fine_graph.graph, fine_graph.vmap
    end
end

""""""
function get_edge_key(node::Tuple, nbr::Tuple)
    ii = 1
    sz = 0
    for ii = 1:min(length(node), length(nbr))
        eq = (node[ii] == nbr[ii])
        sz += Int(eq)
        if !eq
            break
        end
    end
    return node[1:sz]
end

""""""
function build_subgraphs_by_level!(
    graphs_by_level::Vector{Dict{Tuple{Vararg{String}}, SimpleWeightedGraph}},
    vmaps::Vector{Dict{Tuple{Vararg{String}}, Vector{Int}}},
    graph::MultiLevelGraph,
    node_set::Dict{Tuple{Vararg{String}},Any},
    level::Int = 1,
    key::Tuple = ()
)
    @assert level <= graph.num_levels

    sub_nodes = Vector{Int}(undef, 0)
    for node in keys(node_set)
        push!(sub_nodes, graph.partition_to_ids[level][node])
    end

    if level == 1
        simple_graph = graph.graphs_by_level[level].simple_graph
        subgraph, vmap = induced_subgraph(simple_graph, sub_nodes)
    else
        key_id = graph.partition_to_ids[level-1][key]
        sub_simple_graph = graph.coarse_to_fine_graphs[level-1][key_id]
        simple_graph = sub_simple_graph.graph
        vmap = sub_simple_graph.vmap
        sub_nodes = [id for id = 1:nv(simple_graph)
                     if vmap[id] in sub_nodes] # not ideal; will be called a lot
        subgraph, svmap = induced_subgraph(simple_graph, sub_nodes)
        # map sub partition to global node index:
        vmap = [vmap[v] for v in svmap]
    end

    graphs_by_level[level][key] = subgraph
    vmaps[level][key] = vmap

    for (node, sub_nodes) in node_set
        if sub_nodes == nothing || length(sub_nodes) == 0
            continue
        end
        build_subgraphs_by_level!(graphs_by_level, vmaps, graph,
                                  sub_nodes, level+1, node)
    end
end

""""""
function get_populations!(
    modified_node_populations::Vector{Dict{Tuple{Vararg{String}}, Int}},
    graphs_by_level::Vector{Dict{Tuple{Vararg{String}}, SimpleWeightedGraph}},
    vvmaps::Vector{Dict{Tuple{Vararg{String}}, Vector{Int}}},
    graph::MultiLevelGraph,
    level::Int=1,
    key::Tuple=()
)
    num_levels = graph.num_levels
    levels = graph.levels
    pop_col = graph.graphs_by_level[1].pop_col

    modified_node_populations[level][key] = 0

    fine_graph = graph.graphs_by_level[level]
    fine_node_attributes = fine_graph.node_attributes

    for node_id in vvmaps[level][key]
        n_atr = fine_node_attributes[node_id]
        if level == num_levels
            modified_node_populations[level][key] += n_atr[pop_col]
        else
            node = Tuple([n_atr[l] for l in levels[1:level]])
            if haskey(graphs_by_level[level+1], node)
                get_populations!(modified_node_populations, graphs_by_level,
                                 vvmaps, graph, level+1, node)
                modified_node_populations[level][key] +=
                    modified_node_populations[level+1][node]
            else
                modified_node_populations[level][key] += n_atr[pop_col]
            end
        end
    end
end

""""""
function get_edge_weight(
    graph::MultiLevelGraph,
    node_set::Union{Dict{Tuple{Vararg{String}}, Any}, Nothing},
    nbr_set::Union{Dict{Tuple{Vararg{String}}, Any}, Nothing},
    node::Tuple{Vararg{String}},
    nbr::Tuple{Vararg{String}}
)
    # fine_neighbors::Vector{Vector{Dict{Int,Any}}}   # for each level, for each node, have fine nbrs
    new_weight = 0
    if !intact_huh(node_set)
        fine_nodes_ids_near_nbr = get_fine_neighbors(graph, nbr, node)
        for fine_node_id_near_nbr in fine_nodes_ids_near_nbr
            fine_node_near_nbr = graph.id_to_partitions[length(node)+1][
                fine_node_id_near_nbr]
            if fine_node_near_nbr ∉ keys(node_set)
                continue
            end
            new_weight += get_edge_weight(graph,
                                          node_set[fine_node_near_nbr],
                                          nbr_set, fine_node_near_nbr, nbr)
        end
    elseif intact_huh(nbr_set) # && intact_huh(node_set)
        edge_key = Set([node, nbr])
        new_weight = edge_weight(graph, edge_key)
    else # !intact_huh(nbr_set) && intact_huh(node_set)
        fine_nbr_ids_near_node = get_fine_neighbors(graph, node, nbr)
        for fine_nbr_id_near_node in fine_nbr_ids_near_node
            fine_nbr_near_node = graph.id_to_partitions[length(nbr)+1][
                fine_nbr_id_near_node]
            if fine_nbr_near_node ∉ keys(nbr_set)
                continue
            end
            new_weight += get_edge_weight(graph, node_set,
                                          nbr_set[fine_nbr_near_node],
                                          node, fine_nbr_near_node)
        end
    end
    return new_weight
end

""""""
function get_edge_weight(
    subgraph::MultiLevelSubGraph,
    node_set::Union{Dict{Tuple{Vararg{String}}, Any}, Nothing},
    nbr_set::Union{Dict{Tuple{Vararg{String}}, Any}, Nothing},
    node::Tuple{Vararg{String}},
    nbr::Tuple{Vararg{String}}
)
    edge_key = Set([node, nbr])
    if haskey(subgraph.modified_nbr_weights, edge_key)
        return subgraph.modified_nbr_weights[edge_key]
    # elseif length(node)!=length(nbr)
    #     return subgraph.parent.mixed_nbr_weights[edge_key]
    end
    graph = subgraph.parent

    return get_edge_weight(graph, node_set, nbr_set, node, nbr)
end

""""""
function get_edge_weight(
    subgraph::MultiLevelSubGraph,
    node::Tuple{Vararg{String}},
    nbr::Tuple{Vararg{String}}
)
    node_set = get_node_set(subgraph, node)
    nbr_set = get_node_set(subgraph, nbr)
    return get_edge_weight(subgraph, node_set, nbr_set, node, nbr)
end

""""""
function modify_edge_weights!(
    modified_nbr_weights::Dict{Set{Tuple{Vararg{String}}}, Any},
    graph::MultiLevelGraph,
    node_set::Union{Dict{Tuple{Vararg{String}}, Any}, Nothing},
    nbr_set::Union{Dict{Tuple{Vararg{String}}, Any}, Nothing},
    node::Tuple{Vararg{String}},
    nbr::Tuple{Vararg{String}}
)
    # @show "recursing on modifying edge weights", node, nbr
    # fine_neighbors::Vector{Vector{Dict{Int,Any}}}   # for each level, for each node, have fine nbrs
    new_weight = 0
    if !intact_huh(node_set)
        # @show "branch 1", node, nbr
        fine_nodes_ids_near_nbr = get_fine_neighbors(graph, nbr, node)
        # @show fine_nodes_ids_near_nbr
        for fine_node_id_near_nbr in fine_nodes_ids_near_nbr
            fine_node_near_nbr = graph.id_to_partitions[length(node)+1][
                fine_node_id_near_nbr]
            if fine_node_near_nbr ∉ keys(node_set)
                continue
            end
            weight = modify_edge_weights!(modified_nbr_weights, graph,
                                          node_set[fine_node_near_nbr], nbr_set,
                                          fine_node_near_nbr, nbr)
            new_weight += weight
            # @show "branch 1 -- back", node, fine_node_near_nbr, weight, new_weight
            edge_key = Set([fine_node_near_nbr, nbr])
            if edge_weight(graph, edge_key) != weight
                modified_nbr_weights[edge_key] = weight
            end
        end
    elseif intact_huh(nbr_set) # && intact_huh(node_set)
        edge_key = Set([node, nbr])
        new_weight = edge_weight(graph, edge_key)
        # @show "branch 2", node, nbr, new_weight
    else # !intact_huh(nbr_set) && intact_huh(node_set)
        # @show "branch 3", node, nbr
        fine_nbr_ids_near_node = get_fine_neighbors(graph, node, nbr)
        # @show fine_nbr_ids_near_node
        for fine_nbr_id_near_node in fine_nbr_ids_near_node
            fine_nbr_near_node = graph.id_to_partitions[length(nbr)+1][
                fine_nbr_id_near_node]
            if fine_nbr_near_node ∉ keys(nbr_set)
                continue
            end
            weight = modify_edge_weights!(modified_nbr_weights, graph,
                                          node_set, nbr_set[fine_nbr_near_node],
                                          node, fine_nbr_near_node)
            new_weight += weight
            edge_key = Set([node, fine_nbr_near_node])
            # @show "branch 3 -- back", node, fine_nbr_near_node, weight, new_weight
            if edge_weight(graph, edge_key) != weight
                modified_nbr_weights[edge_key] = weight
            end
        end
    end
    return new_weight
end

""""""
function modify_edge_weights!(
    modified_nbr_weights::Dict{Set{Tuple{Vararg{String}}}, Any},
    node_set::Dict{Tuple{Vararg{String}}, Any},
    graphs_by_level::Vector{Dict{Tuple{Vararg{String}}, SimpleWeightedGraph}},
    vvmaps::Vector{Dict{Tuple{Vararg{String}}, Vector{Int}}},
    graph::MultiLevelGraph,
    level::Int=1,
    key::Tuple=()
)
    num_levels = graph.num_levels
    levels = graph.levels

    full_graph = graph.graphs_by_level[level]
    node_attributes = full_graph.node_attributes
    sub_graph = graphs_by_level[level][key]
    sub_vmap = vvmaps[level][key]

    for (ii, node_id) in enumerate(sub_vmap)
        name = graph.id_to_partitions[level][node_id]
        n_atr = node_attributes[node_id]
        node = Tuple([n_atr[l] for l in levels[1:level]])
        # @show "In modify_edge_weights", node, level
        if intact_huh(node_set[node])
            continue
        end
        nbrs = copy(neighbors(sub_graph, ii))
        for jj in nbrs
            nbr_id = sub_vmap[jj]
            nbr_atr = node_attributes[nbr_id]
            nbr = Tuple([nbr_atr[l] for l in levels[1:level]])
            wght = modify_edge_weights!(modified_nbr_weights, graph,
                                        node_set[node], node_set[nbr],
                                        node, nbr)
            if wght != sub_graph.weights[ii, jj]
                sub_graph.weights[ii, jj] = wght
                sub_graph.weights[jj, ii] = wght
                # if wght == 0
                #     rem_edge!(sub_graph, ii, jj)
                # else
                #     add_edge!(sub_graph, ii, jj, wght)
                # end
            end
        end
        if level < num_levels-1
            modify_edge_weights!(modified_nbr_weights, node_set[node],
                                 graphs_by_level, vvmaps, graph, level+1,
                                 node)
        end
    end
end

""""""
function MultiLevelSubGraph(
    graph::MultiLevelGraph,
    node_set::Dict{Tuple{Vararg{String}},Any}
)::MultiLevelSubGraph
    num_levels = graph.num_levels

    graphs_by_level = Vector{Dict{Tuple{Vararg{String}},
        SimpleWeightedGraph}}(undef, num_levels)
    vmaps = Vector{Dict{Tuple{Vararg{String}}, Vector{Int}}}(undef, num_levels)
    modified_node_populations = Vector{Dict{Tuple{Vararg{String}}, Int}}(undef,
        num_levels)
    modified_nbr_weights = Dict{Set{Tuple{Vararg{String}}}, Any}()

    for ii = 1:num_levels
        graphs_by_level[ii] = Dict{Tuple{Vararg{String}}, SimpleWeightedGraph}()
        vmaps[ii] = Dict{Tuple{Vararg{String}}, Vector{Int}}()
        modified_node_populations[ii] = Dict{Tuple{Vararg{String}}, Int}()
    end

    build_subgraphs_by_level!(graphs_by_level, vmaps, graph, node_set)
    get_populations!(modified_node_populations, graphs_by_level, vmaps, graph)
    modify_edge_weights!(modified_nbr_weights, node_set, graphs_by_level,
                         vmaps, graph)

    return MultiLevelSubGraph(
        node_set,
        graphs_by_level,
        vmaps,
        modified_node_populations,
        modified_nbr_weights,
        graph,
    )
end
