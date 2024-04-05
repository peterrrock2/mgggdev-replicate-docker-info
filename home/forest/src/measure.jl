mutable struct Measure
    gamma::Float64
    weights::Vector{Float64}
    scores::Vector{Function}
    descriptions::Vector{String}
end

""""""
function Measure(gamma::Real)
    scores = Vector{Function}(undef, 0)
    weights = Vector{Float64}(undef, 0)
    descriptions = Vector{String}(undef, 0)
    return Measure(gamma, weights, scores, descriptions)
end

""""""
function log_measure(
    partition::MultiLevelPartition,
    measure::Measure;
    gamma::Union{Real,Nothing}=nothing,
    weights::Union{Vector,Nothing}=nothing
)
    log_p = 0
    if gamma == nothing
        gamma = measure.gamma
    end
    if gamma != 1
        log_linking_edges = get_log_linking_edges(partition)
        log_forests = get_log_spanning_forests(partition)
        log_p += (1-gamma)*(log_forests + log_linking_edges)
    end

    log_p -= get_log_energy(partition, measure, weights=weights)

    return log_p
end


""""""
function push_measure!(
    measure::Measure,
    score::Function,
    weight::Real;
    desc::String=""
)
    @assert length(measure.weights) == length(measure.scores)
    push!(measure.weights, weight)
    push!(measure.scores, score)
    if desc == ""
        desc = string(score)
    end
    push!(measure.descriptions, desc)
end

""""""
function get_delta_energy(
    partition::MultiLevelPartition, 
    measure::Measure,
    update::Tuple
)
    score = 0

    if length(update) == 1 && typeof(update[1])==MultiLevelPartition
        proposed_partition = update[1]
        score += get_log_energy(partition, measure)
        #println("score",score)
        score -= get_log_energy(proposed_partition, measure)
        #println("score",score)
        return exp(score)
    end

    changed_districts, node_sets_w_pops, edge = update
    for ii = 1:length(measure.weights)
        weight = measure.weights[ii]
        if weight == 0
            continue
        end
        energy = measure.scores[ii]
        score += weight*energy(partition, changed_districts)
    end

    old_dists = [partition.district_to_nodes[cd] for cd in changed_districts]
    old_pops = [partition.dist_populations[cd] for cd in changed_districts]
    for (ii, cd) in enumerate(changed_districts)
        partition.district_to_nodes[cd] = node_sets_w_pops[ii][1][()]
        partition.dist_populations[cd] = node_sets_w_pops[ii][2]
    end
    # partition.node_to_district = 
    # set_cross_district_edges!(partition, changed_districts)
    for ii = 1:length(measure.weights)
        weight = measure.weights[ii]
        if weight == 0
            continue
        end
        energy = measure.scores[ii]
        score -= weight*energy(partition, changed_districts)
    end
    for (ii, cd) in enumerate(changed_districts)
        partition.district_to_nodes[cd] = old_dists[ii]
        partition.dist_populations[cd] = old_pops[ii]
    end
    # partition.node_to_district = 
    # set_cross_district_edges!(partition, changed_districts)
    return exp(score)
end

""""""
function get_log_energy(
    partition::MultiLevelPartition, 
    measure::Measure;
    weights::Union{Vector, Nothing}=nothing
)
    score = 0
    if weights == nothing
        weights = measure.weights
    end
    for ii = 1:length(measure.weights)
        weight = weights[ii]
        if weight == 0
            continue
        end
        energy = measure.scores[ii]
        score += weight*energy(partition)
    end
    return score
end

""""""
function get_cut_edge_count(partition::MultiLevelPartition)
    return get_cut_edge_weights(partition, "connections")
end


""""""
function get_cut_edge_perimeter(partition::MultiLevelPartition)
    edge_perimeter_col = partition.graph.graphs_by_level[1].edge_perimeter_col
    return get_cut_edge_weights(partition, edge_perimeter_col)
end


""""""
function get_cut_edge_weights(
    partition::MultiLevelPartition,
    node::Tuple{Vararg{String}},
    district::Int,
    fine_nbr_ids_near_node::Dict{Int,Any},
    column::String,
    level::Int
)
    cut_edge_sum = 0
    graph = partition.graph
    for nbr_id in keys(fine_nbr_ids_near_node)
        nbr = graph.id_to_partitions[level][nbr_id]
        nbr_district = get_district(partition, nbr)
        if nbr_district == nothing
            cut_edge_sum += get_cut_edge_weights(partition, node, district,
                                                 fine_nbr_ids_near_node[nbr_id],
                                                 column, level+1)
        elseif nbr_district < district
            node_level = length(node)
            full_graph = graph.graphs_by_level[level]
            simple_graph = full_graph.simple_graph
            for nbr_nbr_id in neighbors(simple_graph, nbr_id)
                nbr_nbr = graph.id_to_partitions[level][nbr_nbr_id]
                if nbr_nbr[1:node_level]==node
                    edge = Set([nbr_id, nbr_nbr_id])
                    cut_edge_sum += full_graph.edge_attributes[edge][column]
                end
            end
        end
    end
    return cut_edge_sum
end

""""""
function get_cut_edge_weights(
    partition::MultiLevelPartition,
    column::String
)
    cut_edge_sum = 0
    graph = partition.graph

    for (node, district) in partition.node_to_district
        level = length(node)
        full_graph = graph.graphs_by_level[level]
        simple_graph = full_graph.simple_graph
        node_id = graph.partition_to_ids[level][node]
        for nbr_id in neighbors(simple_graph, node_id)
            nbr = graph.id_to_partitions[level][nbr_id]
            nbr_district = get_district(partition, nbr)
            if nbr_district == nothing
                fine_nbr_ids_near_node = graph.fine_neighbors[level][node_id]
                fine_nbr_ids_near_node = fine_nbr_ids_near_node[nbr_id]
                cut_edge_sum += get_cut_edge_weights(partition, node, district,
                                                     fine_nbr_ids_near_node,
                                                     column, level+1)
            elseif nbr_district < district
                edge = Set([node_id, nbr_id])
                cut_edge_sum += full_graph.edge_attributes[edge][column]
            end
        end
    end
    return cut_edge_sum
end


""""""
function get_log_spanning_forests(
    partition::MultiLevelPartition,
    districts::Vector{Int}=collect(1:partition.num_dists))
    return sum(get_log_spanning_trees(partition, districts))
end


""""""
function get_log_spanning_trees(
    partition::MultiLevelPartition,
    districts::Vector{Int}=collect(1:partition.num_dists)
)
    log_spanning_trees = Vector{Float64}(undef, 0)
    precompute_node_tree_counts!(partition)

    graph = partition.graph
    node_tree_counts = partition.extensions[node_trees::EXTENSIONS]

    for di in districts
        subgraph = MultiLevelSubGraph(graph, partition.district_to_nodes[di])
        lsptr = get_log_spanning_trees(subgraph, node_tree_counts)
        push!(log_spanning_trees, lsptr)
    end
    return log_spanning_trees
end


""""""
function get_log_spanning_trees(
    subgraph::MultiLevelSubGraph,
    node_tree_counts::Vector{Vector{Float64}},
    node_set::Union{Nothing,Dict{Tuple{Vararg{String}},Any}}=subgraph.node_set,
    key::Tuple=(),
    level::Int=1,
)
    log_spanning_trees = 0
    graph = subgraph.parent

    if intact_huh(node_set)
        node_id = graph.partition_to_ids[level-1][key]
        log_spanning_trees += node_tree_counts[level-1][node_id]
        node_set = get_intact_node_set(graph, key)
    else
        graph_at_level = subgraph.graphs_by_level[level][key]
        log_spanning_trees += log_nspanning(graph_at_level)
    end


    if level < subgraph.parent.num_levels
        for node in keys(node_set)
            log_spanning_trees += get_log_spanning_trees(subgraph,
                                                         node_tree_counts,
                                                         node_set[node], node,
                                                         level+1)
        end
    end
    return log_spanning_trees
end

""""""
function get_log_linking_edges(
    partition::MultiLevelPartition,
    districts::Vector{Int}=collect(1:partition.num_dists)
)
    cross_district_edges = partition.cross_district_edges
    district_to_nodes = partition.district_to_nodes

    log_linking_edges = 0

    num_levels = partition.graph.num_levels

    for di in 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            if di ∉ districts && dj ∉ districts
                continue
            else
                max_connections = maximum(cross_district_edges[di, dj, ii]
                                          for ii = 1:num_levels)
                adjacent_huh = max_connections!=0
                if !adjacent_huh
                    continue
                elseif !exists_hierarchical_tree(district_to_nodes[di],
                                                  district_to_nodes[dj],
                                                  num_levels)
                    continue
                end
                level = maximum(ii for ii = 1:num_levels
                                if cross_district_edges[di, dj, ii] > 0)
                choices = cross_district_edges[di, dj, level]
                log_linking_edges += log(choices)
            end
        end
    end
    return log_linking_edges
end
