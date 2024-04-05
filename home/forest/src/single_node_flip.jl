""""""
function single_node_flip!(
    partition::MultiLevelPartition,
    measure::Measure,
    constraints::Dict,
    rng::AbstractRNG,
    idx_dists::Vector{Int64}=collect(1:partition.num_dists),
    level::Int64=partition.graph.num_levels
)
    conflicted_edge_data = get_conflicted_edge(partition, rng, idx_dists, level)
    chosen_conflicted_edge, distpair, prob_edge = conflicted_edge_data
    coin = rand(rng,1:2)
    node_to_flip = chosen_conflicted_edge[coin]
    from_dist_idx = distpair[coin]
    to_dist_idx = distpair[mod(coin, 2)+1]
    if coin == 2
        distpair = [distpair[2], distpair[1]]
    end
##############
    # (node_to_flip, to_dist_idx, from_dist_idx) = (("FORSYTH", "042"), 1, 5)
    # distpair = [1,5]
    # @show get_district(partition, ("FORSYTH", "042"))
    # @show node_to_flip, to_dist_idx, from_dist_idx
##############

    to_dist_connections = get_neighbors_in_dist(partition, level, node_to_flip,
                                                to_dist_idx)
    prob_move_forward = prob_edge*to_dist_connections

    from_dist = deepcopy(partition.district_to_nodes[from_dist_idx])
    to_dist = deepcopy(partition.district_to_nodes[to_dist_idx])
    update_dists!(from_dist, to_dist, node_to_flip, partition.graph)
    correct_intact_nodes(partition.graph, to_dist)

##########
    level = length(node_to_flip)
    node_id = partition.graph.partition_to_ids[level][node_to_flip]
    graph = partition.graph.graphs_by_level[level]
    node_pop = graph.node_attributes[node_id][graph.pop_col]

    node_set1 = Dict{Tuple{Vararg{String}}, Any}(()=>from_dist)
    node_set2 = Dict{Tuple{Vararg{String}}, Any}(()=>to_dist)
    pop1 = partition.dist_populations[from_dist_idx]-node_pop
    pop2 = partition.dist_populations[to_dist_idx]+node_pop
    node_sets_w_pops = ((node_set1, pop1), (node_set2, pop2))
    update = (distpair, node_sets_w_pops, chosen_conflicted_edge)
##########
    # @show update
    if !satisfies_constraints(partition, constraints, update)
        # println("not sf const")
        return 0, nothing
    end

    if !check_connectivity(partition.graph, from_dist)
        # println("not concect from")
        return 0, nothing
    end
    if !check_connectivity(partition.graph, to_dist)
        # println("not concect to")
        return 0, nothing
    end

    from_dist_connections = get_neighbors_in_dist(partition, level, 
                                                  node_to_flip, from_dist_idx)
    cur_edges = (1.0/prob_edge)
    prob_edge_back = cur_edges - to_dist_connections + from_dist_connections
    prob_edge_back = 1.0/prob_edge_back
    prob_move_backward =  prob_edge_back * from_dist_connections
    p = prob_move_backward/prob_move_forward

    if measure.gamma != 1
        log_linking_edge_ratio = get_log_linking_edge_ratio_tree_space(
                                                            partition, measure,
                                                            distpair,
                                                            node_sets_w_pops)
        old_edge = [chosen_conflicted_edge]
        graph = partition.graph
        simple_graph = graph.graphs_by_level[level].simple_graph
        node_to_flip_id = graph.partition_to_ids[level][node_to_flip]
        new_nbr = ()
        cur_agreement = -1
        for nbr_id in neighbors(simple_graph, node_to_flip_id)
            nbr = partition.graph.id_to_partitions[level][nbr_id]
            if !intersects_huh(from_dist, nbr)
                continue
            end
            agrmnt = sum(nbr[ii]==node_to_flip[ii] for ii = 1:level)
            if agrmnt > cur_agreement
                cur_agreement = agrmnt
                new_nbr = nbr
            end
        end
        node_sets = [to_dist, from_dist]
        node_set = merge_nodesets(partition.graph, node_sets)
        subgraph = MultiLevelSubGraph(partition.graph, node_set)
        proposed_cut = node_sets_w_pops, (node_to_flip, new_nbr), 0
        # @show (node_to_flip, new_nbr)
        log_tree_count_ratio = get_log_tree_count_ratio(partition, measure,
                                                        subgraph,
                                                        distpair,
                                                        proposed_cut, old_edge)
        p*=exp((1-measure.gamma)*(log_linking_edge_ratio + log_tree_count_ratio))
        
        adjacent_edge_ratio = get_log_linking_edge_ratio_adjacent(
                                                            partition, measure,
                                                            distpair, 
                                                            node_sets_w_pops,
                                                            rng)
        p*=exp((1-measure.gamma)*adjacent_edge_ratio)
    end
    return p, update
end

function get_dist_idx(
    partition::MultiLevelPartition,
    chosen_node::Tuple{Vararg{String}}
)
    if haskey(partition.node_to_district, chosen_node)
        idx = partition.node_to_district[chosen_node]
    else
        idx = get_dist_idx(partition, chosen_node[1:length(chosen_node)-1])
    end
    return idx
end

function sample_adj_dists(
    partition::MultiLevelPartition, 
    rng::AbstractRNG
)
    cross_district_edges = partition.cross_district_edges
    choice = rand(rng)*sum(cross_district_edges[:,])/2
    cum_sum = 0
    for di = 1:partition.num_dists
        for dj = di+1:partition.num_dists
            for l = 1:partition.graph.num_levels
                cum_sum += cross_district_edges[di, dj, l]
                if cum_sum > choice
                    return (di, dj)
                end
            end
        end
    end
end

function get_conflicted_edge(
    partition::MultiLevelPartition, 
    rng::AbstractRNG,
    idx_dists::Vector{Int64}=collect(1:partition.num_dists), 
    level::Int64=partition.graph.num_levels
)
    (d1, d2) = sample_adj_dists(partition, rng)
    dist_pair_weight = sum(partition.cross_district_edges[d1, d2, :])
    choice = rand(rng)*dist_pair_weight
    district_to_nodes = partition.district_to_nodes
    edge = get_conflicted_edge_dist(partition, district_to_nodes[d1], 
                                    district_to_nodes[d2], choice)[1]
    # @show edge
    @assert edge != nothing
    edge_prob = 2.0/(sum(partition.cross_district_edges))
    return edge, [d1, d2], edge_prob
end       

function get_conflicted_edge_dist(
    partition::MultiLevelPartition,
    node_set1::Dict{Tuple{Vararg{String}},Any},
    node_set2::Dict{Tuple{Vararg{String}},Any},
    choice::Real,
    level::Int=1,
    cum_sum::Real=0,
)
    graph = partition.graph
    simple_graph = graph.graphs_by_level[level].simple_graph
    # conflicted_edges_dist = Vector{Dict{Union{Tuple{Vararg{String}},Int64},Union{Tuple{Vararg{String}}, Int64}}}()
    # conflicted_edges_dist = Vector{Dict{Any, Any}}()

    for (node, node_set) in node_set1
        # @show "before intersect", node, cum_sum, choice
        if haskey(node_set2, node)
            edge, cum_sum = get_conflicted_edge_dist(partition, node_set1[node],
                                                     node_set2[node], 
                                                     choice, level+1, 
                                                     cum_sum)
            if edge != nothing
                return edge, cum_sum
            end
        end
        # @show "afrter intersect", node, cum_sum, choice
        node_id = graph.partition_to_ids[level][node]
        #exterior -- across coarse nodes
        for nbr_id in neighbors(simple_graph, node_id)
            nbr = graph.id_to_partitions[level][nbr_id]
            if !intersects_huh(node_set2, nbr)
                continue
            end
            edge_weight = get_edge_weight(graph, node_set1[node], 
                                          node_set2[nbr], node, nbr)
            if cum_sum + edge_weight < choice
                cum_sum += edge_weight
                continue
            end
            if level == graph.num_levels
                return (node, nbr), cum_sum
            end 
            # @show node, nbr
            # @show graph.fine_neighbors[nbr_id]
            fine_nodes_near_nbr_ids = graph.fine_neighbors[level][nbr_id][node_id]
            fine_nodes_near_nbr_ids = keys(fine_nodes_near_nbr_ids)
            bndry_nodeset1 = Dict{Tuple{Vararg{String}},Any}()
            for fine_nodes_near_nbr_id in fine_nodes_near_nbr_ids
                fine_nodes_near_nbr = graph.id_to_partitions[level+1][fine_nodes_near_nbr_id]
                if intact_huh(node_set1[node])
                    bndry_nodeset1[fine_nodes_near_nbr] = nothing
                elseif haskey(node_set1[node], fine_nodes_near_nbr)
                    bndry_nodeset1[fine_nodes_near_nbr] = node_set1[node][fine_nodes_near_nbr]
                end
            end
            fine_nbrs_near_node_ids = graph.fine_neighbors[level][node_id][nbr_id]
            fine_nbrs_near_node_ids = keys(fine_nbrs_near_node_ids)
            bndry_nodeset2 = Dict{Tuple{Vararg{String}},Any}()
            for fine_nbrs_near_node_id in fine_nbrs_near_node_ids
                fine_nbrs_near_node = graph.id_to_partitions[level+1][fine_nbrs_near_node_id]
                if intact_huh(node_set2[nbr])
                    bndry_nodeset2[fine_nbrs_near_node] = nothing
                elseif haskey(node_set2[nbr], fine_nbrs_near_node)
                    bndry_nodeset2[fine_nbrs_near_node] = node_set2[nbr][fine_nbrs_near_node]
                end
            end
            return get_conflicted_edge_dist(partition, bndry_nodeset1, 
                                            bndry_nodeset2, choice, level+1,
                                            cum_sum)
        end

        # if !intact_huh(node_set1[node]) && haskey(node_set2, node)
        #     fine_simple_graph = graph.graphs_by_level[level+1].simple_graph
        #     cum_edge_weight = 0
        #     bndry_nodeset1 = Dict{Tuple{Vararg{String}},Any}()
        #     bndry_nodeset2 = Dict{Tuple{Vararg{String}},Any}()
        #     for fine_node in keys(node_set1[node])
        #         fine_node_id = graph.partition_to_ids[level+1][fine_node]
        #         for fine_nbr_id in neighbors(fine_simple_graph, fine_node_id)
        #             fine_nbr = graph.id_to_partitions[level+1][fine_nbr_id]
        #             if fine_nbr ∉ keys(node_set2[node])
        #                 continue
        #             end
        #             edge_weight = get_edge_weight(graph, 
        #                                           node_set1[node][fine_node],
        #                                           node_set2[node][fine_nbr],
        #                                           fine_node, fine_nbr)
        #             if cum_sum + edge_weight < choice
        #                 cum_sum += edge_weight
        #                 continue
        #             end
        #             if level+1 == graph.num_levels
        #                 return (fine_node, fine_nbr)
        #             end
        #             @show node, fine_node, fine_nbr, node_set1[node][fine_node]==nothing, node_set2[node][fine_nbr]==nothing
        #             return get_conflicted_edge_dist(partition, 
        #                                             node_set1[node][fine_node], 
        #                                             node_set2[node][fine_nbr],
        #                                             edge_weight, rng, level+=1)
        #         end
        #     end
        # end
    end
    return nothing, cum_sum
end

function get_neighbors_in_dist(
    partition::MultiLevelPartition, 
    level::Int64, 
    node_to_flip::Tuple{Vararg{String}}, 
    dist_idx::Int64
)
    same_node_color = 0
    graph = partition.graph
    simple_graph = graph.graphs_by_level[level].simple_graph
    node_id = graph.partition_to_ids[level][node_to_flip]
    for nbr_id in neighbors(simple_graph, node_id)
        nbr = graph.id_to_partitions[level][nbr_id]##same level neighbor.
        if get_district(partition, nbr) == dist_idx
            same_node_color += 1
        end
    end
    return same_node_color
end



function update_dists!(
    from_dist::Dict{Tuple{Vararg{String}}, Any},
    to_dist::Dict{Tuple{Vararg{String}}, Any},
    node_to_flip::Tuple{Vararg{String}},
    graph::MultiLevelGraph
)
    add_node_to_district!(to_dist, node_to_flip)
    remove_node_from_district!(from_dist, graph, node_to_flip)
end


function check_connectivity(
    graph::MultiLevelGraph,
    node_set::Dict{Tuple{Vararg{String}}, Any}
)
    subgraph = MultiLevelSubGraph(graph, node_set)
    for level = 1:graph.num_levels
        graphs_by_level = subgraph.graphs_by_level[level]
        for (key, graph) in graphs_by_level
            if !is_connected_bf(graph)
                return false
            end
        end
    end
    return true
end

""""""
function build_single_node_flip(
    constraints::Dict
)
    f(p, m, r) = single_node_flip!(p, m, constraints, r)
    return f
end
