""""""
function get_adjacent_tree_dists(
    partition::MultiLevelPartition,
    constraints::Dict,
)::Vector{Pair{Int,Int}}
    adj_dists = Vector{Pair{Int,Int}}(undef, 0)

    for di = 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            adjacent = false
            for lv = 1:partition.graph.num_levels
                if partition.cross_district_edges[di, dj, lv] > 0
                    adjacent = true
                end
            end
            if adjacent
                if !exists_hierarchical_tree(partition.district_to_nodes[di],
                                             partition.district_to_nodes[dj],
                                             partition.graph.num_levels)
                    continue
                end
                # println("searching on ", di, " ", dj, " ", total_adjacent_dists)
                push!(adj_dists, Pair(di, dj))
            end
        end
    end

    return adj_dists
end

""""""
function get_total_adjacent_dists(
    partition::MultiLevelPartition,
    constraints::Dict,
)::Int
    total_adjacent_dists = 0

    for di = 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            adjacent = false
            for lv = 1:partition.graph.num_levels
                if partition.cross_district_edges[di, dj, lv] > 0
                    adjacent = true
                end
            end
            if adjacent
                # println("searching on ", di, " ", dj, " ", total_adjacent_dists)
                total_adjacent_dists += 1
            end
        end
    end

    return total_adjacent_dists
end

""""""
function sample_adjacent_districts_randomly(
    partition::MultiLevelPartition,
    constraints::Dict,
    rng::AbstractRNG
)
    adj_dists = get_adjacent_tree_dists(partition, constraints)
    dist_pair = rand(rng, adj_dists)
    return [dist_pair.first, dist_pair.second], 1.0/length(adj_dists)
    # total_adjacent_dists = get_total_adjacent_dists(partition, constraints)

    # cum_sum = 0
    # choice = rand(rng)*total_adjacent_dists

    # for di = 1:partition.num_dists-1
    #     for dj = di+1:partition.num_dists
    #         adjacent = false
    #         for lv = 1:partition.graph.num_levels
    #             if partition.cross_district_edges[di, dj, lv] > 0
    #                 adjacent = true
    #             end
    #         end
    #         if adjacent
    #             # println("testing on ", di, " ", dj, " ", total_adjacent_dists, " ", choice)
    #             cum_sum += 1
    #             if cum_sum > choice
    #                 prob_of_dists = 1/total_adjacent_dists
    #                 return [di, dj], prob_of_dists
    #             end
    #         end
    #     end
    # end
end

""""""
function sample_subgraph2(
    partition::MultiLevelPartition,
    constraints::Dict,
    rng::AbstractRNG)
    proposed_dists = sample_adjacent_districts_randomly(partition, constraints, rng)

    changed_districts, prob_of_dists = proposed_dists
    node_sets = [partition.district_to_nodes[di] for di in changed_districts]

    node_set = merge_nodesets(partition.graph, node_sets)
    subgraph = MultiLevelSubGraph(partition.graph, node_set)
    return changed_districts, prob_of_dists, subgraph
end

""""""
function exists_hierarchical_tree(
    district1::Dict{Tuple{Vararg{String}}, Any},
    district2::Dict{Tuple{Vararg{String}}, Any},
    num_levels::Int,
    cur_level::Int=1
)
    num_of_shared = 0
    for node in keys(district1)
        if intact_huh(district1[node])
            continue
        end
        if haskey(district2, node)
            num_of_shared += 1
            if num_of_shared > 1
                return false
            end
            if cur_level != num_levels-1
                check_finer = exists_hierarchical_tree(district1[node],
                                                        district2[node],
                                                        num_levels,
                                                        cur_level+1)
                if !check_finer
                    return false
                end
            end
        end
    end
    return true
end

""""""
function check_prestep_constraints(
    subgraph::MultiLevelSubGraph,
    partition::MultiLevelPartition,
    constraints::Dict,
    changed_districts::Vector{Int}
)
    if haskey(constraints, ConstrainDiscontinuousTraversals)
        if !satisfies_constraint(constraints[ConstrainDiscontinuousTraversals],
                                 partition, subgraph)
            return false
        end
    end

    num_levels = partition.graph.num_levels
    district1 = partition.district_to_nodes[changed_districts[1]]
    district2 = partition.district_to_nodes[changed_districts[2]]
    if !exists_hierarchical_tree(district1, district2, num_levels)
        return false
    end

    return true
end

""""""
function get_prob_of_new_dists(
    partition::MultiLevelPartition,
    changed_districts::Vector{Int},
    constraints::Dict,
    node_sets_w_pops,
)
    old_dists = Dict(cd=>partition.district_to_nodes[cd]
                     for cd in changed_districts)
    partition.district_to_nodes[changed_districts[1]] =
        node_sets_w_pops[1][1][()]
    partition.district_to_nodes[changed_districts[2]] =
        node_sets_w_pops[2][1][()]
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)
    # total_adjacent_dists = get_total_adjacent_dists(partition, constraints)
    # prob_of_new_dists = 1.0/total_adjacent_dists
    adj_dists = get_adjacent_tree_dists(partition, constraints)
    prob_of_new_dists = 1.0/length(adj_dists)

    for (od, ons) in old_dists
        partition.district_to_nodes[od] = ons
    end
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)
    return prob_of_new_dists
end

""""""
function get_log_linking_edge_prob_tree_space(
    partition::MultiLevelPartition,
    districts::Vector{Int}=collect(1:partition.num_dists)
)
    cross_district_edges = partition.cross_district_edges
    district_to_nodes = partition.district_to_nodes

    log_linking_edge_prob = 0

    num_levels = partition.graph.num_levels

    for di in 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            di_nin = (di ∉ districts)
            dj_nin = (dj ∉ districts)
            if di_nin && dj_nin
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
                log_linking_edge_prob -= log(choices)
            end
        end
    end
    return log_linking_edge_prob
end

""""""
function get_log_linking_edge_ratio_tree_space(
    partition::MultiLevelPartition,
    measure::Measure,
    changed_districts::Vector{Int},
    node_sets_w_pops,
)
    if measure.gamma == 0
        return 0
    end
    num_levels = partition.graph.num_levels

    log_l_edge_prob_cur = get_log_linking_edge_prob_tree_space(partition,
                                                              changed_districts)

    old_dists = Dict(cd=>partition.district_to_nodes[cd]
                     for cd in changed_districts)
    for ii = 1:length(changed_districts)
        partition.district_to_nodes[changed_districts[ii]] =
            node_sets_w_pops[ii][1][()]
    end
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)
    log_l_edge_prob_proposed = get_log_linking_edge_prob_tree_space(partition,
                                                              changed_districts)

    for (od, ons) in old_dists
        partition.district_to_nodes[od] = ons
    end
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)

    return log_l_edge_prob_proposed - log_l_edge_prob_cur
end

""""""
function get_log_linking_edge_prob_partition_space(
    partition::MultiLevelPartition,
    districts::Vector{Int}=collect(1:partition.num_dists)
)
    cross_district_edges = partition.cross_district_edges
    district_to_nodes = partition.district_to_nodes

    log_linking_edge_prob = 0

    num_levels = partition.graph.num_levels

    for di in 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            if di ∉ districts && dj ∉ districts
                continue
            elseif di in districts && dj in districts
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
                log_linking_edge_prob += log(choices)
            end
        end
    end
    return log_linking_edge_prob
end

""""""
function get_log_linking_edge_ratio_partition_space(
    partition::MultiLevelPartition,
    measure::Measure,
    changed_districts::Vector{Int},
    node_sets_w_pops,
)
    if measure.gamma == 1
        return 0
    end
    num_levels = partition.graph.num_levels

    log_choices_cur = get_log_linking_edge_prob_partition_space(partition,
                                                              changed_districts)

    old_dists = Dict(cd=>partition.district_to_nodes[cd]
                     for cd in changed_districts)
    for ii = 1:length(changed_districts)
        partition.district_to_nodes[changed_districts[ii]] =
            node_sets_w_pops[ii][1][()]
    end
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)
    log_choices_proposed = get_log_linking_edge_prob_partition_space(partition,
                                                              changed_districts)

    for (od, ons) in old_dists
        partition.district_to_nodes[od] = ons
    end
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)

    return log_choices_proposed - log_choices_cur
end

""""""
function get_log_linking_edge_prob_adjacent(
    partition::MultiLevelPartition,
    measure::Measure,
    rng::AbstractRNG,
    districts::Vector{Int}=collect(1:partition.num_dists)
)
    cross_district_edges = partition.cross_district_edges
    district_to_nodes = partition.district_to_nodes

    log_linking_edge_prob = 0

    num_levels = partition.graph.num_levels

    for di in 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            if di ∉ districts && dj ∉ districts
                continue
            elseif di in districts && dj in districts
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
                # edge, prob? = sample_linking_edge(partition, di, dj, rng)
                # refined_edge
                # edge_weight = get_edge_weight(partition.graph, edge)
                # log_linking_edge_prob += log(prob)

                level = maximum(ii for ii = 1:num_levels
                                if cross_district_edges[di, dj, ii] > 0)
                choices = cross_district_edges[di, dj, level]
                log_linking_edge_prob -= log(choices)
            end
        end
    end
    return log_linking_edge_prob
end

""""""
function get_log_linking_edge_ratio_adjacent(
    partition::MultiLevelPartition,
    measure::Measure,
    changed_districts::Vector{Int},
    node_sets_w_pops,
    rng::AbstractRNG
)
    if measure.gamma == 0
        return 0
    end
    num_levels = partition.graph.num_levels

    log_prob_cur = get_log_linking_edge_prob_adjacent(partition, measure,
                                                      rng, changed_districts)

    old_dists = Dict(cd=>partition.district_to_nodes[cd]
                     for cd in changed_districts)
    for ii = 1:length(changed_districts)
        partition.district_to_nodes[changed_districts[ii]] =
            node_sets_w_pops[ii][1][()]
    end
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)
    log_prob_proposed = get_log_linking_edge_prob_adjacent(partition, measure,
                                                           rng,
                                                           changed_districts)

    for (od, ons) in old_dists
        partition.district_to_nodes[od] = ons
    end
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)

    return log_prob_proposed - log_prob_cur
end

""""""
function get_log_tree_count_ratio(
    district_graph₁::MultiLevelSubGraph,
    district_graph₂::MultiLevelSubGraph,
    subgraph::MultiLevelSubGraph,
    node_tree_counts::Vector{Vector{Float64}},
    edge
)
    log_tree_count_ratio = 0
    n1, n2 = edge
    for ii = 1:min(length(n1), length(n2))
        coarse_node = n1[1:ii-1]
        if coarse_node != n2[1:ii-1]
            break
        end
        sg1 = district_graph₁.graphs_by_level[ii][coarse_node]
        sg2 = district_graph₂.graphs_by_level[ii][coarse_node]
        log_count1 = log_nspanning(sg1)
        log_count2 = log_nspanning(sg2)
        log_tree_count_ratio -= log_count1+log_count2

        if ii > 1
            joint_node_set = get_node_set(subgraph.node_set, coarse_node)
            if joint_node_set == nothing
                graph = subgraph.parent
                coarse_node_id = graph.partition_to_ids[ii-1][coarse_node]
                log_tree_count_ratio += node_tree_counts[ii-1][coarse_node_id]
            else
                sgj = subgraph.graphs_by_level[ii][coarse_node]
                log_countj = log_nspanning(sgj)
                log_tree_count_ratio += log_countj
            end
        end
    end
    return log_tree_count_ratio
end


""""""
function get_log_tree_count_ratio(
    partition::MultiLevelPartition,
    measure::Measure,
    subgraph::MultiLevelSubGraph,
    changed_districts::Vector{Int},
    proposed_cut,
    old_edges
)
    if measure.gamma == 0
        return 0
    end

    log_tree_count_ratio = 0

    graph = partition.graph
    # TODO: add extensions
    extensions = partition.extensions
    # tracked_trees = extensions[tracked_trees::EXTENSIONS]
    node_tree_counts = extensions[node_trees::EXTENSIONS]
    # proposed_extension[tracked_trees::EXTENSIONS] =
    #                                 Vector{Dict{Tuple{Vararg{String}}}, Float64}

    node_sets_w_pops, edge, prob_edge = proposed_cut
    subgraph₁ = MultiLevelSubGraph(graph, node_sets_w_pops[1][1][()])
    subgraph₂ = MultiLevelSubGraph(graph, node_sets_w_pops[2][1][()])

    log_tree_count_ratio += get_log_tree_count_ratio(subgraph₁, subgraph₂,
                                                     subgraph, node_tree_counts,
                                                     edge)

    subgraph₁ = partition.subgraphs[changed_districts[1]]
    subgraph₂ = partition.subgraphs[changed_districts[2]]

    old_edge = collect(old_edges[1])
    log_tree_count_ratio -= get_log_tree_count_ratio(subgraph₁, subgraph₂,
                                                     subgraph, node_tree_counts,
                                                     old_edge)

    return log_tree_count_ratio
end

""""""
function set_node_sets_order( #currently unused
    partition::MultiLevelPartition,
    subgraph::MultiLevelSubGraph,
    changed_districts::Vector{Int},
    node_sets_w_pops
)
    subgraph_pop = get_total_graph_population(subgraph)

    cur_node_set₁ = partition.district_to_nodes[changed_districts[1]]
    prop_node_set₁ = node_sets_w_pops[1][1][()]

    cur_pop_1 = partition.dist_populations[changed_districts[1]]
    prop_pop_1 = node_sets_w_pops[1][2]

    pop_col = partition.graph.graphs_by_level[1].pop_col
    intersection11 = intersect_node_sets(cur_node_set₁, prop_node_set₁)
    common_pop_11 = sum_node_data(partition.graph, intersection11, pop_col)

    h11 = cur_pop_1 + prop_pop_1 - 2*common_pop_11
    h12 = subgraph_pop - h11
    if h12 < h11
        return h12, (node_sets_w_pops[2], node_sets_w_pops[1])
    else
        return h11, node_sets_w_pops
    end
end

""""""
function set_node_sets_order_general(
    partition::MultiLevelPartition,
    subgraph::MultiLevelSubGraph,
    changed_districts::Vector{Int},
    node_sets_w_pops::Array
)
    k = length(changed_districts)
    subgraph_pop = get_total_graph_population(subgraph)
    cur_node_sets = [partition.district_to_nodes[changed_districts[i]] for i = 1:k]
    prop_node_sets = [node_sets_w_pops[i][1][()] for i = 1:k]
    pop_col = partition.graph.graphs_by_level[1].pop_col
    pop_intersections = Array{Float64,2}(undef,k,k)
    for i = 1:k
        for j = 1:k
            intersection = intersect_node_sets(cur_node_sets[i], prop_node_sets[j])
            pop_intersections[i,j] = sum_node_data(partition.graph, intersection, pop_col)
        end
    end
    neg_pop_intersections = pop_intersections.*-1 #flip sign since hungarian() minimizes
    assignment = hungarian(neg_pop_intersections)
    return subgraph_pop+assignment[2], Tuple(node_sets_w_pops[assignment[1][i]] for i = 1:k)
end

""""""
function forest_recom2!(
    partition::MultiLevelPartition,
    measure::Measure,
    constraints::Dict,
    rng::AbstractRNG
)
    # println("running fr2")
    graph = partition.graph

    proposed_subregion = sample_subgraph2(partition, constraints, rng)
    changed_districts, prob_of_dists, subgraph = proposed_subregion
    # println("changed_districts ", changed_districts, " ", prob_of_dists)
    # @show changed_districts

    if !check_prestep_constraints(subgraph, partition, constraints,
                                 changed_districts)
        return 0, nothing
    end

    multiscale_cuttable_tree = MultiScaleCuttableTree(graph.num_levels)

    construct_cuttable_tree!(multiscale_cuttable_tree, subgraph, constraints,
                             rng=rng)

    proposed_cut = cut_edge(multiscale_cuttable_tree, subgraph, rng)
    node_sets_w_pops, edge, prob_edge = proposed_cut

    if edge == nothing
        return 0, nothing
    end
    update = (changed_districts, node_sets_w_pops, edge)
    if !satisfies_constraints(partition, constraints, update)
        return 0, nothing
    end
    # check hamming distances
    hamming, node_sets_w_pops = set_node_sets_order_general(partition, subgraph,
                                                            changed_districts,
                                                            node_sets_w_pops)
    if hamming == 0
        return 0, nothing
    end

    old_tree = sample_merged_tree2(changed_districts, subgraph, partition, 
                                   constraints, rng)
    old_multiscale_cuttable_tree, prob_old_edge = old_tree

    if prob_old_edge == 0
        println("couldn't draw old tree -- problem")
        @assert false
    end

    prob_of_new_dists = get_prob_of_new_dists(partition, changed_districts,
                                              constraints, node_sets_w_pops)

    partialFwdPrpProb = prob_of_dists*prob_edge
    partialBckPrpProb = prob_of_new_dists*prob_old_edge
    p = partialBckPrpProb/partialFwdPrpProb
    # println("partialBckPrpProb = ", prob_of_new_dists, " * ", prob_old_edge)
    # println("partialFwdPrpProb = ", prob_of_dists, " * ", prob_edge)

    if measure.gamma != 0
        log_linking_edge_ratio = get_log_linking_edge_ratio_tree_space(
                                                            partition, measure,
                                                            changed_districts,
                                                            node_sets_w_pops)
        old_edge = old_multiscale_cuttable_tree[2]
        log_tree_count_ratio = get_log_tree_count_ratio(partition, measure,
                                                        subgraph,
                                                        changed_districts,
                                                        proposed_cut, old_edge)
        p*=exp(measure.gamma*(log_linking_edge_ratio + log_tree_count_ratio))

        # absorbe into the log_linkin_edge_ratio for now; can do until move to 
        # weighted graphs
        # adjacent_edge_ratio = get_log_linking_edge_ratio_adjacent(
        #                                                     partition, measure,
        #                                                     changed_districts,
        #                                                     node_sets_w_pops,
        #                                                     rng)
        # p*=exp(measure.gamma*adjacent_edge_ratio)
    end

    # if measure.gamma != 1
    #     log_linking_edge_ratio = get_log_linking_edge_ratio_partition_space(
    #                                                         partition, measure,
    #                                                         changed_districts,
    #                                                         node_sets_w_pops)
    #     # println("p_cur = ", p, " and log_linking_edge_ratio = ", log_linking_edge_ratio)
    #     # println("modifying with linked edge measure on partitions ", p, " x ", exp((1-measure.gamma)*log_linking_edge_ratio))
    #     p*=exp((1-measure.gamma)*log_linking_edge_ratio)
    # end

    # println("p and info ", p, " ", log_linking_edge_ratio, " ", log_tree_count_ratio)
    # P(x'|x)pi(x) = P(x|x')pi(x')
    # A(x'|x)Q(x'|x)pi(x) = A(x|x')Q(x|x')pi(x')
    # A(x'|x) = min(1, Q(x|x')pi(x'))
    #                  Q(x'|x)pi(x))

    return p, (changed_districts, node_sets_w_pops, edge)
end

""""""
function specified_edge_probability(
    edge::Set{Tuple},
    subgraph::MultiLevelSubGraph
)
    n1, n2 = collect(edge)
    if length(n1) > length(n2)
        tmp = n1
        n1 = n2
        n2 = tmp
    end

    edge_weight = get_edge_weight(subgraph, n1, n2)
    possible_edges = get_edge_weight(subgraph, n1, n2[1:length(n1)])
    return edge_weight/possible_edges
end

""""""
function compute_proposed_tree_ratio(
    subgraph::MultiLevelSubGraph,
    old_multiscale_cuttable_tree::Tuple{Vector{MultiScaleCuttableTree},
                                         Vector{Set{T}}},
    new_multiscale_cuttable_tree::MultiScaleCuttableTree,
    rng::AbstractRNG
) where T<: Tuple{Vararg{String}}
    old_trees = old_multiscale_cuttable_tree[1]
    old_edges = old_multiscale_cuttable_tree[2]
    new_cuttable_trees = new_multiscale_cuttable_tree.cuttable_trees
    new_specified_edges = new_multiscale_cuttable_tree.specified_edges

    prob_forward = 1
    prob_backward = 1

    for edge in old_edges
        n1, n2 = collect(edge)
        edge_weight = get_edge_weight(subgraph, n1, n2)
        prob_backward *= edge_weight # /possible_edges;
                                     # cancels in merged tree
        # this is now the same as above, as we have specified the edge to
        # the finest level
        # refined_edge = refined_edge([n1, n2], subgraph, rng)
        # edge_weight = get_edge_weight(subgraph, refined_edge[1],
        #                               refined_edge[2])
        # prob_backward *= edge_weight # /possible_edges;
        #                              # cancels in merged tree

    end


    for level = 1:subgraph.parent.num_levels
        for key in keys(new_cuttable_trees[level])
            old_tree_ids = [jj for jj=1:length(old_trees)
                            if haskey(old_trees[jj].cuttable_trees[level], key)]
            if length(old_tree_ids) == 0
                continue
            end

            new_tree = new_cuttable_trees[level][key].tree
            prob_forward *= prod([new_tree.weights[src(e), dst(e)]
                                  for e in edges(new_tree)])
            for old_tree_id in old_tree_ids
                old_cuttable_trees = old_trees[old_tree_id].cuttable_trees
                old_tree = old_cuttable_trees[level][key].tree
                prob_backward *= prod([old_tree.weights[src(e), dst(e)]
                                       for e in edges(old_tree)])
            end

            # this is what cancels on the backward probability above in first loop
            # for edge in old_edges
            #     n1, n2 = [e[1:level] for e in collect(edge)]
            #     if n1[1:level-1] != key || n2[1:level-1] != key
            #         continue
            #     end
            #     prob_backward *= get_edge_weight(subgraph, n1, n2)
            # end

            for edge in values(new_specified_edges[level][key])
                # prob_forward *= specified_edge_probability(edge, subgraph)
                vec_edge = collect(edge)
                refined_edge = Set(refine_edge(vec_edge, subgraph, rng))
                prob_forward *= specified_edge_probability(refined_edge,
                                                           subgraph)
            end
            for old_tree_id in old_tree_ids
                old_tree = old_trees[old_tree_id]
                for edge in values(old_tree.specified_edges[level][key])
                    # prob_backward *= specified_edge_probability(edge, subgraph)
                    vec_edge = collect(edge)
                    refined_edge = Set(refine_edge(vec_edge, subgraph, rng))
                    prob_backward *= specified_edge_probability(refined_edge,
                                                                subgraph)
                end
            end
        end
    end
    return prob_backward/prob_forward
end

""""""
function sample_linking_edge(
    node_set₁::Dict{Tuple{Vararg{String}}, Any},
    node_set₂::Dict{Tuple{Vararg{String}}, Any},
    graph::MultiLevelGraph,
    weight::Real,
    edge_level::Int,
    rng::AbstractRNG,
    cur_level::Int=1
)
    if cur_level == edge_level
        simple_graph = graph.graphs_by_level[cur_level].simple_graph
        choice = rand(rng)*weight
        cum_sum = 0
        for node in keys(node_set₁)
            node_id = graph.partition_to_ids[cur_level][node]
            for nbr_id in neighbors(simple_graph, node_id)
                nbr = graph.id_to_partitions[cur_level][nbr_id]
                if haskey(node_set₂, nbr)
                    edge_weight = get_edge_weight(graph, node_set₁[node],
                                                  node_set₂[nbr], node, nbr)
                    cum_sum += edge_weight
                    if cum_sum > choice
                        return [node, nbr]
                    end
                end
            end
        end
    else # node_level < edge_level
        nodes = Vector{Tuple}(undef, 0)
        for node in keys(node_set₁)
            if !intact_huh(node_set₁[node]) && haskey(node_set₂, node)
                push!(nodes, node)
            end
        end
        if length(nodes) > 1
            return nothing
        end
        node = nodes[1]
        return sample_linking_edge(node_set₁[node], node_set₂[node], graph,
                                   weight, edge_level, rng, cur_level+1)
    end
end

""""""
function sample_linking_edge(
    partition::MultiLevelPartition,
    D₁::Integer, D₂::Integer, rng::AbstractRNG
)
    cross_district_edges = partition.cross_district_edges
    num_levels = partition.graph.num_levels
    edge_level = maximum([ii for ii = 1:num_levels
                          if cross_district_edges[D₁, D₂, ii] > 0])
    node_set₁ = partition.district_to_nodes[D₁]
    node_set₂ = partition.district_to_nodes[D₂]

    weight = partition.cross_district_edges[D₁, D₂, edge_level]
    edge = sample_linking_edge(node_set₁, node_set₂, partition.graph,
                               weight, edge_level, rng)
    return edge
end

""""""
function refine_linking_edge(
    edge::Vector{T},
    subgraph::MultiLevelSubGraph,
    subgraph₁::MultiLevelSubGraph,
    subgraph₂::MultiLevelSubGraph,
    rng::AbstractRNG
) where T <: Tuple
    @assert length(edge[1]) == length(edge[2])

    level = length(edge[1])
    graph = subgraph₁.parent
    num_levels = graph.num_levels

    if level == num_levels
        return edge
    end

    simple_graph = graph.graphs_by_level[level+1].simple_graph
    ns1 = get_node_set(subgraph₁, edge[1])
    ns2 = get_node_set(subgraph₂, edge[2])

    # suppose dict of nodes
    edge1_id = graph.partition_to_ids[level][edge[1]]
    edge2_id = graph.partition_to_ids[level][edge[2]]
    cum_weight = 0
    edge_weights = Dict{Vector{Tuple}, Real}()

    for fine_edge2_id in keys(graph.fine_neighbors[level][edge1_id][edge2_id])
        fine_edge2 = graph.id_to_partitions[level+1][fine_edge2_id]
        if !intact_huh(ns2) && fine_edge2 ∉ keys(ns2)
            continue
        end
        for fine_edge1_id in neighbors(simple_graph, fine_edge2_id)
            fine_edge1 = graph.id_to_partitions[level+1][fine_edge1_id]
            if fine_edge1[1:level] != edge[1]
                continue
            end
            if !intact_huh(ns1) && fine_edge1 ∉ keys(ns1)
                continue
            end
            weight = get_edge_weight(subgraph, fine_edge1, fine_edge2)
            cum_weight += weight
            edge_weights[[fine_edge1, fine_edge2]] = weight
        end
    end

    choice = rand(rng)*cum_weight
    cum_sum = 0
    for (edge_key, weight) in edge_weights
        cum_sum += weight
        if cum_sum > choice
            refined_edge = edge_key
            return refine_linking_edge(refined_edge, subgraph, subgraph₁, 
                                   subgraph₂, rng)
        end
    end
end

""""""
function refine_edge(
    edge::Vector{T},
    subgraph::MultiLevelSubGraph,
    rng::AbstractRNG,
    depth::Int=-1
) where T <: Tuple
    if length(edge[1]) > length(edge[2])
        edge = [edge[2], edge[1]]
    end
    level1 = length(edge[1])
    level2 = length(edge[2])
    graph = subgraph.parent
    num_levels = graph.num_levels

    if depth == -1
        depth = num_levels
    end
    if length(edge[1]) == length(edge[2]) == depth
        return edge
    end

    # simple_graph = graph.graphs_by_level[level+1].simple_graph
    # ns1 = get_node_set(subgraph, edge[1])
    # ns2 = get_node_set(subgraph, edge[2])

    # suppose dict of nodes
    edge1_id = graph.partition_to_ids[level1][edge[1]]
    edge2_id = graph.partition_to_ids[level2][edge[2]]
    cum_weight = 0
    edge_weights = Dict{Vector{Tuple}, Real}()

    if level1 == level2
        fine_neighbors = graph.fine_neighbors[level1][edge1_id][edge2_id]
        for fine_edge2_id in keys(fine_neighbors)
            fine_edge2 = graph.id_to_partitions[level2+1][fine_edge2_id]
            if !intersects_huh(subgraph, fine_edge2)
                continue
            end
            weight = get_edge_weight(subgraph, edge[1], fine_edge2)
            cum_weight += weight
            edge_weights[[edge[1], fine_edge2]] = weight
        end
    else # level1 < level2
        simple_graph = graph.graphs_by_level[level2].simple_graph
        for nbr_id in neighbors(simple_graph, edge2_id)
            nbr = graph.id_to_partitions[level2][nbr_id]
            if !intersects_huh(subgraph, nbr) || nbr[1:level1] != edge[1]
                continue
            end
            weight = get_edge_weight(subgraph, nbr, edge[2])
            cum_weight += weight
            edge_weights[[nbr, edge[2]]] = weight
        end
    end
    choice = rand(rng)*cum_weight
    cum_sum = 0
    for (edge_key, weight) in edge_weights
        cum_sum += weight
        if cum_sum > choice
            refined_edge = edge_key
            return refine_edge(refined_edge, subgraph, rng, depth)
        end
    end

end

""""""
function sample_merged_tree2(
    changed_districts::Vector{Int},
    subgraph₁₂::MultiLevelSubGraph,
    partition::MultiLevelPartition,
    constraints::Dict,
    rng::AbstractRNG
)
    # prob_selecting_old_edge_to_remove
    D₁, D₂ = changed_districts
    edge = sample_linking_edge(partition, D₁, D₂, rng)
    if edge == nothing
        println("returng nothing!")
        return nothing, 0
    end

    graph = partition.graph

    subgraph₁ = partition.subgraphs[D₁]
    subgraph₂ = partition.subgraphs[D₂]

    specified_edge = refine_linking_edge(edge, subgraph₁₂, subgraph₁, subgraph₂, 
                                         rng)

    multiscale_cuttable_tree₁ = MultiScaleCuttableTree(graph.num_levels)
    multiscale_cuttable_tree₂ = MultiScaleCuttableTree(graph.num_levels)
    pop₁ = partition.dist_populations[D₁]
    pop₂ = partition.dist_populations[D₂]
    exterior_decorators₁ = Dict{Tuple,Dict{Tuple,Int}}(specified_edge[1]=>
                                                  Dict(specified_edge[2]=>
                                                       pop₂))
    exterior_decorators₂ = Dict{Tuple,Dict{Tuple,Int}}(specified_edge[2]=>
                                                  Dict(specified_edge[1]=>
                                                       pop₁))
    construct_cuttable_tree!(multiscale_cuttable_tree₁, subgraph₁, constraints,
                             ext_decorators=exterior_decorators₁, rng=rng)
    # println("construct_cuttable_tree on D2")
    construct_cuttable_tree!(multiscale_cuttable_tree₂, subgraph₂, constraints,
                             ext_decorators=exterior_decorators₂, rng=rng)
    cuttable_edges = Set{Tuple}([Tuple(edge)])
    get_all_cuttable_edges!(cuttable_edges, multiscale_cuttable_tree₁)
    get_all_cuttable_edges!(cuttable_edges, multiscale_cuttable_tree₂)
    edge_prob = get_probability_of_edge(cuttable_edges) #, edge)

    old_tree = ([multiscale_cuttable_tree₁, multiscale_cuttable_tree₂],
                [Set(specified_edge)])
    return old_tree, edge_prob
end


""""""
function build_forest_recom2(
    constraints::Dict
)
    f(p, m, r) = forest_recom2!(p, m, constraints, r)
    return f
end
