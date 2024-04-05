""""""
function sample_three_adjacent_districts_randomly(
    partition::MultiLevelPartition,
    constraints::Dict,
    rng::AbstractRNG
)
    first_two_dists, pair_prob = sample_adjacent_districts_randomly(partition, constraints, rng)
    D1 = first_two_dists[1]
    D2 = first_two_dists[2]
    D1_nbrs = Set([])
    D2_nbrs = Set([])
    for di = 1:partition.num_dists
        if di in first_two_dists
            continue
        end
        for lv = 1:partition.graph.num_levels
            if partition.cross_district_edges[D1, di, lv] > 0
                push!(D1_nbrs, di)
            end
            if partition.cross_district_edges[D2, di, lv] > 0
                push!(D2_nbrs, di)
            end
        end
    end
    adjacent_dists = union(D1_nbrs,D2_nbrs)
    D3 = rand(rng,adjacent_dists)
    D3_nbrs = Set([])
    for di = 1:partition.num_dists
        for lv = 1:partition.graph.num_levels
            if partition.cross_district_edges[D3, di, lv] > 0
                push!(D3_nbrs, di)
            end
        end
    end
    trio_prob = 1/length(union(D1_nbrs,D2_nbrs))
    if D1 in D3_nbrs
        trio_prob += 1/(length(union(D1_nbrs,D3_nbrs))-1)
    end
    if D2 in D3_nbrs
        trio_prob += 1/(length(union(D2_nbrs,D3_nbrs))-1)
    end
    trio_prob *= pair_prob
    return [D1,D2,D3],trio_prob
end

""""""
function sample_subgraph3(
    partition::MultiLevelPartition,
    constraints::Dict,
    rng::AbstractRNG
)
    proposed_dists = sample_three_adjacent_districts_randomly(partition, constraints, rng)
    changed_districts, prob_of_dists = proposed_dists
############
    # changed_districts = [14, 7, 6]
    # println("manually setting changed districts ", changed_districts)
############
    node_sets = [partition.district_to_nodes[di] for di in changed_districts]
    node_set = merge_nodesets(partition.graph, node_sets)
    subgraph = MultiLevelSubGraph(partition.graph, node_set)
    return changed_districts, prob_of_dists, subgraph
end

""""""
function get_prob_of_new_dists3(
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
    total_adjacent_dists = get_total_adjacent_dists(partition, constraints)
    prob_of_pair_dists = 1.0/total_adjacent_dists

    D1 = changed_districts[1]
    D2 = changed_districts[2]
    D3 = changed_districts[3]
    D1_nbrs = Set([])
    D2_nbrs = Set([])
    D3_nbrs = Set([])
    for di = 1:partition.num_dists
        for lv = 1:partition.graph.num_levels
            if partition.cross_district_edges[D1, di, lv] > 0
                push!(D1_nbrs, di)
            end
            if partition.cross_district_edges[D2, di, lv] > 0
                push!(D2_nbrs, di)
            end
            if partition.cross_district_edges[D3, di, lv] > 0
                push!(D3_nbrs, di)
            end
        end
    end
    trio_prob = 0
    if D1 in D2_nbrs
        trio_prob += 1/(length(union(D1_nbrs,D2_nbrs))-2)
    end
    if D1 in D3_nbrs
        trio_prob += 1/(length(union(D1_nbrs,D3_nbrs))-2)
    end
    if D2 in D3_nbrs
        trio_prob += 1/(length(union(D2_nbrs,D3_nbrs))-2)
    end
    trio_prob *= prob_of_pair_dists

    for (od, ons) in old_dists
        partition.district_to_nodes[od] = ons
    end
    partition.node_to_district = construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)
    return trio_prob
end

""""""
function get_log_tree_count_ratio3(
    district_graphs::Vector{MultiLevelSubGraph},
    subgraph::MultiLevelSubGraph,
    node_tree_counts::Vector{Vector{Float64}},
    edges
)
    log_tree_count_ratio = 0
    for ii = 1:subgraph.parent.num_levels
        coarse_edges = [[e[1][1:ii-1], e[2][1:ii-1]] for e in edges]
        coarse_nodes = Set([e[1] for e in coarse_edges if e[1]==e[2]])
        for coarse_node in coarse_nodes
            sgs = [dg.graphs_by_level[ii][coarse_node] for dg in district_graphs
                   if intersects_huh(dg, coarse_node)]
            tree_sum = sum(log_nspanning(sg) for sg in sgs)
            log_tree_count_ratio -= tree_sum

            if ii > 1
                if intact_huh(get_node_set(subgraph.node_set, coarse_node))
                    graph = subgraph.parent
                    c_node_id = graph.partition_to_ids[ii-1][coarse_node]
                    log_tree_count_ratio += node_tree_counts[ii-1][c_node_id]
                else
                    sgj = subgraph.graphs_by_level[ii][coarse_node]
                    log_countj = log_nspanning(sgj)
                    log_tree_count_ratio += log_countj
                end
            end
        end
    end
    return log_tree_count_ratio
end

""""""
function get_log_tree_count_ratio3(
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

    node_sets_w_pops, edges, prob_edge = proposed_cut
    edges = collect(values(edges))
    distgraphs = [MultiLevelSubGraph(graph, node_sets_w_pops[ii][1][()])
                 for ii=1:length(node_sets_w_pops)]

    log_tree_count_ratio += get_log_tree_count_ratio3(distgraphs, subgraph,
                                                     node_tree_counts, edges)

    node_sets = [partition.district_to_nodes[cd] for cd in changed_districts]
    distgraphs = [MultiLevelSubGraph(graph, node_set) for node_set in node_sets]

    log_tree_count_ratio -= get_log_tree_count_ratio3(distgraphs, subgraph,
                                                     node_tree_counts,
                                                     old_edges)

    return log_tree_count_ratio
end

""""""
function forest_recom3!(
    partition::MultiLevelPartition,
    measure::Measure,
    constraints::Dict,
    rng::AbstractRNG
)
    # println("running fr3")

    graph = partition.graph

    changed_districts, prob_of_dists, subgraph =
                                   sample_subgraph3(partition, constraints, rng)

    # println("changed_districts ", changed_districts)
    ### TODO-- test this
    if !check_prestep_constraints3(subgraph, partition, constraints,
                                   changed_districts)
        # println("failed check_prestep_constraints3")
        return 0, nothing
    end
    # println("passed check_prestep_constraints3")

    multiscale_cuttable_tree = MultiScaleCuttableTree(graph.num_levels)
    construct_cuttable_tree!(multiscale_cuttable_tree, subgraph, constraints,
                             rng=rng, balance=(1,2))

    proposed_cut = cut_edge2(multiscale_cuttable_tree, subgraph, constraints,
                             rng)
    node_sets_w_pops, edges, prob_edge = proposed_cut

    if edges == nothing
        if haskey(partition.extensions, rejection_counter::EXTENSIONS)
            partition.extensions[rejection_counter::EXTENSIONS]["rejection_noedge"]+=1
        end
        return 0, nothing
    end
    update = (changed_districts, node_sets_w_pops, edges)
    if !satisfies_constraints(partition, constraints, update)
        partition.extensions[rejection_counter::EXTENSIONS]["rejection_fr3_constraint"]+=1
        return 0, nothing
    end
    # println("found edges and statisfies constraints")

    # check hamming distances
    hamming, node_sets_w_pops = set_node_sets_order_general(partition, subgraph,
                                                            changed_districts,
                                                            node_sets_w_pops)
    if hamming == 0
        if haskey(partition.extensions, rejection_counter::EXTENSIONS)
            partition.extensions[rejection_counter::EXTENSIONS]["rejection_proposingSameDists"]+=1
        end
        # println("rejection on Hamming")
        return 0, nothing
    end

    old_tree = sample_merged_tree3(changed_districts, partition, constraints,
                                  rng)
    old_multiscale_cuttable_tree, prob_old_edge = old_tree
    # println("after sampled merge tree ", prob_old_edge)

    if prob_old_edge == 0
        return 0, nothing
    end

    prob_of_new_dists = get_prob_of_new_dists3(partition, changed_districts,
                                               constraints, node_sets_w_pops)

    partialFwdPrpProb = prob_of_dists*prob_edge
    partialBckPrpProb = prob_of_new_dists*prob_old_edge
    p = partialBckPrpProb/partialFwdPrpProb

    if measure.gamma != 0
        # won't work with weighted edges
        log_linking_edge_ratio = get_log_linking_edge_ratio_tree_space(
                                                            partition, measure,
                                                            changed_districts,
                                                            node_sets_w_pops)
        old_edges = values(old_multiscale_cuttable_tree[2])
        log_tree_count_ratio = get_log_tree_count_ratio3(partition, measure,
                                                        subgraph,
                                                        changed_districts,
                                                        proposed_cut, old_edges)
        p*=exp(measure.gamma*(log_linking_edge_ratio + log_tree_count_ratio))

        adjacent_edge_ratio = get_log_linking_edge_ratio_adjacent(
                                                            partition, measure,
                                                            changed_districts,
                                                            node_sets_w_pops,
                                                            rng)
        # @show adjacent_edge_ratio
        p*=exp(measure.gamma*adjacent_edge_ratio)
    end

    # if measure.gamma != 1
    #     log_linking_edge_ratio = get_log_linking_edge_ratio_partition_space(
    #                                                         partition, measure,
    #                                                         changed_districts,
    #                                                         node_sets_w_pops)
    #     p*=exp((1-measure.gamma)*log_linking_edge_ratio)
    # end

    # println("p and info ", p, " ", log_linking_edge_ratio, " ", log_tree_count_ratio)
    # P(x'|x)pi(x) = P(x|x')pi(x')
    # A(x'|x)Q(x'|x)pi(x) = A(x|x')Q(x|x')pi(x')
    # A(x'|x) = min(1, Q(x|x')pi(x'))
    #                  Q(x'|x)pi(x))

    return p, (changed_districts, node_sets_w_pops, edges)
end

""""""
function check_prestep_constraints3(
    subgraph::MultiLevelSubGraph,
    partition::MultiLevelPartition,
    constraints::Dict,
    changed_districts::Vector{Int}
)
    if ConstrainDiscontinuousTraversals in keys(constraints)
        if !satisfies_constraint(constraints[ConstrainDiscontinuousTraversals],
                                 partition, subgraph)
            return false
        end
    end

    num_levels = partition.graph.num_levels
    district1 = partition.district_to_nodes[changed_districts[1]]
    district2 = partition.district_to_nodes[changed_districts[2]]
    district3 = partition.district_to_nodes[changed_districts[3]]
    if !exists_hierarchical_tree3(district1, district2, district3, num_levels)
        # println("failing exists_hierarchical_tree3")
        return false
    end

    return true
end

""""""
function exists_hierarchical_tree3(
    district1::Dict{Tuple{Vararg{String}}, Any},
    district2::Dict{Tuple{Vararg{String}}, Any},
    district3::Dict{Tuple{Vararg{String}}, Any},
    num_levels::Int,
    cur_level::Int=1
)

    # check there is a tree between any two districts
    shared12 = hierarchical_tree_shared_nodes(district1, district2, num_levels)
    if shared12 == nothing
        return false
    end
    shared13 = hierarchical_tree_shared_nodes(district1, district3, num_levels)
    if shared13 == nothing
        return false
    end
    shared23 = hierarchical_tree_shared_nodes(district2, district3, num_levels)
    if shared23 == nothing
        return false
    end

    # println("looking at shared in exists_hierarchical_tree3 ", shared12, shared23, shared13)
    if minimum([length(shared12), length(shared13), length(shared23)]) == 0
        return true
    end
    #if one shared node, then recurse down and repeat
    return exists_hierarchical_tree3(shared12, shared13, shared23, num_levels)
end

""""""
function exists_hierarchical_tree3(
    shared12::Dict{Int, Tuple{Vararg{String}}},
    shared13::Dict{Int, Tuple{Vararg{String}}},
    shared23::Dict{Int, Tuple{Vararg{String}}},
    num_levels::Int,
    level::Int=1
)
    shared_nodes = Set()
    count = 0
    if haskey(shared12,level)
        push!(shared_nodes, shared12[level])
        count += 1
    end
    if haskey(shared13,level)
        push!(shared_nodes, shared13[level])
        count += 1
    end
    if haskey(shared23,level)
        push!(shared_nodes, shared23[level])
        count += 1
    end
    if length(shared_nodes) == 3
        return false
    end
    if length(shared_nodes) == 1 && count == 3 && level < num_levels-1
        println()
        return exists_hierarchical_tree3(shared12, shared13, shared23,
                                         num_levels, level+1)
    end
    return true
end

""""""
function hierarchical_tree_shared_nodes(
    district1::Dict{Tuple{Vararg{String}}, Any},
    district2::Dict{Tuple{Vararg{String}}, Any},
    num_levels::Int,
    level::Int=1,
    shared_nodes::Dict{Int,Tuple{Vararg{String}}}=Dict{Int,Tuple{Vararg{String}}}()
)
    num_of_shared = 0
    for node in keys(district1)
        if intact_huh(district1[node])
            continue
        end
        if node in keys(district2)
            shared_nodes[level] = node
            num_of_shared += 1
            if num_of_shared > 1
                return nothing
            end
            if level != num_levels-1
                check_finer = hierarchical_tree_shared_nodes(district1[node],
                                                             district2[node],
                                                             num_levels,
                                                             cur_level+1,
                                                             shared_nodes)
                if !check_finer
                    return nothing
                end
            end
        end
    end
    return shared_nodes
end


""""""
function edge_in_set(
    node_set::Dict{Tuple{Vararg{String}}, Any},
    edge::Union{Tuple, Vector},
    ns_d::Int,
    edge_d::Vector{Int},
    dists_to_edges::Dict
)
    ind = findall(x->x==ns_d, edge_d)
    if length(ind) == 2
        return strictly_contained_huh(node_set, Set(edge))
    elseif length(ind) == 1
        return strictly_contained_huh(node_set, Set([edge[ind[1]]]))
    end

    g = SimpleGraph(partition.num_dists)
    for dpair in keys(dists_to_edges)
        add_edge!(g, dpair[1], dpair[2])
    end
    nbrs = neighbors(g, ns_d)
    if length(nbrs) == 1
        new_pair = [dpair for dpair in keys(dists_to_edges) if ns_d in dpair]
        new_pair = new_pair[1]
        ind = findall(x->x==ns_d, new_pair)
        connecting_node = dists_to_edges[new_pair][ind]
        return strictly_contained_huh(node_set, Set([connecting_node]))
    elseif length(nbrs) == 2
        @assert edge_d[1] == edge_d[2] # o.w. we couldn't have a tree
                                       # since we know the edge_d would span the other two dists
        new_pair = [ns_d, edge_d[1]]
        if new_pair in keys(dists_to_edges)
            connecting_node = dists_to_edges[new_pair][1]
        else
            connecting_node = dists_to_edges[new_pair][2]
        end
        return strictly_contained_huh(node_set, Set([connecting_node]))
    end
    @assert false
end

""""""
function get_all_cuttable_edge_pairs!(
    cuttable_edgepairs::Set{Set{Tuple}},
    partition::MultiLevelPartition,
    multiscale_cuttable_trees::Dict{Int,MultiScaleCuttableTree},
    subgraphs::Dict{Int, MultiLevelSubGraph},
    constraints::Dict,
    dists_to_edges::Dict,
    exterior_decorators::Dict
)
    e1, e2 = collect(values(dists_to_edges))
    e1, e2 = Tuple(e1), Tuple(e2)
    push!(cuttable_edgepairs, Set([e1, e2]))

    first_cut_edges = Dict{Int, Set{Tuple}}()
    for (di, mstree) in multiscale_cuttable_trees
        first_cut_edges[di] = get_all_cuttable_edges(mstree)
    end

    min_pop = constraints[PopulationConstraint].min_pop
    max_pop = constraints[PopulationConstraint].max_pop
    total_pop = sum(get_total_graph_population(sg) for sg in values(subgraphs))

    edge_cuts = Dict{Tuple,Array}()
    for di in keys(first_cut_edges)
        for edge in first_cut_edges[di]
            multiscale_cuttable_tree = multiscale_cuttable_trees[di]
            subgraph = subgraphs[di]
            edge_cuts[edge] = get_cut_node_sets_w_pop(edge,
                                                      multiscale_cuttable_tree,
                                                      subgraph,
                                                      exterior_decorators[di])
        end
    end
    g = SimpleGraph(partition.num_dists)
    for dpair in keys(dists_to_edges)
        add_edge!(g, dpair[1], dpair[2])
    end
    for di in keys(first_cut_edges)
        for edge in first_cut_edges[di]
            to_cut_info = edge_cuts[edge]
            first_dist, first_dist_pop = to_cut_info[1]
            to_cut_set, to_cut_pop = to_cut_info[2]
            for dj in keys(first_cut_edges)
                for edge2 in first_cut_edges[dj]
                    if edge_in_set(first_dist, edge2, di, [dj, dj],
                                   dists_to_edges)
                        continue
                    end
                    if Set(edge) == Set(edge2)
                        continue
                    end
                    to_cut_info2 = edge_cuts[edge2]
                    second_dist, second_dist_pop = to_cut_info2[1]
                    if edge_in_set(second_dist, edge, dj, [di, di],
                                   dists_to_edges)
                        continue
                    end
                    third_dist_pop = total_pop - first_dist_pop
                    third_dist_pop -= second_dist_pop
                    if third_dist_pop > max_pop || third_dist_pop < min_pop
                        continue
                    end
                    push!(cuttable_edgepairs,Set([edge, edge2]))
                end
            end
            for (dpair, edge2) in dists_to_edges
                if edge_in_set(first_dist, edge2, di, dpair, dists_to_edges)
                    continue
                end
                if Set(edge) == Set(edge2)
                    continue
                end
            ########
                rem_edge!(g, dpair[1], dpair[2])
                pop_component1 = get_total_graph_population(subgraphs[dpair[1]])
                nbrs = neighbors(g, dpair[1])
                if length(nbrs) > 0
                    pop_component1 += sum(
                                       get_total_graph_population(subgraphs[dj])
                                       for dj in neighbors(g, dpair[1]))
                end
                pop_component2 = total_pop - pop_component1
                if pop_component1 < pop_component2
                    second_dist_pop = pop_component1
                    second_dist = dpair[1]
                else
                    second_dist_pop = pop_component2
                    second_dist = dpair[2]
                end
                add_edge!(g, dpair[1], dpair[2])
            ########
                if di == second_dist
                    continue
                end
                third_dist_pop = total_pop - first_dist_pop
                third_dist_pop -= second_dist_pop
                if third_dist_pop > max_pop || third_dist_pop < min_pop
                    continue
                end
                push!(cuttable_edgepairs,Set([edge, Tuple(edge2)]))
            end
        end
    end
end


""""""
function sample_merged_tree3(
    changed_districts::Vector{Int},
    partition::MultiLevelPartition,
    constraints::Dict,
    rng::AbstractRNG
)
    # prob_selecting_old_edge_to_remove
    D₁, D₂, D₃ = changed_districts
    dists_to_edges = sample_linking_edges(partition, D₁, D₂, D₃, rng)
    if nothing in values(dists_to_edges)
        return nothing, 0
    end

    graph = partition.graph

    node_sets = Dict(D₁=>partition.district_to_nodes[D₁],
                     D₂=>partition.district_to_nodes[D₂],
                     D₃=>partition.district_to_nodes[D₃])
    subgraphs = Dict(D₁=>MultiLevelSubGraph(partition.graph, node_sets[D₁]),
                     D₂=>MultiLevelSubGraph(partition.graph, node_sets[D₂]),
                     D₃=>MultiLevelSubGraph(partition.graph, node_sets[D₃]))

    dists_to_edges = Dict(p=>refine_linking_edge(e, subgraphs[p[1]],
                                                 subgraphs[p[2]], rng)
                          for (p,e) in dists_to_edges)

    num_levels = graph.num_levels
    multiscale_cuttable_trees = Dict(D₁=>MultiScaleCuttableTree(num_levels),
                                     D₂=>MultiScaleCuttableTree(num_levels),
                                     D₃=>MultiScaleCuttableTree(num_levels))
    pops = Dict(D₁=>partition.dist_populations[D₁],
                D₂=>partition.dist_populations[D₂],
                D₃=>partition.dist_populations[D₃])
    exterior_decorators = Dict(Dᵢ=>Dict{Tuple,Dict{Tuple,Int}}()
                               for Dᵢ in changed_districts)
    for (pair, edge) in dists_to_edges
        if edge[1] in keys(exterior_decorators[pair[1]])
            exterior_decorators[pair[1]][edge[1]][edge[2]] = pops[pair[2]]
        else
            exterior_decorators[pair[1]][edge[1]] = Dict(edge[2]=>pops[pair[2]])
        end
        if edge[2] in keys(exterior_decorators[pair[2]])
            exterior_decorators[pair[2]][edge[2]][edge[1]] = pops[pair[1]]
        else
            exterior_decorators[pair[2]][edge[2]] = Dict(edge[1]=>pops[pair[1]])
        end
    end
    # println("here1")
    construct_cuttable_tree!(multiscale_cuttable_trees[D₁], subgraphs[D₁],
                             constraints,
                             ext_decorators=exterior_decorators[D₁],
                             rng=rng, balance = (1,2))
    construct_cuttable_tree!(multiscale_cuttable_trees[D₂], subgraphs[D₂],
                             constraints,
                             ext_decorators=exterior_decorators[D₂],
                             rng=rng, balance = (1,2))
    construct_cuttable_tree!(multiscale_cuttable_trees[D₃], subgraphs[D₃],
                             constraints,
                             ext_decorators=exterior_decorators[D₃],
                             rng=rng, balance = (1,2))
    cuttable_edgepairs = Set{Set{Tuple}}()
    get_all_cuttable_edge_pairs!(cuttable_edgepairs, partition,
                                 multiscale_cuttable_trees, subgraphs,
                                 constraints, dists_to_edges,
                                 exterior_decorators)
    edges_prob = 1.0/length(cuttable_edgepairs)
    old_tree = ([multiscale_cuttable_trees[D₁], multiscale_cuttable_trees[D₂],
                 multiscale_cuttable_trees[D₃]],
                 dists_to_edges)
    return old_tree, edges_prob
end

""""""
function build_forest_recom3(
    constraints::Dict
)
    f(p, m, r) = forest_recom3!(p, m, constraints, r)
    return f
end


""""""
function sample_linking_edg3(
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
                if nbr in keys(node_set₂)
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
            if !intact_huh(node_set₁[node]) && node in keys(node_set₂)
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
function get_shared_features(partition::MultiLevelPartition,
    D₁::Integer,
    D₂::Integer,
    D₃::Integer
)
    district1 = partition.district_to_nodes[D₁]
    district2 = partition.district_to_nodes[D₂]
    district3 = partition.district_to_nodes[D₃]
    num_levels = partition.graph.num_levels
    shared_features = Dict{Set{Int}, Dict{Int,Tuple{Vararg{String}}}}()
    shared12 = hierarchical_tree_shared_nodes(district1, district2, num_levels)
    shared_features[Set([1,2])] = shared12
    shared13 = hierarchical_tree_shared_nodes(district1, district3, num_levels)
    shared_features[Set([1,3])] = shared13
    shared23 = hierarchical_tree_shared_nodes(district2, district3, num_levels)
    shared_features[Set([2,3])] = shared23
    return [shared12, shared13, shared23]
end

""""""
function get_finest_shared_features(
    sharednodelist::Vector{Dict{Int,Tuple{Vararg{String}}}}
)
    finest_shared_features = Dict{Tuple{Vararg{String}},Set{Int}}()
    for i = 1:length(sharednodelist)
        l = maximum(keys(sharednodelist[i]))
        for node in sharednodelist[i][l]
            if !haskey(finest_shared_features,node)
                finest_shared_features[node] = Set([i])
            else
                push!(finest_shared_features[node],i)
            end
        end
    end
    return finest_shared_features
end

""""""
function sample_linking_edge(
    partition::MultiLevelPartition,
    district_pairs::Vector{Vector{Int}},
    rng::AbstractRNG
)
    num_levels = partition.graph.num_levels
    lvl_to_wght = Dict{Int, Real}(l=>0 for l=1:num_levels)
    for level = 1:num_levels
        for pair in district_pairs
            di, dj = collect(pair)
            lvl_to_wght[level] += partition.cross_district_edges[di, dj, level]
        end
    end
    finest_level = maximum([l for l=1:num_levels if lvl_to_wght[l]>0])
    choice = rand(rng)*lvl_to_wght[finest_level]
    cum_sum = 0
    chosen_pair = nothing
    for pair in district_pairs
        di, dj = collect(pair)
        cum_sum += partition.cross_district_edges[di, dj, finest_level]
        if cum_sum > choice
            edge = sample_linking_edge(partition, di, dj, rng)
            chosen_pair = pair
            return edge, pair
        end
    end
    @assert false
end
""""""
function sample_linking_edges(
    partition::MultiLevelPartition,
    D₁::Integer, D₂::Integer, D₃::Integer,
    rng::AbstractRNG
)
    districts_ps = Combinatorics.powerset([D₁, D₂, D₃])
    district_pairs = [e for e in districts_ps if length(e) == 2]
    edge1, pair1 = sample_linking_edge(partition, district_pairs, rng)

    district_pairs = [e for e in district_pairs if e != pair1]
    edge2, pair2 = sample_linking_edge(partition, district_pairs, rng)

    return Dict(pair1=>edge1, pair2=>edge2)
end
