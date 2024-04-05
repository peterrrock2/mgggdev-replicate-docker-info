""""""
function attemptConstructPartition(
    graph::MultiLevelGraph,
    constraints::Dict,
    num_dists::Int, rng::AbstractRNG,
    max_attempts::Int,
    district_to_nodes::Vector{Dict{Tuple{Vararg{String}}, Any}} = [],
    remaining_nodes::Dict{Tuple{Vararg{String}}, Any} = 
                                              Dict{Tuple{Vararg{String}}, Any}()
)
    if length(remaining_nodes) == 0
        for node in graph.id_to_partitions[1]
            remaining_nodes[node] = nothing # Dict{Tuple{Vararg{String}}}()
        end
    end
    # println("remaining_nodes length ", length(remaining_nodes))

    for district = length(district_to_nodes)+1:num_dists-1
        found_district = false
        for attempt = 1:max_attempts
            subgraph = MultiLevelSubGraph(graph, remaining_nodes)
            bal = (1, num_dists-district)
            
            multiscale_cuttable_tree = MultiScaleCuttableTree(graph.num_levels)

            # try
            construct_cuttable_tree!(multiscale_cuttable_tree, subgraph,
                                     constraints, rng=rng, balance=bal)
            proposed_cut = cut_edge(multiscale_cuttable_tree, subgraph, rng)
            node_sets_w_pops, edge, prob_edge = proposed_cut

            if edge == nothing
                continue
            end
            # district_to_nodes[district] = node_sets_w_pops[1][1][()]
            push!(district_to_nodes, node_sets_w_pops[1][1][()])
            if district == num_dists-1
                # district_to_nodes[district+1] = node_sets_w_pops[2][1][()]
                push!(district_to_nodes, node_sets_w_pops[2][1][()])
                tot_dists = num_dists
            else
                remaining_nodes_tmp = node_sets_w_pops[2][1][()]
                # remaining_nodes = node_sets_w_pops[2][1][()]
                tot_dists = district
            end

            pop_const = constraints[PopulationConstraint]
            ideal_pop = (pop_const.min_pop + pop_const.max_pop)*0.5
            new_dist_pop_dev = round(100*(node_sets_w_pops[1][2]-ideal_pop)/ideal_pop, digits = 3)
            rem_dists = num_dists-district
            rem_avg_pop_dev = round(100*(node_sets_w_pops[2][2]/rem_dists-ideal_pop)/(ideal_pop), digits = 3)
            
            cur_districts = view(district_to_nodes, 1:tot_dists)
            if !satisfies_constraints(graph, cur_districts, num_dists,
                                      constraints)
                continue
            elseif district < num_dists-1
                remaining_nodes = remaining_nodes_tmp
            end

            # if success continue with remaining graph
            found_district = true
            break
        end
        if !found_district
            return nothing, false
        end
    end

    success = true::Bool

    node_to_district = construct_node_map(district_to_nodes)

    dist_populations = Vector{Int}(undef, num_dists)
    subgraphs = Vector{MultiLevelSubGraph}(undef, num_dists)
    for di = 1:num_dists
        subgraph = MultiLevelSubGraph(graph, district_to_nodes[di])
        dist_populations[di] = get_total_graph_population(subgraph)
        subgraphs[di] = subgraph
    end
    num_levels = graph.num_levels
    cross_district_edges = Array{Real}(undef,
                                       (num_dists, num_dists, num_levels))
    parent = nothing
    extensions = Dict{String,Any}()
    proposed_extensions = Dict{String,Any}()

    partition = MultiLevelPartition(num_dists, cross_district_edges,
                                    node_to_district, district_to_nodes,
                                    subgraphs, dist_populations, graph, parent,
                                    extensions, proposed_extensions)

    set_cross_district_edges!(partition)

    return partition, success
end

""""""
function constructMultiLevelPartition(
    multi_level_graph::MultiLevelGraph,
    constraints::Dict,
    num_dists::Int,
    rng::AbstractRNG,
    max_attempts::Int,
    district_to_nodes::Vector{Dict{Tuple{Vararg{String}}, Any}} = 
                             Vector{Dict{Tuple{Vararg{String}}, Any}}(undef, 0),
    remaining_nodes::Dict{Tuple{Vararg{String}}, Any} = 
                                              Dict{Tuple{Vararg{String}}, Any}()
)::MultiLevelPartition

    d2ns = deepcopy(district_to_nodes)
    rem_n = deepcopy(remaining_nodes)
    partition, success = attemptConstructPartition(multi_level_graph,
                                                   constraints, num_dists,
                                                   rng, max_attempts,
                                                   d2ns, rem_n)
    for attempt = 1:max_attempts
        if success
            break
        end
        d2ns = deepcopy(district_to_nodes)
        rem_n = deepcopy(remaining_nodes)
        partition, success = attemptConstructPartition(multi_level_graph,
                                                       constraints,
                                                       num_dists, rng,
                                                       max_attempts,
                                                       d2ns, rem_n)
    end

    if !success
        throw(
        DomainError(
            success,
            "Failed to construct a random partition after "*
            string(max_attempts)*" attempts."
            )
        )
    end

    return partition
end

""""""
function MultiLevelPartition(
    multi_level_graph::MultiLevelGraph,
    constraints::Dict,
    num_dists::Int;
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64),
    max_attempts::Int=100
)::MultiLevelPartition
        partition = constructMultiLevelPartition(multi_level_graph,
                                                 constraints, num_dists, rng,
                                                 max_attempts)
        return partition
end

""""""
function MultiLevelPartition(
    multi_level_graph::MultiLevelGraph,
    constraints::Dict,
    clusters::Dict{String, Any};
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64),
    max_attempts::Int=100,
    node_key::String="nodes",
    cluster_key::String="clusters",
    districts_key::String="districts"
)::MultiLevelPartition
    num_levels = multi_level_graph.num_levels
    levels = multi_level_graph.levels
    base_graph = multi_level_graph.graphs_by_level[num_levels]
    node_to_district = Dict{Tuple{Vararg{String}}, Int}()
    dist_index_shift = 0
    for cluster in clusters[cluster_key]
        println(cluster)
        single_cluster = Dict(cluster_key=>[cluster])
        single_cluster_graph = cluster_base_graph(base_graph, single_cluster,
                                                  node_key=node_key,
                                                  cluster_key=cluster_key)
        multi_level_cluster_graph = MultiLevelGraph(single_cluster_graph,
                                                    levels)
        num_dists = cluster[districts_key]
        update_node_to_district = Dict{Tuple{Vararg{String}}, Int}()
        if num_dists == 1
            graph = multi_level_cluster_graph.graphs_by_level[1]
            fine_nodes = multi_level_cluster_graph.id_to_partitions[num_levels]
            for node in fine_nodes
                update_node_to_district[node] = dist_index_shift+1
            end
        else
            single_cluster_partition = constructMultiLevelPartition(
                                            multi_level_cluster_graph,
                                            constraints, num_dists, rng,
                                            max_attempts)
            n2d = single_cluster_partition.node_to_district
            update_node_to_district = Dict{Tuple{Vararg{String}}, Int}(
                                         n=>d+dist_index_shift for (n,d) in n2d)
        end

        merge!(node_to_district, update_node_to_district)
        dist_index_shift += num_dists
    end

    return MultiLevelPartition(multi_level_graph, node_to_district)
end

""""""
function MultiLevelPartition(
    multi_level_graph::MultiLevelGraph,
    constraints::Dict,
    assignment_col::AbstractString,
    fixed_dists::Vector;
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64),
    max_attempts::Int=100
)::MultiLevelPartition
    partition = MultiLevelPartition(multi_level_graph, assignment_col)
    district_to_nodes = [partition.district_to_nodes[d] 
                         for d=1:partition.num_dists if d in fixed_dists]
    
    fixed_nodes = merge_nodesets(multi_level_graph, district_to_nodes)
    
    remaining_nodes = Dict{Tuple{Vararg{String}}, Any}()
    for node in multi_level_graph.id_to_partitions[1]
        remaining_nodes[node] = nothing # Dict{Tuple{Vararg{String}}}()
    end 
    remaining_nodes = get_node_set_complement(remaining_nodes, fixed_nodes, 
                                              multi_level_graph)
    
    num_dists = partition.num_dists
    partition = constructMultiLevelPartition(multi_level_graph, constraints, 
                                             num_dists, rng, max_attempts,
                                             district_to_nodes, remaining_nodes)
    return partition

end
