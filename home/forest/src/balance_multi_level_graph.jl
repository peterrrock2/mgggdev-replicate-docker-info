mutable struct MultiScaleForest
    trees::Vector{Dict{Tuple{Vararg{String}},SimpleWeightedGraph}}
    vmaps::Vector{Dict{Tuple{Vararg{String}},Vector{Int}}}
end

struct CuttableTree
    tree::SimpleWeightedGraph
    vmap::Vector{Int}
    enter::Vector{Int}
    exit::Vector{Int}
    parent::Vector{Int}
    name_to_index::Dict{Tuple{Vararg{String}}, Int}
    pop_weight::Vector{Int}
    cuttable_edges::Vector{Tuple}
    cuttable_nodes::Vector{Tuple}
end

mutable struct MultiScaleCuttableTree
    cuttable_trees::Vector{Dict{Tuple,CuttableTree}}
    specified_edges::Vector{Dict{Tuple,Dict{Set{Int},Set{Tuple}}}}
end


""""""
function MultiScaleCuttableTree(num_levels::Int)
    cuttable_trees = [Dict{Tuple,CuttableTree}() for ii = 1:num_levels]
    specified_edges = [Dict{Tuple,Dict{Set{Int},Set{Tuple}}}()
                       for ii = 1:num_levels]
    return MultiScaleCuttableTree(cuttable_trees, specified_edges)
end


"""
    Update curr's pop to indicate pop of entire subtree rooted here.
"""
function determine_subtree_pop!(
    pop_weight::Vector{Int},
    parent::Vector{Int},
    subtree::SimpleWeightedGraph,
    vmap::Vector{Int},
    subgraph::MultiLevelSubGraph,
    edge_decorators::Dict{Tuple,Dict{Tuple,Int}},
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}},
    curr::Int,
    level::Int
)
    curr_name = subgraph.parent.id_to_partitions[level][vmap[curr]]

    pop_w = get_node_population(subgraph, curr_name)

    if haskey(edge_decorators, curr_name)
        for pop in values(edge_decorators[curr_name])
            pop_w += pop
        end
    end

    if haskey(ext_decorators, curr_name)
        for pop in values(ext_decorators[curr_name])
            pop_w += pop
        end
    end

    for nbr in neighbors(subtree, curr)
        if nbr != parent[curr]
            pop_w += pop_weight[nbr]
        end
    end
    pop_weight[curr] = pop_w
end

""""""
function get_neighbors_pop(
    subtree::SimpleWeightedGraph,
    vmap::Vector{Int},
    subgraph::MultiLevelSubGraph,
    edge_decorators::Dict{Tuple,Dict{Tuple,Int}},
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}},
    pop_weight::Vector{Int},
    parent::Vector{Int},
    curr::Int,
    level::Int
)::Dict{Tuple,Int}
    total_pop = get_total_population(subgraph, ext_decorators)

    nbr_cut_pops = Dict{Tuple,Int}()
    curr_name = subgraph.parent.id_to_partitions[level][vmap[curr]]
    if haskey(edge_decorators, curr_name)
        merge!(nbr_cut_pops, edge_decorators[curr_name])
    end
    if haskey(ext_decorators, curr_name)
        merge!(nbr_cut_pops, ext_decorators[curr_name])
    end
    for nbr in neighbors(subtree, curr)
        nbr_name = subgraph.parent.id_to_partitions[level][vmap[nbr]]
        if parent[curr] == nbr
            nbr_cut_pops[nbr_name] = total_pop - pop_weight[curr]
        else
            nbr_cut_pops[nbr_name] = pop_weight[nbr]
        end
    end
    return nbr_cut_pops
end

""""""
function check_cuttable_node(
    subtree::SimpleWeightedGraph,
    vmap::Vector{Int},
    subgraph::MultiLevelSubGraph,
    edge_decorators::Dict{Tuple,Dict{Tuple,Int}},
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}},
    constraints::Dict,
    pop_weight::Vector{Int},
    parent::Vector{Int},
    curr::Int,
    level::Int,
    balance::Tuple{Int,Int}
)
    if level == subgraph.parent.num_levels
        return false
    end

    nbr_cut_pops = get_neighbors_pop(subtree, vmap, subgraph, edge_decorators,
                                     ext_decorators, pop_weight, parent, curr,
                                     level)

    # get possible node combinations -- will be based on orient neighbors
    nbrs = collect(keys(nbr_cut_pops))
    sub_nbrs = view(nbrs, 1:length(nbrs)-1)
    nbr_combinations = Combinatorics.powerset(sub_nbrs)

    min_pop = constraints[PopulationConstraint].min_pop
    max_pop = constraints[PopulationConstraint].max_pop

    for nbr in nbrs
        if nbr_cut_pops[nbr] > max(balance[1], balance[2])*max_pop
            return false
        end
    end

    curr_pop = get_node_population(subgraph, level, vmap[curr])
    min_pops_1, min_pops_2 = 0, 0
    min1, min2, max1, max2 = false, false, false, false
    cuttable = false
    tot_nbr_pop = 0
    for n in nbrs
        tot_nbr_pop += nbr_cut_pops[n]
    end
    for nbr_cmb in nbr_combinations
        min_pops_1 = 0
        for n in nbr_cmb
            min_pops_1 += nbr_cut_pops[n]
        end
        min_pops_2 = tot_nbr_pop - min_pops_1
        # min_pops_1 = sum([nbr_cut_pops[n] for n in nbr_cmb])
        # min_pops_2 = sum([nbr_cut_pops[n] for n in nbrs if n ∉ nbr_cmb])

        min1 = (balance[1]*min_pop <= min_pops_1 + curr_pop)
        min2 = (balance[2]*min_pop <= min_pops_2 + curr_pop)
        max1 = (balance[1]*max_pop >  min_pops_1)
        max2 = (balance[2]*max_pop >  min_pops_2)
        cuttable = (min1 && max1 && min2 && max2)

        if balance[1] != balance[2]
            min1 = (balance[2]*min_pop <= min_pops_1 + curr_pop)
            min2 = (balance[1]*min_pop <= min_pops_2 + curr_pop)
            max1 = (balance[2]*max_pop >  min_pops_1)
            max2 = (balance[1]*max_pop >  min_pops_2)
            cuttable = cuttable || (min1 && max1 && min2 && max2)
        end

        if cuttable
            # curr_name = subgraph.parent.id_to_partitions[level][vmap[curr]]
            return cuttable # true
        end
    end
    return cuttable # false
end


""""""
function check_cuttable_edge_huh(
    curr::Int,
    pop_weight::Vector{Int},
    total_pop::Int,
    constraints::Dict,
    balance::Tuple{Int, Int}
)
    pop_1 = pop_weight[curr]
    pop_2 = total_pop - pop_1

    min_pop = constraints[PopulationConstraint].min_pop
    max_pop = constraints[PopulationConstraint].max_pop

    check_11 = balance[1]*min_pop <= pop_1
    check_11 = check_11 && pop_1 <= balance[1]*max_pop
    check_22 = balance[2]*min_pop <= pop_2
    check_22 = check_22 && pop_2 <= balance[2]*max_pop

    check_12 = balance[2]*min_pop <= pop_1
    check_12 = check_12 && pop_1 <= balance[2]*max_pop
    check_21 = balance[1]*min_pop <= pop_2
    check_21 = check_21 && pop_2 <= balance[1]*max_pop

    return (check_11 && check_22) || (check_12 && check_21)
end

""""""
function get_total_population(
    subgraph::MultiLevelSubGraph,
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}}
)
    total_pop = get_total_graph_population(subgraph)
    for (k,d) in ext_decorators
        total_pop += sum(values(d))
    end
    return total_pop
end

""""""
function find_cuttable_components(
    subtree::SimpleWeightedGraph,
    vmap::Vector{Int},
    subgraph::MultiLevelSubGraph,
    constraints::Dict,
    edge_decorators::Dict{Tuple,Dict{Tuple,Int}},
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}},
    level::Int,
    balance::Tuple{Int,Int}=(1,1)
)::CuttableTree
    total_pop = get_total_population(subgraph, ext_decorators)
    topological_sort_ind = 0

    visited = [false for ii = 1:nv(subtree)]
    parent = [0 for ii = 1:nv(subtree)]

    enter = [-1 for ii = 1:nv(subtree)]
    exit = [-1 for ii = 1:nv(subtree)]
    pop_weight = [-1 for ii = 1:nv(subtree)]

    name_to_index = Dict{Tuple{Vararg{String}}, Int}()
    cuttable_edges = Vector{Tuple}(undef, 0)
    cuttable_nodes = Vector{Tuple}(undef, 0)

    stack = [1] # arbitrarily root at 1

    while length(stack) > 0
        curr = pop!(stack);

        if !visited[curr] # first visit
            visited[curr] = true
            enter[curr] = topological_sort_ind
            topological_sort_ind += 1

            push!(stack, curr) # Revisit after children
            for nbr in neighbors(subtree, curr)
                if nbr != parent[curr]
                    parent[nbr] = curr
                    push!(stack, nbr)
                end
            end
        else # second visit
            exit[curr] = topological_sort_ind
            topological_sort_ind += 1

            determine_subtree_pop!(pop_weight, parent, subtree, vmap,
                                   subgraph, edge_decorators, ext_decorators,
                                   curr, level)

            cuttable_node_huh = check_cuttable_node(subtree, vmap, subgraph,
                                                    edge_decorators,
                                                    ext_decorators,
                                                    constraints, pop_weight,
                                                    parent, curr, level,
                                                    balance)
            if cuttable_node_huh
                graph = subgraph.parent
                curr_name = graph.id_to_partitions[level][vmap[curr]]
                name_to_index[curr_name] = curr
                for nbr in neighbors(subtree, curr)
                    nbr_name = graph.id_to_partitions[level][vmap[nbr]]
                    name_to_index[nbr_name] = nbr
                end
                push!(cuttable_nodes, curr_name)
            end

            cuttable_edge_huh = check_cuttable_edge_huh(curr, pop_weight,
                                                        total_pop, constraints,
                                                        balance)
            if cuttable_edge_huh && curr != 1
                graph = subgraph.parent
                curr_name = graph.id_to_partitions[level][vmap[curr]]
                parent_name = graph.id_to_partitions[level][vmap[parent[curr]]]
                cut = (curr_name, parent_name)
                name_to_index[curr_name] = curr
                push!(cuttable_edges, cut)
            end
        end
    end

    return CuttableTree(subtree, vmap, enter, exit, parent, name_to_index,
                        pop_weight, cuttable_edges, cuttable_nodes)
end

""""""
function get_nbr_name_in_edge(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelSubGraph,
    node::Tuple,
    edge_tuple::Tuple{Int, Int}
)
    edge_key = Set(edge_tuple)
    level = length(node)
    key = node[1:level-1]
    specified_edges = multiscale_cuttable_tree.specified_edges[level][key]

    if haskey(specified_edges, edge_key)
        edge = specified_edges[edge_key]
        specified_edges = multiscale_cuttable_tree.specified_edges[level][key]
        return collect(setdiff(edge,Set([node])))[1]
    else
        nbr_id = edge_tuple[2]
        graph = subgraph.parent
        cuttable_tree = multiscale_cuttable_tree.cuttable_trees[level][key]
        return graph.id_to_partitions[level][cuttable_tree.vmap[nbr_id]]
    end
end

""""""
function refine_node(
    subgraph::MultiLevelSubGraph,
    node::Tuple,
    nbr::Tuple,
    rng::AbstractRNG,
)::Tuple
    level = length(node)
    graph = subgraph.parent
    edge_weight = get_edge_weight(subgraph, node, nbr)
    rnd_num = rand(rng)*edge_weight
    cum_edge_weight = 0
    fine_node_ids_near_nbr = get_fine_neighbors(graph, nbr, node)
    for fine_node_id_near_nbr in fine_node_ids_near_nbr
        fine_node = graph.id_to_partitions[level+1][fine_node_id_near_nbr]
        if !intersects_huh(subgraph, fine_node)
            continue
        end
        cum_edge_weight += get_edge_weight(subgraph, fine_node, nbr)
        if cum_edge_weight > rnd_num
            return fine_node
        end
    end
end

""""""
function add_decorator!(
    sub_decorators::Dict{Tuple,Dict{Tuple,Int}},
    cuttable_tree::CuttableTree,
    subgraph::MultiLevelSubGraph,
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}},
    fine_node::Tuple,
    nbr::Tuple,
    node_id::Int,
    nbr_id::Int
)
    pop = 0
    if nbr_id == cuttable_tree.parent[node_id]
        pop = get_total_population(subgraph, ext_decorators)
        pop -= cuttable_tree.pop_weight[node_id]
    else
        pop = cuttable_tree.pop_weight[nbr_id]
    end
    if haskey(sub_decorators, fine_node)
        sub_decorators[fine_node][nbr] = pop
    else
        sub_decorators[fine_node] = Dict(nbr=>pop)
    end
end


""""""
function add_decorator!(
    sub_decorators::Dict{Tuple,Dict{Tuple,Int}},
    cuttable_tree::CuttableTree,
    fine_node::Tuple,
    nbr::Tuple,
    pop::Int
)
    if haskey(sub_decorators, fine_node)
        sub_decorators[fine_node][nbr] = pop
    else
        sub_decorators[fine_node] = Dict(nbr=>pop)
    end
end


""""""
function specify_edges!(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelSubGraph,
    vmap_t2sg::Vector{Int},
    edge_decorators::Dict{Tuple,Dict{Tuple,Int}},
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}},
    node::Tuple,
    rng::AbstractRNG
)
    sub_decorators = Dict{Tuple,Dict{Tuple,Int}}()
    graph = subgraph.parent
    level = length(node)
    key = node[1:level-1]
    cuttable_tree = multiscale_cuttable_tree.cuttable_trees[level][key]
    specified_edges = multiscale_cuttable_tree.specified_edges[level][key]
    node_id = cuttable_tree.name_to_index[node]

    for nbr_id in neighbors(cuttable_tree.tree, node_id)
        nbr = get_nbr_name_in_edge(multiscale_cuttable_tree, subgraph, node,
                                   (node_id, nbr_id))
        # fine_node = nothing
        # try 
            fine_node = refine_node(subgraph, node, nbr, rng)
        # catch e 
        #     level = length(node)
        #     graph = subgraph.parent
        #     edge_weight = get_edge_weight(subgraph, node, nbr)
        #     @show node, nbr
        #     @show level, edge_weight
        #     rnd_num = rand(rng)*edge_weight
        #     @show rnd_num
        #     cum_edge_weight = 0
        #     fine_node_ids_near_nbr = get_fine_neighbors(graph, nbr, node)
        #     @show fine_node_ids_near_nbr
        #     for fine_node_id_near_nbr in fine_node_ids_near_nbr
        #         fine_node = graph.id_to_partitions[level+1][fine_node_id_near_nbr]
        #         @show fine_node
        #         if !intersects_huh(subgraph, fine_node)
        #             continue
        #         end
        #         edge_weight = get_edge_weight(subgraph, fine_node, nbr)
        #         @show fine_node, nbr, edge_weight
        #         cum_edge_weight += get_edge_weight(subgraph, fine_node, nbr)
        #         if cum_edge_weight > rnd_num
        #             break
        #         end
        #     end
        #     @show subgraph.modified_nbr_weights
        #     @assert false
        # end
        specified_edges[Set([node_id, nbr_id])] = Set([fine_node, nbr])
        add_decorator!(sub_decorators, cuttable_tree, subgraph, ext_decorators,
                       fine_node, nbr, node_id, nbr_id)
    end

    extended_decorators = []
    if haskey(edge_decorators, node)
        extended_decorators = edge_decorators[node]
    end

    for (nbr, pop) in extended_decorators
        key = get_edge_key(node, nbr)
        level = length(key)+1
        cuttable_tree = multiscale_cuttable_tree.cuttable_trees[level][key]
        specified_edges = multiscale_cuttable_tree.specified_edges[level][key]
        fine_node = refine_node(subgraph, node, nbr, rng)
        node_id = cuttable_tree.name_to_index[node[1:level]]
        nbr_id = cuttable_tree.name_to_index[nbr[1:level]]
        specified_edges[Set([node_id, nbr_id])] = Set([fine_node, nbr])
        add_decorator!(sub_decorators, cuttable_tree, fine_node, nbr, pop)
    end

    return sub_decorators
end

""""""
function get_ext_decorators_lvl(
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}},
    level::Int
)
    if length(ext_decorators) == 0
        return ext_decorators
    end
    ext_decorators_lvl = Dict{Tuple,Dict{Tuple,Int}}()
    for (k,v) in ext_decorators
        key = k[1:level]
        if key ∉ keys(ext_decorators_lvl)
            ext_decorators_lvl[key] = v
        else
            merge!(ext_decorators_lvl[key], v)
        end
    end
    return ext_decorators_lvl
end


""""""
function sample_spanning_tree(
    simple_graph::SimpleWeightedGraph,
    rng::AbstractRNG
)
    if nv(simple_graph) == 1
        subtree, vmap_t2sg = induced_subgraph(simple_graph, [1])
    else
        subtree_edges = wilson_rst(simple_graph, rng)
        subtree, vmap_t2sg = induced_subgraph(simple_graph, subtree_edges)
    end
    return subtree, vmap_t2sg
end

""""""
function construct_cuttable_tree!(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelSubGraph,
    constraints::Dict,
    edge_decorators::Dict{Tuple,Dict{Tuple,Int}}=Dict{Tuple,Dict{Tuple,Int}}();
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}}=Dict{Tuple,Dict{Tuple,Int}}(),
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64),
    key::Tuple=(),
    level::Int=1,
    balance::Tuple=(1,1)
)
    simple_graph, vmap_sg2g = get_simple_subgraph(subgraph, level, key)
    subtree, vmap_t2sg = sample_spanning_tree(simple_graph, rng)
    vmap_t2g = [vmap_sg2g[vmap_t2sg[ii]] for ii = 1:nv(subtree)] # global index
    ext_decorators_lvl = get_ext_decorators_lvl(ext_decorators, level)
    cuttable_tree = find_cuttable_components(subtree, vmap_t2g, subgraph,
                                             constraints, edge_decorators,
                                             ext_decorators_lvl, level,
                                             balance)
    multiscale_cuttable_tree.cuttable_trees[level][key] = cuttable_tree
    multiscale_cuttable_tree.specified_edges[level][key] =
                                                     Dict{Set{Int},Set{Tuple}}()

    for node in cuttable_tree.cuttable_nodes
        sub_edge_decorators = specify_edges!(multiscale_cuttable_tree, subgraph,
                                             vmap_t2sg, edge_decorators,
                                             ext_decorators, node, rng)
        construct_cuttable_tree!(multiscale_cuttable_tree, subgraph,
                                 constraints, sub_edge_decorators,
                                 ext_decorators=ext_decorators, rng=rng,
                                 key=node, level=level+1, balance=balance)
    end
end

""""""
function get_all_cuttable_edges!(
    cuttable_edges::Set{Tuple},
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    level::Int = 1
)
    for key in keys(multiscale_cuttable_tree.cuttable_trees[level])
        ce = multiscale_cuttable_tree.cuttable_trees[level][key].cuttable_edges
        union!(cuttable_edges, ce)
        if level < length(multiscale_cuttable_tree.cuttable_trees)
            get_all_cuttable_edges!(cuttable_edges, multiscale_cuttable_tree,
                                    level+1)
        end
    end
end

""""""
function get_all_cuttable_edges(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    level::Int = 1
)
    cuttable_edges = Set{Tuple}()
    for key in keys(multiscale_cuttable_tree.cuttable_trees[level])
        ce = multiscale_cuttable_tree.cuttable_trees[level][key].cuttable_edges
        union!(cuttable_edges, ce)
        if level < length(multiscale_cuttable_tree.cuttable_trees)
            get_all_cuttable_edges!(cuttable_edges, multiscale_cuttable_tree,
                                    level+1)
        end
    end
    return cuttable_edges
end

""""""
function get_all_cuttable_edge_pairs!(
    cuttable_edges::Set{Set{Tuple}},
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelSubGraph,
    constraints::Dict
)
    first_cuttable_edges = Set{Tuple}()
    get_all_cuttable_edges!(first_cuttable_edges, multiscale_cuttable_tree)
    min_pop = constraints[PopulationConstraint].min_pop
    max_pop = constraints[PopulationConstraint].max_pop
    total_pop = get_total_graph_population(subgraph)
    edge_cuts = Dict{Tuple,Array}()
    for edge in first_cuttable_edges
        edge_cuts[edge] = get_cut_node_sets_w_pop(edge, 
                                                  multiscale_cuttable_tree, 
                                                  subgraph)
    end
    for edge in first_cuttable_edges
        to_cut_info = edge_cuts[edge]
        first_dist, first_dist_pop = to_cut_info[1]
        to_cut_set, to_cut_pop = to_cut_info[2]
        for edge2 in first_cuttable_edges  # first_cuttable_edges(edge+1:end)
            if strictly_contained_huh(first_dist, Set(edge2))
                continue
            end
            if Set(edge) == Set(edge2)
                continue
            end
            to_cut_info2 = edge_cuts[edge2]
            second_dist, second_dist_pop = to_cut_info2[1]
            if strictly_contained_huh(second_dist,Set(edge))
                continue
            end
            third_dist_pop = total_pop - first_dist_pop - second_dist_pop
            if third_dist_pop > max_pop || third_dist_pop < min_pop
                continue
            end
            push!(cuttable_edges, Set([edge, edge2]))
        end
    end
end

""""""
function get_probability_of_edge(
    cuttable_edges::Set{Tuple},
    # linking_edge::Set{Tuple}
)
    return 1.0/length(cuttable_edges)
end

""""""
function choose_random_cuttable_edge(
    cuttable_edges::Set,
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    rng::AbstractRNG
)
    if length(cuttable_edges)==0
        return nothing, nothing
    end
    rnd_ind = Int(ceil(rand(rng)*length(cuttable_edges)))
    return collect(cuttable_edges)[rnd_ind], 1.0/length(cuttable_edges)
end

""""""
function get_coarsened_node_sets(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelSubGraph,
    edge::Tuple,
    node_set1::Dict{Tuple{Vararg{String}}, Any},
    node_set2::Dict{Tuple{Vararg{String}}, Any},
    key::Tuple
)
    if length(key) == 0
        return node_set1, node_set2
    end
    level = length(key)
    key = key[1:level-1]
    cuttable_tree = multiscale_cuttable_tree.cuttable_trees[level][key]
    n1, n2 = edge
    n1_id = cuttable_tree.name_to_index[n1[1:level]]
    n2_id = cuttable_tree.parent[n1_id]

    specified_edges = multiscale_cuttable_tree.specified_edges[level][key]

###### determine conditions -- should be in a function
    node_set1_conditions = Dict{Int, Tuple{Int, Int}}()
    node_set2_conditions = Dict{Int, Tuple{Int, Int}}()
    for nbr_id in neighbors(cuttable_tree.tree, n1_id)
        specified_edge = specified_edges[Set([n1_id, nbr_id])]
        specified_node = [n for n in specified_edge
                          if n[1:level] == n1[1:level]][1]
        specified_nbr = [n for n in specified_edge
                         if n[1:level] != n1[1:level]][1]
        if nbr_id == cuttable_tree.parent[n1_id]
            enter = cuttable_tree.enter[n1_id]
            exit = cuttable_tree.exit[n1_id]
        else
            enter = cuttable_tree.enter[nbr_id]
            exit = cuttable_tree.exit[nbr_id]
        end

        if intersects_huh(node_set1, specified_node)
            node_set1_conditions[nbr_id] = (enter, exit)
        elseif intersects_huh(node_set2, specified_node)
            node_set2_conditions[nbr_id] = (enter, exit)
        else
            throw(
            DomainError(
                false,
                "Failed to find "*string(specified_node)*" in either node set."

                )
            )
        end
    end
#######

    node_set1 = Dict{Tuple{Vararg{String}}, Any}(key => node_set1)
    node_set2 = Dict{Tuple{Vararg{String}}, Any}(key => node_set2)

    for n_id = 1:nv(cuttable_tree.tree)
        n_enter = cuttable_tree.enter[n_id]
        n_exit = cuttable_tree.exit[n_id]
        g_n_id = cuttable_tree.vmap[n_id]
        node = subgraph.parent.id_to_partitions[level][g_n_id]
        in1 = false
########### determine if in node set -- should be in a function
        for (nbr_id, bnds) in node_set1_conditions
            if nbr_id == cuttable_tree.parent[n1_id]
                if n_enter < bnds[1] || n_exit > bnds[2]
                    in1 = true
                    break
                end
            else
                if n_enter >= bnds[1] && n_exit <= bnds[2]
                    in1 = true
                    break
                end
            end
        end
###########
        node = subgraph.parent.id_to_partitions[level][cuttable_tree.vmap[n_id]]
        if node == n1[1:level]
            continue
        end
        if in1
            node_set1[key][node] = get_node_set(subgraph, node)
        else
            node_set2[key][node] = get_node_set(subgraph, node)
        end
    end

    return get_coarsened_node_sets(multiscale_cuttable_tree, subgraph, edge,
                                   node_set1, node_set2, key)
end

""""""
function get_cut_node_sets_w_pop(
    edge::Tuple,
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelSubGraph,
    ext_decorators::Dict{Tuple,Dict{Tuple,Int}}=Dict{Tuple,Dict{Tuple,Int}}()
)
    n1, n2 = edge
    key = get_edge_key(n1, n2)
    level = length(key) + 1
    cuttable_tree = multiscale_cuttable_tree.cuttable_trees[level][key]
    n1_id = cuttable_tree.name_to_index[n1[1:level]]
    pop1 = cuttable_tree.pop_weight[n1_id]
    pop2 = get_total_population(subgraph, ext_decorators) - pop1

    node_set1 = Dict{Tuple{Vararg{String}}, Any}()
    node_set2 = Dict{Tuple{Vararg{String}}, Any}()
    node_set1[key] = Dict{Tuple{Vararg{String}}, Any}()
    node_set2[key] = Dict{Tuple{Vararg{String}}, Any}()

    n1_enter = cuttable_tree.enter[n1_id]
    n1_exit = cuttable_tree.exit[n1_id]

    for n_id = 1:nv(cuttable_tree.tree)
        n_enter = cuttable_tree.enter[n_id]
        n_exit = cuttable_tree.exit[n_id]
        g_n_id = cuttable_tree.vmap[n_id]
        node = subgraph.parent.id_to_partitions[level][g_n_id]
        if n_enter >= n1_enter && n_exit <= n1_exit
            node_set1[key][node] = get_node_set(subgraph, node)
        else
            node_set2[key][node] = get_node_set(subgraph, node)
        end
    end

    node_sets = get_coarsened_node_sets(multiscale_cuttable_tree, subgraph,
                                        edge, node_set1, node_set2, key)
    node_set1, node_set2 = node_sets

    if pop1 < pop2
        return [(node_set1, pop1), (node_set2, pop2)]
    else
        return [(node_set2, pop2), (node_set1, pop1)]
    end
end

""""""
function cut_edge(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelSubGraph,
    rng::AbstractRNG
)
    cuttable_edges = Set{Tuple}()
    get_all_cuttable_edges!(cuttable_edges, multiscale_cuttable_tree)
    edge, prob_edge = choose_random_cuttable_edge(cuttable_edges,
                                                  multiscale_cuttable_tree, rng)
    if edge === nothing
        return [], edge, nothing
    end
    node_sets_w_pops = get_cut_node_sets_w_pop(edge, multiscale_cuttable_tree,
                                               subgraph)
    return node_sets_w_pops, edge, prob_edge
end



""""""
function cut_edge2(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelSubGraph,
    constraints::Dict,
    rng::AbstractRNG
)
    cuttable_edges = Set{Set{Tuple}}()
    get_all_cuttable_edge_pairs!(cuttable_edges, multiscale_cuttable_tree,
                                 subgraph, constraints)
    edges, prob_edges = choose_random_cuttable_edge(cuttable_edges,
                                multiscale_cuttable_tree, rng)
    if edges === nothing
        return [], edges, 0
    end

########## put in function
    edges = [e for e in edges]
    edge1 = edges[1]
    edge2 = edges[2]
    node_sets_w_pops_1 = get_cut_node_sets_w_pop(edge1,multiscale_cuttable_tree,
                                                 subgraph)
    node_sets_w_pops_2 = get_cut_node_sets_w_pop(edge2,multiscale_cuttable_tree,
                                                 subgraph)
    district1 = node_sets_w_pops_1[1]
    district2 = node_sets_w_pops_2[1]
    complement_district1 = node_sets_w_pops_1[2]
    complement_district2 = node_sets_w_pops_2[2]
    district3_nodes = intersect_node_sets(complement_district1[1], 
                                          complement_district2[1])
    dist12_pop = district1[2] + district2[2]
    district3_pop = get_total_graph_population(subgraph) - dist12_pop
    node_sets_w_pops = [district1, district2, (district3_nodes, district3_pop)]
############

    return node_sets_w_pops, edges, prob_edges
end
