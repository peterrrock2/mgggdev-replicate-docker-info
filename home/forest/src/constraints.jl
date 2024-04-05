# """
#     ContiguityConstraint()

# Initializes and returns a `ContiguityConstraint` object.
# """
# struct ContiguityConstraint <: AbstractConstraint
#     # No metadata (for now); implements the `AbstractConstraint` interface.
# end

""""""
function PopulationConstraint(
    graph::BaseGraph,
    num_dists::Int,
    tolerance::Float64,
    ideal_pop::Real = 0,
)::PopulationConstraint
    if ideal_pop == 0
        ideal_pop = graph.total_pop / num_dists
    end
    # no particular reason to not use floor() instead of ceil()
    min_pop = Int(ceil((1 - tolerance) * ideal_pop))
    max_pop = Int(floor((1 + tolerance) * ideal_pop))
    return PopulationConstraint(min_pop, max_pop)
end

""""""
function PopulationConstraint(
    graph::MultiLevelGraph,
    num_dists::Int,
    tolerance::Float64,
    ideal_pop::Real = 0,
)::PopulationConstraint
    graph = graph.graphs_by_level[1]
    return PopulationConstraint(graph, num_dists, tolerance, ideal_pop)
end

""""""
function PackNodeConstraint(
    graph::MultiLevelGraph,
    unpack::Int=0;
    num_dists::Int=0,
    ideal_pop::Real=0
)::PackNodeConstraint
    if ideal_pop == 0 && num_dists == 0
        throw(
            ArgumentError("Need to specify either ideal_pop or num_dists",)
        )
    elseif ideal_pop == 0
        ideal_pop = graph.graphs_by_level[1].total_pop / num_dists
    end
    packed_nodes = Dict{Tuple{Vararg{String}}, Int}()
    pop_col = graph.graphs_by_level[1].pop_col
    for level = 1:graph.num_levels
        g = graph.graphs_by_level[level]
        for node_id = 1:g.num_nodes
            node_pop = g.node_attributes[node_id][pop_col]
            packed_districts = floor(node_pop/ideal_pop) - unpack
            if packed_districts > 0
                node = graph.id_to_partitions[level][node_id]
                packed_nodes[node] = packed_districts
            end
        end
    end
    return PackNodeConstraint(packed_nodes, ideal_pop)
end

""""""
function ConstrainDiscontinuousTraversals(
    graph::MultiLevelGraph, max_line_segments::Int=1
)::ConstrainDiscontinuousTraversals
    return ConstrainDiscontinuousTraversals(max_line_segments)
end

""""""
function AllowedExcessDistsInCoarseNodes(
    graph::MultiLevelGraph,
    num_dists::Int,
    allowable_excess::Int=0;
    ideal_pop::Real = 0
)::AllowedExcessDistsInCoarseNodes
    if ideal_pop == 0
        ideal_pop = graph.graphs_by_level[1].total_pop / num_dists
    end
    return AllowedExcessDistsInCoarseNodes(allowable_excess, ideal_pop)
end

""""""
function MaxHammingDistance(
    graph::MultiLevelGraph,
    partition::MultiLevelPartition, #initial partition
    max_distance::Real, #maximum fraction of deviation in a single distrct
    norm_type::String = "infinity"
)::MaxHammingDistance
    return MaxHammingDistance(partition,max_distance,norm_type)
end

""""""
function add_constraint!(
    constraints::Dict,
    constraint::T
) where T <: AbstractConstraint
    constraints[T] = constraint
end

""""""
function initialize_constraints()
    return Dict{Type{T} where T<:AbstractConstraint, AbstractConstraint}()
end

""""""
function satisfies_constraint(
    constraint::PopulationConstraint,
    partition::MultiLevelPartition,
    districts::Vector{Int} = collect(1:partition.num_dists)
)
    for di in districts
        pop = partition.dist_populations[di]
        if pop < constraint.min_pop || pop > constraint.max_pop
            return false
        end
    end
    return true
end


""""""
function satisfies_constraint(
    constraint::ConstrainDiscontinuousTraversals,
    partition::MultiLevelPartition,
    subgraph::MultiLevelSubGraph,
    node_set::Union{Nothing, Dict{Tuple{Vararg{String}},Any}}=nothing
)
    @assert constraint.max_line_segments == 1 # not built yet otherwise

    if node_set == nothing
        node_set = subgraph.node_set
    end

    level = length(collect(keys(node_set))[1])+1

    for node in keys(node_set)
        if intact_huh(node_set[node])
            continue
        end
        if !is_connected_bf(subgraph.graphs_by_level[level][node])
            return false
        end
        if !satisfies_constraint(constraint, partition, subgraph,
                                 node_set[node])
            return false
        end
    end
    return true
end


""""""
function satisfies_constraint(
    constraint::PackNodeConstraint,
    graph::MultiLevelGraph,
    district_to_nodes::AbstractArray{Dict{Tuple{Vararg{String}},Any}, 1},
    num_dists::Int
)
    pop_col = graph.graphs_by_level[1].pop_col
    for (node, target_packed_dists) in constraint.nodes
        # println("loop head: node ", node, target_packed_dists)
        level = length(node)
        if node ∉ keys(graph.partition_to_ids[level])
            continue
        end
        node_id = graph.partition_to_ids[level][node]
        packed_dists_in_node = 0
        for di = 1:length(district_to_nodes)
            node_set = district_to_nodes[di]
            for ii = 1:level
                if length(node_set) > 1 || node[1:ii] ∉ keys(node_set)
                    break
                elseif ii == level
                    packed_dists_in_node += 1
                end
                if node_set[node[1:ii]] == nothing
                    packed_dists_in_node += 1
                    break
                end
                node_set = node_set[node[1:ii]]
            end
        end
        if packed_dists_in_node >= target_packed_dists
            continue
        end
        if num_dists == length(district_to_nodes)
            return false
        else
            claimed_pop = 0
            for di = 1:length(district_to_nodes)
                if !intersects_huh(district_to_nodes[di], node)
                    continue
                end
                node_set = get_node_set(district_to_nodes[di], node)
                if node_set == nothing
                    node_set = Dict(node=>nothing)
                end
                lvl_ns = length(collect(keys(node_set))[1])
                claimed_pop += sum_node_data(graph, node_set, pop_col, lvl_ns)
            end
            n_atr = graph.graphs_by_level[level].node_attributes
            total_pop = n_atr[node_id][pop_col]
            remaining_pop = total_pop-claimed_pop
            packed_dists_in_node += floor(remaining_pop/constraint.ideal_pop)
            if packed_dists_in_node < target_packed_dists
                return false
            end
        end
    end
    return true
end


""""""
function satisfies_constraint(
    constraint::MaxCoarseNodeSplits,
    graph::MultiLevelGraph,
    district_to_nodes::AbstractArray{Dict{Tuple{Vararg{String}},Any}, 1}
)
    coarse_graph = graph.graphs_by_level[1]
    coarse_node_splits = -coarse_graph.num_nodes
    for node_set in district_to_nodes
        coarse_node_splits += length(node_set)
    end
    return coarse_node_splits <= constraint.max_coarse_node_splits
end


""""""
function satisfies_constraint(
    constraint::MaxSharedNodes,
    partition::MultiLevelPartition
)
    levels = partition.graph.num_levels
    for di = 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            cross_district_edges = partition.cross_district_edges
            inner_coarse_edges = sum(cross_district_edges[di, dj, 2:levels])
            if inner_coarse_edges == 0
                continue
            end

            println("satisfies_constraint for MaxSharedNodes is currently broken... exiting.")
            @assert false
            # get_max_shared()

            # for node in keys(partition.district_to_nodes[di])
            #     if haskey(partition.district_to_nodes[dj], node)
            #         num_shared += 1


            di_coarse_nodes = Set(keys(partition.district_to_nodes[di]))
            dj_coarse_nodes = Set(keys(partition.district_to_nodes[dj]))
            shared_coarse_nodes = intersect(di_coarse_nodes, dj_coarse_nodes)
            num_shared = length(shared_coarse_nodes)
            if num_shared > constraint.max_shared_coarse_nodes
                return false
            end
        end
    end
    return true
end

""""""
function satisfies_constraint(
    constraint::MaxSharedCoarseNodes,
    partition::MultiLevelPartition
)
    levels = partition.graph.num_levels
    for di = 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            cross_district_edges = partition.cross_district_edges
            inner_coarse_edges = sum(cross_district_edges[di, dj, 2:levels])
            if inner_coarse_edges == 0
                continue
            end
            di_coarse_nodes = Set(keys(partition.district_to_nodes[di]))
            dj_coarse_nodes = Set(keys(partition.district_to_nodes[dj]))
            shared_coarse_nodes = intersect(di_coarse_nodes, dj_coarse_nodes)
            num_shared = length(shared_coarse_nodes)
            if num_shared > constraint.max_shared_coarse_nodes
                return false
            end
        end
    end
    return true
end


""""""
function satisfies_constraint(
    constraint::AllowedExcessDistsInCoarseNodes,
    graph::MultiLevelGraph,
    district_to_nodes::AbstractArray{Dict{Tuple{Vararg{String}},Any}, 1},
    num_dists::Int
)
    coarse_graph = graph.graphs_by_level[1]
    districts_in_node = zeros(coarse_graph.num_nodes)
    for node_set in district_to_nodes
        for node in keys(node_set)
            node_id = graph.partition_to_ids[1][node]
            districts_in_node[node_id] += 1
        end
    end
    # @show districts_in_node

    pop_col = coarse_graph.pop_col
    node_attributes = coarse_graph.node_attributes
    ideal_pop = constraint.ideal_pop

    if num_dists > length(district_to_nodes)
        claimed_pop = zeros(coarse_graph.num_nodes)
        for node_set in district_to_nodes
            for (node, sub_node_set) in node_set
                node_id = graph.partition_to_ids[1][node]
                if sub_node_set == nothing
                    sns = Dict{Tuple{Vararg{String}}, Any}(node => nothing)
                else
                    sns = sub_node_set
                end
                lvl_ns = length(collect(keys(sns))[1])
                node_pop_in_dist = sum_node_data(graph, sns, pop_col, lvl_ns)
                claimed_pop[node_id] += node_pop_in_dist
            end
        end

        for node_id = 1:coarse_graph.num_nodes
            pop_tot = node_attributes[node_id][pop_col]
            unclaimed_pop = pop_tot - claimed_pop[node_id]
            districts_in_node[node_id] += ceil(unclaimed_pop/ideal_pop)
        end
    end

    excess = constraint.excess_splitting
    node_pop = 0
    for node_id = 1:coarse_graph.num_nodes
        node_pop = node_attributes[node_id][pop_col]
        dists = districts_in_node[node_id]
        if node_pop < ideal_pop
            if dists > excess + 2
                # node = graph.id_to_partitions[1][node_id]
                # @show "here1", node, node_pop, ideal_pop, dists, excess+2
                return false
            end
        else
            needed_dists = ceil(node_pop/ideal_pop)
            if dists > needed_dists + excess
                # node = graph.id_to_partitions[1][node_id]
                # @show "here2", node, node_pop, ideal_pop, dists, excess, needed_dists
                return false
            end
        end
    end

    return true
end

""""""
function satisfies_constraint(
    constraint::MaxHammingDistance,
    graph::MultiLevelGraph,
    district_to_nodes::AbstractArray{Dict{Tuple{Vararg{String}},Any}, 1}
)
    init_partition = constraint.partition
    init_node_sets = init_partition.district_to_nodes
    prop_node_sets = district_to_nodes
    init_num_dists = length(init_node_sets)
    prop_num_dists = length(prop_node_sets)
    pop_col = init_partition.graph.graphs_by_level[1].pop_col
    pop_intersections = Array{Float64,2}(undef,init_num_dists,prop_num_dists)
    for i = 1:init_num_dists
        for j = 1:prop_num_dists
            intersection = intersect_node_sets(init_node_sets[i], prop_node_sets[j])
            pop_intersections[i,j] = sum_node_data(init_partition.graph, intersection, pop_col)
        end
    end
    neg_pop_intersections = pop_intersections.*-1 #flip sign since hungarian() minimizes
    assignment = hungarian(neg_pop_intersections)
    matched_districts = [dist for dist = 1:init_num_dists if assignment[1][dist] != 0] #indices of original districts that are matched to proposed districts
    fraction_moved = Vector{Float64}()
    for dist in matched_districts
        fraction_moved_out = 1 - pop_intersections[dist,assignment[1][dist]]/sum(pop_intersections[dist,:])
        fraction_moved_in = 1 - pop_intersections[dist,assignment[1][dist]]/sum(pop_intersections[:,assignment[1][dist]])
        push!(fraction_moved,(fraction_moved_out+fraction_moved_in)/2)
    end
    if constraint.norm_type == "1" #l_1 norm
        return sum(fraction_moved)/length(matched_districts) <= constraint.max_distance
    else #l_infinity norm
        return maximum(fraction_moved) <= constraint.max_distance
    end
end


""""""
function satisfies_constraints(
    partition::MultiLevelPartition,
    constraints::Dict,
    proposed_cut::Union{Tuple, Nothing}=nothing
)
    if proposed_cut != nothing
        changed_districts, node_sets_w_pops, edge = proposed_cut
        old_dists = [partition.district_to_nodes[cd] 
                     for cd in changed_districts]
        old_pops = [partition.dist_populations[cd] for cd in changed_districts]
        for (ii,cd) in enumerate(changed_districts)
            partition.district_to_nodes[cd] = node_sets_w_pops[ii][1][()]
            partition.dist_populations[cd] = node_sets_w_pops[ii][2]
        end
        partition.node_to_district = construct_node_map(partition.district_to_nodes)
        set_cross_district_edges!(partition, changed_districts)
##########
        # tmp = deepcopy(partition.cross_district_edges)
        # set_cross_district_edges!(partition)
        # @assert tmp == partition.cross_district_edges
##########
    end

    satisfies_constraints = true
    graph = partition.graph
    district_to_nodes = partition.district_to_nodes
    num_dists = partition.num_dists

    if satisfies_constraints && haskey(constraints, PopulationConstraint)
        constraint = constraints[PopulationConstraint]
        if proposed_cut == nothing
            s = satisfies_constraint(constraint, partition)
        else
            s = satisfies_constraint(constraint, partition, changed_districts)
        end
        satisfies_constraints = satisfies_constraints && s
    end

    if satisfies_constraints && haskey(constraints, PackNodeConstraint)
        constraint = constraints[PackNodeConstraint]
        s = satisfies_constraint(constraint, graph, district_to_nodes,
                                 num_dists)
        satisfies_constraints = satisfies_constraints && s
    end

    if satisfies_constraints && haskey(constraints, MaxCoarseNodeSplits)
        constraint = constraints[MaxCoarseNodeSplits]
        s = satisfies_constraint(constraint, graph, district_to_nodes)
        satisfies_constraints = satisfies_constraints && s
    end

    if satisfies_constraints &&
       haskey(constraints, MaxSharedCoarseNodes)
        constraint = constraints[MaxSharedCoarseNodes]
        s = satisfies_constraint(constraint, partition)
        satisfies_constraints = satisfies_constraints && s
    end

    if satisfies_constraints && haskey(constraints, AllowedExcessDistsInCoarseNodes)
        constraint = constraints[AllowedExcessDistsInCoarseNodes]
        s = satisfies_constraint(constraint, graph, district_to_nodes,
                                 num_dists)
        satisfies_constraints = satisfies_constraints && s
    end

    if satisfies_constraints &&
       haskey(constraints, MaxHammingDistance)
        constraint = constraints[MaxHammingDistance]
        s = satisfies_constraint(constraint, graph, district_to_nodes)
        satisfies_constraints = satisfies_constraints && s
    end

    if proposed_cut != nothing
        for (ii,cd) in enumerate(changed_districts)
            partition.district_to_nodes[cd] = old_dists[ii]
            partition.dist_populations[cd] = old_pops[ii]
        end
        partition.node_to_district = construct_node_map(partition.district_to_nodes)
        set_cross_district_edges!(partition, changed_districts)
    end

    return satisfies_constraints
end

function satisfies_constraints(
    graph::MultiLevelGraph,
    district_to_nodes::AbstractArray{Dict{Tuple{Vararg{String}},Any}, 1},
    num_dists::Int,
    constraints::Dict
)
    if haskey(constraints, PackNodeConstraint)
        constraint = constraints[PackNodeConstraint]
        if !satisfies_constraint(constraint, graph, district_to_nodes,
                                 num_dists)
            return false
        end
    end

    if haskey(constraints, MaxCoarseNodeSplits)
        constraint = constraints[MaxCoarseNodeSplits]
        if !satisfies_constraint(constraint, graph, district_to_nodes)
            return false
        end
    end

    # skip for now -- only called in the constructor
    # if haskey(constraints, MaxSharedCoarseNodes)
    #     constraint = constraints[MaxSharedCoarseNodes]
    #     if !satisfies_constraint(constraint, graph, district_to_nodes)
    #         return false
    #     end
    # end

    if haskey(constraints, AllowedExcessDistsInCoarseNodes)
        constraint = constraints[AllowedExcessDistsInCoarseNodes]
        if !satisfies_constraint(constraint, graph, district_to_nodes,
                                 num_dists)
            return false
        end
    end

    if MaxHammingDistance in keys(constraints)
        constraint = constraints[MaxHammingDistance]
        if !satisfies_constraint(constraint, graph, district_to_nodes,
                                 num_dists)
            return false
        end
    end

    return true
end
