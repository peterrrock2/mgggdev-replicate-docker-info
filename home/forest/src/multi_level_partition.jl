@enum EXTENSIONS changed_districts node_trees tracked_trees big_coarse_nodes tracked_step replica_id rejection_counter bath_swaps del_dists measure_id

# tracked_trees => Vect[Dicts]j
# index is district, dict node to spanning tree count
# on accept
# tracked_trees <- proposed_tracked_trees
# clear proposed_tracked_trees

abstract type AbstractPartition end

mutable struct MultiLevelPartition <: AbstractPartition
    num_dists::Int
    cross_district_edges::Array{Real}          # real in case weighted edges
    node_to_district::Dict{Tuple{Vararg{String}}, Int}
    district_to_nodes::Vector{Dict{Tuple{Vararg{String}}, Any}}
    subgraphs::Vector{MultiLevelSubGraph}
    dist_populations::Vector{Int}              # of length(num_districts)
    graph::MultiLevelGraph
    parent::Union{MultiLevelPartition,Nothing} # optional parent partition
    extensions::Dict{EXTENSIONS,Any}
    proposed_extensions::Dict{EXTENSIONS,Any}
end

function reset!(
    partition_dst::MultiLevelPartition, 
    partition_src::MultiLevelPartition;
    ext_huh=true
)
    partition_dst.num_dists = partition_src.num_dists
    partition_dst.cross_district_edges = partition_src.cross_district_edges
    partition_dst.node_to_district = partition_src.node_to_district
    partition_dst.district_to_nodes = partition_src.district_to_nodes
    partition_dst.subgraphs = partition_src.subgraphs
    partition_dst.dist_populations = partition_src.dist_populations
    partition_dst.graph = partition_src.graph
    partition_dst.parent = partition_src.parent
    if ext_huh
        partition_dst.extensions = partition_src.extensions
        partition_dst.proposed_extensions = partition_src.proposed_extensions
    end
end

#district_to_nodes[1]
#(Alamance,) => ()
#(Burlington,) => Dict((Burlington, 01) => (),
#                      (Burlington, 02) => Dict((Burlington, 02, 37012309123) => (), ...)
#                      )
#...

""""""
function construct_districts_node_map!(
    node_to_district::Dict{Tuple{Vararg{String}},Int},
    sub_nodes::Dict{Tuple{Vararg{String}},Any},
    district::Int
)
    for key in keys(sub_nodes)
        if intact_huh(sub_nodes[key])
            node_to_district[key] = district
        else
            sub_sub_nodes = sub_nodes[key]
            construct_districts_node_map!(node_to_district, sub_sub_nodes,
                                          district)
        end
    end
end

""""""
function construct_node_map(
    district_to_nodes::Vector{Dict{Tuple{Vararg{String}}, Any}}
)::Dict{Tuple{Vararg{String}},Int}
    node_to_district = Dict{Tuple{Vararg{String}},Int}()
    for district = 1:length(district_to_nodes)
        for key in keys(district_to_nodes[district])
            if intact_huh(district_to_nodes[district][key])
                node_to_district[key] = district
            else
                sub_nodes = district_to_nodes[district][key]
                construct_districts_node_map!(node_to_district, sub_nodes,
                                              district)
            end
        end
    end
    return node_to_district
end

""""""
function reset_cross_district_edges!(
    partition::MultiLevelPartition,
    changed_districts::Vector{Int} = collect(1:partition.num_dists)
)
    cross_district_edges = partition.cross_district_edges

    num_dists = partition.num_dists
    num_levels = partition.graph.num_levels
    cds = changed_districts

    to_reset = view(cross_district_edges, cds, 1:num_dists, 1:num_levels)
    fill!(to_reset, 0)
    to_reset = view(cross_district_edges, 1:num_dists, cds, 1:num_levels)
    fill!(to_reset, 0)
end

""""""
function reset_cross_district_edges!(
    partition::MultiLevelPartition,
    changed_districts::Tuple{Vararg{Int}}
)
    cross_district_edges = partition.cross_district_edges

    num_dists = partition.num_dists
    num_levels = partition.graph.num_levels
    cds = [cd for cd in changed_districts]

    to_reset = view(cross_district_edges, cds, 1:num_dists, 1:num_levels)
    fill!(to_reset, 0)
    to_reset = view(cross_district_edges, 1:num_dists, cds, 1:num_levels)
    fill!(to_reset, 0)
end

""""""
function set_cross_district_edges!(
    partition::MultiLevelPartition,
    node_set::Dict{Tuple{Vararg{String}}, Any},
    node_dist::Int,
    changed_districts::Union{Vector{Int}, Tuple{Vararg{Int}}},
    level::Int=1
)
    graph = partition.graph
    num_levels = graph.num_levels
    levels = graph.levels
    simple_graph = graph.graphs_by_level[level].simple_graph

    cross_district_edges = partition.cross_district_edges
    district_to_nodes = partition.district_to_nodes

    for node in keys(node_set)
        node_id = graph.partition_to_ids[level][node]

        for nbr_id in neighbors(simple_graph, node_id)
            nbr = graph.id_to_partitions[level][nbr_id]

            if nbr[1:level-1]!=node[1:level-1]
                continue
            end

            districts_in_nbr = get_districts(partition, nbr)

            for di in districts_in_nbr
                if di == node_dist
                    continue
                end
                if di in changed_districts && nbr_id > node_id
                    continue
                end
                nbr_set = get_node_set(district_to_nodes[di], nbr)
                edge_weight = get_edge_weight(graph, node_set[node], nbr_set,
                                              node, nbr)
                cross_district_edges[node_dist, di, level] += edge_weight
                cross_district_edges[di, node_dist, level] += edge_weight
                # if di == 1 || di == 4
                #     @show di, edge_weight, node, nbr, level
                # end
            end
        end
        if level < num_levels
            if !intact_huh(node_set[node])
                # println("recursing on node ", node)
                set_cross_district_edges!(partition, node_set[node], node_dist,
                                          changed_districts, level+1)
            end
        end
    end
end

""""""
function set_cross_district_edges!(
    partition::MultiLevelPartition,
    changed_districts::Vector{Int} = collect(1:partition.num_dists)
)
    reset_cross_district_edges!(partition, changed_districts)

    districts = partition.district_to_nodes
    for cd in changed_districts
        set_cross_district_edges!(partition, districts[cd], cd,
                                  changed_districts)
    end
end

""""""
function precompute_node_tree_counts!(
    partition::MultiLevelPartition
)
    if haskey(partition.extensions, node_trees::EXTENSIONS)
        return
    end

    node_tree_counts = Vector{Vector{Float64}}(undef, 0)
    graph = partition.graph
    for level = 1:graph.num_levels-1
        full_graph = graph.graphs_by_level[level]
        simple_graph = full_graph.simple_graph
        node_tree_counts_level = Vector{Float64}(undef, nv(simple_graph))
        for node_id = 1:nv(simple_graph)
            subgraph = graph.coarse_to_fine_graphs[level][node_id]
            sg = subgraph.graph
            if nv(sg) == 1
                count = 0
            else
                count = log_nspanning(sg)
            end
            node_tree_counts_level[node_id] = count
        end
        push!(node_tree_counts, node_tree_counts_level)
    end
    partition.extensions[node_trees::EXTENSIONS] = node_tree_counts
end

""""""
function add_node_to_district!(
    district_to_nodes::Dict{Tuple{Vararg{String}},Any},
    node::Tuple{Vararg{String}},
    level::Int=1
)
    if length(node) == level
        district_to_nodes[node] = nothing
    else
        coarse_node = node[1:level]
        if coarse_node ∉ keys(district_to_nodes)
            district_to_nodes[coarse_node] = Dict{Tuple{Vararg{String}},Any}()
        end
        add_node_to_district!(district_to_nodes[coarse_node], node, level+1)
    end
end

function remove_node_from_district!(
    district_to_nodes::Dict{Tuple{Vararg{String}},Any},
    graph::MultiLevelGraph,
    node::Tuple{Vararg{String}},
    level::Int=1
)
    if level == length(node)
        delete!(district_to_nodes, node)
    elseif haskey(district_to_nodes, node[1:level])
        if intact_huh(district_to_nodes[node[1:level]])
            district_to_nodes[node[1:level]] = get_intact_node_set(graph,
                                                                  node[1:level])
        end
        remove_node_from_district!(district_to_nodes[node[1:level]], graph,
                                   node, level+1)
        if length(district_to_nodes[node[1:level]]) == 0
            delete!(district_to_nodes, node[1:level])
        end
    end
end

""""""
function MultiLevelPartition(
    graph::MultiLevelGraph,
    node_to_district::Dict{Tuple{Vararg{String}}, Int};
    copy_huh::Bool=true
)::MultiLevelPartition
    if copy_huh
        node_to_district = deepcopy(node_to_district)
    end

    num_dists = maximum(values(node_to_district))
    district_to_nodes = Vector{Dict{Tuple{Vararg{String}},Any}}(undef,
                                                                num_dists)
    for di = 1:num_dists
        district_to_nodes[di] = Dict{Tuple{Vararg{String}},Any}()
    end

    for (node, dist) in node_to_district
        add_node_to_district!(district_to_nodes[dist], node)
    end
    for di in keys(district_to_nodes)
        district_to_nodes[di] = correct_intact_nodes(graph,
                                                     district_to_nodes[di])
    end

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

    # check compression in node_to_district
    derived_node_to_district = construct_node_map(district_to_nodes)
    if node_to_district != derived_node_to_district && !copy_huh
        println("Warning: node_to_district was not passed in a compressed "*
                "format but MultiLevelPartition requested memory transfer. "*
                "Ignoring request and using different memory with compressed "*
                "format.")
    elseif node_to_district != derived_node_to_district && !copy_huh
        derived_node_to_district = node_to_district
    end

    partition = MultiLevelPartition(num_dists, cross_district_edges,
                                    derived_node_to_district, district_to_nodes,
                                    subgraphs, dist_populations, graph, parent,
                                    extensions, proposed_extensions)

    set_cross_district_edges!(partition)

    return partition
end

""""""
function MultiLevelPartition(
    multi_level_graph::MultiLevelGraph,
    assignment_col::AbstractString
)::MultiLevelPartition
    node_to_district = Dict{Tuple{Vararg{String}}, Int}()
    fine_level = multi_level_graph.num_levels
    fine_graph = multi_level_graph.graphs_by_level[fine_level]
    node_attributes = fine_graph.node_attributes
    for ii = 1:length(node_attributes)
        node = multi_level_graph.id_to_partitions[fine_level][ii]
        district = convert(Int,node_attributes[ii][assignment_col])
        node_to_district[node] = district
    end
    min_district = minimum(values(node_to_district))
    if min_district == 0
        for node in keys(node_to_district)
           node_to_district[node] += 1
        end
    end
    return MultiLevelPartition(multi_level_graph, node_to_district)
end

""""""
function get_district(
    partition::MultiLevelPartition,
    node::Tuple{Vararg{String}}
)
    for ii = 1:length(node)
        if haskey(partition.node_to_district, node[1:ii])
            return partition.node_to_district[node[1:ii]]
        end
    end
    return nothing
end

function get_districts(
    partition::MultiLevelPartition,
    node::Tuple{Vararg{String}}
)
    d = get_district(partition, node)
    if d != nothing
        return [d]
    end
    num_dists = partition.num_dists
    district_to_nodes = partition.district_to_nodes
    return [di for di = 1:num_dists
            if intersects_huh(district_to_nodes[di], node)]
end


""""""
# function get_districts_in_node(
#     partition::MultiLevelPartition,
#     node::Tuple{Vararg{String}}
# )
#     for level = 1:level(node)

# end
