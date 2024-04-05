""""""
@inline function intact_huh(node_set)
    return node_set == nothing || length(node_set) == 0
end

""""""
@inline function intersects_huh(
    node_set::Dict{Tuple{Vararg{String}}, Any},
    node
)
    tmp_set = node_set
    cur_level = Set([length(k) for k in keys(node_set)])
    @assert length(cur_level) == 1
    cur_level = collect(cur_level)[1]

    for ii = cur_level:length(node)
        coarse_key = node[1:ii]
        if coarse_key ∉ keys(tmp_set)
            return false
        elseif intact_huh(tmp_set[coarse_key])
            return true
        end
        tmp_set = tmp_set[coarse_key]
    end
    return true
end

""""""
@inline function strictly_contained_huh(
    node_set::Dict{Tuple{Vararg{String}}, Any},
    node::Tuple{Vararg{String}}
)
    cur_level_node_set = length(collect(keys(node_set))[1]) ##current node_setrict level,
    node_level = length(node) ##the length of the node to be judged,

    if node_level < cur_level_node_set
        return false
    elseif node_level == cur_level_node_set
        if haskey(node_set, node) && intact_huh(node_set[node])
            return true
        else
            return false
        end
    else # node_level > cur_level_node_set
        coarse_node = node[1:cur_level_node_set]
        if haskey(node_set, coarse_node)
            if intact_huh(node_set[coarse_node])
                return true
            else
                return strictly_contained_huh(node_set[coarse_node], node)
            end
        else
            return false
        end
    end
end

""""""
@inline function strictly_contained_huh(
    node_set::Dict{Tuple{Vararg{String}}, Any},
    nodes::Set
)
    for node in nodes
        if !strictly_contained_huh(node_set, node)
            return false
        end
    end
    return true
end

""""""
@inline function strictly_contained_huh(
    super_node_set::Dict{Tuple{Vararg{String}}, Any},
    sub_node_set::Dict{Tuple{Vararg{String}}, Any}
)
    for node in keys(sub_node_set)
        if intact_huh(sub_node_set[node])
            if !strictly_contained_huh(super_node_set, node)
                return false
            end
        else
            for sub_node in keys(sub_node_set[node])
                if !strictly_contained_huh(super_node_set, sub_node)
                    return false
                end
            end
        end
    end
    return true
end

""""""
function intersect_node_sets(
    node_set₁::Union{Nothing, Dict{Tuple{Vararg{String}}, Any}},
    node_set₂::Union{Nothing, Dict{Tuple{Vararg{String}}, Any}},
)
    if node_set₁ == nothing
        return node_set₂
    elseif node_set₂ == nothing
        return node_set₁
    end

    intersecting_set = Dict{Tuple{Vararg{String}}, Any}()
    for node in keys(node_set₁)
        if node ∉ keys(node_set₂)
            continue
        end
        intersection = intersect_node_sets(node_set₁[node], node_set₂[node])
        if typeof(intersection) <: Dict && length(intersection) == 0
            continue
        end
        intersecting_set[node] = intersection
    end
    return intersecting_set
end

""""""
@inline function get_node_set(node_set::Dict{Tuple{Vararg{String}}, Any}, node)
    tmp_set = node_set
    for ii = 1:length(node)
        coarse_key = node[1:ii]
        if intact_huh(tmp_set[coarse_key])
            return nothing
        end
        tmp_set = tmp_set[coarse_key]
    end
    return tmp_set
end

""""""
@inline function get_intact_node_set(graph::MultiLevelGraph, node)
    level = length(node)
    if level == graph.num_levels
        return nothing
    end
    node_id = graph.partition_to_ids[level][node]
    fine_graph = graph.coarse_to_fine_graphs[level][node_id]
    node_set = Dict{Tuple{Vararg{String}},Any}()
    for fine_node_id in fine_graph.vmap
        fine_node = graph.id_to_partitions[level+1][fine_node_id]
        node_set[fine_node] = nothing
    end
    return node_set
end

""""""
function recursive_merge(x::AbstractDict...)
    return merge(recursive_merge, x...)
end

""""""
function correct_intact_nodes(graph, node_set)
    for (node, fine_nodes) in node_set
        if intact_huh(node_set[node])
            continue
        end
        node_set[node] = correct_intact_nodes(graph, fine_nodes)
        level = length(node)
        node_id = graph.partition_to_ids[level][node]
        whole_fine_nodes = sum([Int(intact_huh(node_set[node][fn]))
                                for fn in keys(fine_nodes)])
        num_fine_nodes = nv(graph.coarse_to_fine_graphs[level][node_id].graph)
        if whole_fine_nodes == num_fine_nodes
            node_set[node] = nothing
        end
    end
    return node_set
end

""""""
function merge_nodesets(graph, node_sets)
    node_set = recursive_merge(node_sets...)
    node_set = correct_intact_nodes(graph, node_set)
    return node_set
end

""""""
function get_node_set_complement(
    all_node_set::Dict{Tuple{Vararg{String}}, Any}, 
    sub_node_set::Dict{Tuple{Vararg{String}}, Any}, 
    graph::MultiLevelGraph)
    
    @assert strictly_contained_huh(all_node_set, sub_node_set)
    complement_node_set = Dict{Tuple{Vararg{String}}, Any}()
    for node in keys(all_node_set)
        if !haskey(sub_node_set, node)
            complement_node_set[node] = all_node_set[node]
        elseif !intact_huh(sub_node_set[node])
            if intact_huh(all_node_set[node])
                sub_all = get_intact_node_set(graph, node)
            else
                sub_all = all_node_set[node]
            end
            sub_sub = sub_node_set[node]
            sub_complement = get_node_set_complement(sub_all, sub_sub, graph)
            if length(sub_complement) > 0
                 complement_node_set[node] = sub_complement
            end

        end
    end
    return complement_node_set
end


# Cases to think about:
# (1)=>Dict((1,a), (1,b), (1,c))
# (1)=>Dict((1,a), (1,b), (1,c))
# complement_node_set[(1)] = Dict()


# (1)=>Dict((1,a), (1,b), (1,c))
# (1)=>Dict((1,a), (1,c))

# (1)=>nothing => (1)=>Dict((1,a), (1,b), (1,c), (1,d))
# (1)=>Dict((1,a), (1,b), (1,c))