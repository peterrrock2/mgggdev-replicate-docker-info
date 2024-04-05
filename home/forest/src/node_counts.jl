function get_node_counts(
    partition::MultiLevelPartition,
    districts::Vector{Int} = collect(1:partition.num_dists)
)
    graph = partition.graph
    counts = Vector{Float64}()
    for di in districts
        count = get_node_count(graph, partition.district_to_nodes[di])
        push!(counts, count)
    end
    @assert sum(counts) == graph.graphs_by_level[end].num_nodes
    return counts
end


""""""
function get_node_count(
    graph::MultiLevelGraph,
    node_set::Dict{Tuple{Vararg{String}}, Any}, 
    level::Int=1
)
    count = 0
    full_graph = graph.graphs_by_level[level]
    
    for node in keys(node_set)
        if intact_huh(node_set[node])
            if graph.num_levels == level
                count += 1
            elseif graph.num_levels == level-1
                node_id = graph.partition_to_ids[level][node]
                count += nv(graph.coarse_to_fine_graphs[level][node_id].graph)
            else
                node_id = graph.partition_to_ids[level][node]
                vmap = graph.coarse_to_fine_graphs[level][node_id].vmap
                sub_node_set = Dict{Tuple{Vararg{String}}, Any}()
                for sub_node_id in vmap
                    sub_node = graph.id_to_partitions[level+1][sub_node_id]
                    sub_node_set[sub_node] = nothing
                end
                count += get_node_count(graph, sub_node_set, level+1)
            end
        else
            count += get_node_count(graph, node_set[node], level+1)
        end
    end
    return count
end