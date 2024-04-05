function get_isoperimetric_score(
    partition::MultiLevelPartition,
    districts::Vector{Int} = collect(1:partition.num_dists)
)
    graph = partition.graph
    score = 0
    for di in districts
        area = get_area(graph, partition.district_to_nodes[di])
        pr = get_perimeter(graph, partition.district_to_nodes[di])
        isoperimetric_score = pr^2/area
        score += isoperimetric_score
    end
    return score
end

function get_isoperimetric_scores(
    partition::MultiLevelPartition,
    districts::Vector{Int} = collect(1:partition.num_dists)
)
    graph = partition.graph
    scores = Vector{Float64}()
    for di in districts
        area = get_area(graph, partition.district_to_nodes[di])
        pr = get_perimeter(graph, partition.district_to_nodes[di])
        isoperimetric_score = pr^2/area
        push!(scores, isoperimetric_score)
    end
    return scores
end


function get_polsby_popper_scores(
    partition::MultiLevelPartition,
    districts::Vector{Int} = collect(1:partition.num_dists)
)
    graph = partition.graph
    scores = Vector{Float64}()
    for di in districts
        area = get_area(graph, partition.district_to_nodes[di])
        pr = get_perimeter(graph, partition.district_to_nodes[di])
        polsby_popper_score = (4*pi*area)/(pr^2)
        push!(scores, polsby_popper_score)
    end
    return scores
end


""""""
function get_perimeters(
    partition::MultiLevelPartition,
    districts::Vector{Int} = collect(1:partition.num_dists)
)
    graph = partition.graph
    perimeters = Vector{Float64}()
    for di in districts
        perimeter = get_perimeter(graph, partition.district_to_nodes[di])
        push!(perimeters, perimeter)
    end
    return perimeters
end


""""""
function get_area(
    graph::MultiLevelGraph,
    district::Dict{Tuple{Vararg{String}}, Any}, 
    node_set::Dict{Tuple{Vararg{String}}, Any}=district, 
    level::Int=1
)
    area = 0
    full_graph = graph.graphs_by_level[level]
    area_col = full_graph.area_col

    for node in keys(node_set)
        simple_graph = graph.graphs_by_level[level].simple_graph
        node_id = graph.partition_to_ids[level][node] 
        if !intact_huh(node_set[node])
            area += get_area(graph, district, node_set[node], level+1)
        else
            area += full_graph.node_attributes[node_id][area_col]
        end
    end
               
    return area
end


""""""
function get_perimeter(
    graph::MultiLevelGraph, 
    district::Dict{Tuple{Vararg{String}}, Any}, 
    node_set::Dict{Tuple{Vararg{String}}, Any}=district, 
    level::Int=1
)
    pr = 0
    full_graph = graph.graphs_by_level[level]
    edge_attributes = full_graph.edge_attributes
    node_border_col = full_graph.node_border_col
    edge_perimeter_col = full_graph.edge_perimeter_col

    for node in keys(node_set) 
        simple_graph = graph.graphs_by_level[level].simple_graph 
        node_id = graph.partition_to_ids[level][node] 
        if !intact_huh(node_set[node])
            pr += get_perimeter(graph, district, node_set[node], level+1) 
        else
            pr += full_graph.node_attributes[node_id][node_border_col]
            for nbr_id in neighbors(simple_graph, node_id) 
                nbr = graph.id_to_partitions[level][nbr_id] 
                if !strictly_contained_huh(district, nbr)
                    if intersects_huh(district,nbr) 
                        fine_neighbors = graph.fine_neighbors[level]
                        fine_neighbors = fine_neighbors[node_id][nbr_id]
                        pr += get_intersection_perimeter(graph, district, 
                                                         fine_neighbors, node,
                                                         level)
                    else
                        edge = Set([node_id, nbr_id])
                        pr += edge_attributes[edge][edge_perimeter_col]  
                    end
                end
            end
        end
    end       
    return pr
end


""""""
function get_intersection_perimeter(
    graph::MultiLevelGraph, 
    district::Dict{Tuple{Vararg{String}}, Any}, 
    neighbor_set::Dict{Int64,Any},
    node::Tuple{Vararg{String}},
    level::Int
)
    pr = 0
    for (nbr_id, finer_nbrs) in neighbor_set
        nbr = graph.id_to_partitions[level+1][nbr_id]
        if strictly_contained_huh(district, nbr)
            continue
        elseif intersects_huh(district, nbr)
            pr += get_intersection_perimeter(graph, district, finer_nbrs, node,
                                             level+1)
        else
            pr += compute_intersection_perimeter(graph, district, nbr_id, node,
                                                 level+1)
        end
    end
    return pr
end  


function compute_intersection_perimeter(
    graph::MultiLevelGraph, 
    district::Dict{Tuple{Vararg{String}}, Any}, 
    nbr_id::Int64,
    node::Tuple{Vararg{String}},
    level::Int
)
    pr = 0
    simple_graph = graph.graphs_by_level[level].simple_graph 
    edge_perimeter_col = graph.graphs_by_level[level].edge_perimeter_col
    edge_attributes = graph.graphs_by_level[level].edge_attributes

    for nbrr_id in neighbors(simple_graph, nbr_id) 
        nbrr = graph.id_to_partitions[level][nbrr_id]
        if nbrr[1:length(node)] == node
            edge = Set([nbr_id, nbrr_id])
            pr += edge_attributes[edge][edge_perimeter_col]
        end
    end
    return pr
end  
