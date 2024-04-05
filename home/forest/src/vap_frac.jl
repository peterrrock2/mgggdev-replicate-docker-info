function get_minofracs(
    partition::MultiLevelPartition,
    mino_pop_col::String,
    pop_col::String
)
    popfracs = Vector{Real}(undef, partition.num_dists)

    for di = 1:partition.num_dists
        node_set = partition.subgraphs[di].node_set
        if pop_col == partition.graph.graphs_by_level[1].pop_col
            pop = partition.dist_populations[di]
        else 
            pop = sum_attribute(partition.graph, node_set, pop_col)
        end
        minopop = sum_attribute(partition.graph, node_set, mino_pop_col)
        popfracs[di] = minopop/pop
    end
    return popfracs
end

function get_vra_score(
    partition::MultiLevelPartition,
    target_mino_dists::Int,
    mino_pop_col::String,
    pop_col::String
)
    score = 0
    popfracs = get_minofracs(partition, mino_pop_col, pop_col)
    sort!(popfracs, rev=true)

    curVraDists = findfirst(x->x<0.5, popfracs)-1

    # @show popfracs
    for ii = 1:target_mino_dists
        if popfracs[ii] < 0.5
            score += sqrt(0.51-popfracs[ii])
        end
    end
    # for ii = 1:min(target_mino_dists, curVraDists+1)
    #     score += (0.55-popfracs[ii])^2
    # end
    return score
end

function build_get_vra_score(
    target_mino_dists::Int,
    mino_pop_col::String,
    pop_col::String
)
    f(p, d=collect(1:p.num_dists)) = get_vra_score(p, target_mino_dists, 
                                                   mino_pop_col, pop_col)
end 


