struct MCDSplits
    mcd_names::Vector{String}
    mcd_pops::Vector{Real}
    mcd_ids::Dict{String, Int}
    splits::Array{Real,2}
    ideal_dist_pop::Real
end

function build_mcd_score(
    partition::MultiLevelPartition;
    ideal_pop::Union{Nothing,Real}=nothing
)
    mcds = MCDSplits(partition, ideal_pop=ideal_pop)
    get_mcd_score!(mcds, partition)
    f(p, d=collect(1:p.num_dists)) = get_mcd_score!(mcds, p, d)
    return f
end

function get_mcd_score!(
    mcds::MCDSplits,
    partition::MultiLevelPartition,
    districts::Vector{Int} = collect(1:partition.num_dists)
)
    graph = partition.graph
    ousted_population = 0
    mcd_ids = Set{Int}()
    # @show districts
    for di in districts
        # @show keys(partition.district_to_nodes[di])
        mcds.splits[di,:] .= 0 
        union!(mcd_ids, get_mcd_splits!(mcds, partition, di))
        # @show "after union", mcd_ids
    end
    ideal_pop = mcds.ideal_dist_pop
    # @show ideal_pop
    sorted_perm = zeros(Int, partition.num_dists)
    for mcd_id in mcd_ids
        mcd_pop = mcds.mcd_pops[mcd_id]
        num_dists = mcd_pop/ideal_pop
        whole_dists = Int(floor(num_dists))
        # if mcd_id == 4 
        #     @show mcd_pop, num_dists
        # end
        sortperm!(sorted_perm, mcds.splits[:, mcd_id], rev=true)
        if num_dists > 1
            for ii = 1:whole_dists
                di = sorted_perm[ii]
                ousted_population += partition.dist_populations[di]
                ousted_population -= mcds.splits[di, mcd_id]
                # if mcd_id == 4 
                #     dpop = partition.dist_populations[di]
                #     split = mcds.splits[di, mcd_id]
                #     ousted = partition.dist_populations[di]-mcds.splits[di, mcd_id] 
                #     @show ii, mcd_pop, num_dists, dpop, split, ousted
                # end
            end 
            df = sorted_perm[whole_dists+1] 
            fraction = num_dists - floor(num_dists)
            ousted_frac = ideal_pop*fraction - mcds.splits[df, mcd_id]
            if ousted_frac > 0
                ousted_population += ousted_frac
            end
            # if mcd_id == 4 
            #     ideal_rem = ideal_pop*fraction
            #     @show ousted_frac, ideal_rem
            # end
        else 
            rep_dist = sorted_perm[1]
            ousted_population += mcd_pop - mcds.splits[rep_dist, mcd_id]
        end
    end
    # @show mcds.splits[:,4], ousted_population
    # @show ousted_population
    # @show mcds
    return ousted_population
end

function get_mcd_splits!(
    mcds::MCDSplits,
    partition::MultiLevelPartition,
    di::Int,
    node_set::Dict{Tuple{Vararg{String}}, Any}=partition.district_to_nodes[di],
    level::Int=1,
)::Set{Int}
    mcd_ids = Set{Int}()
    level_graph = partition.graph.graphs_by_level[level]
    node_attributes = level_graph.node_attributes
    mcd_col = level_graph.mcd_col
    for node in keys(node_set)
        # @show node
        if intact_huh(node_set[node])
            node_id = partition.graph.partition_to_ids[level][node]
            for (mcd, pop) in node_attributes[node_id][mcd_col]
                mcd_id = mcds.mcd_ids[mcd]
                mcds.splits[di, mcd_id] += pop
                push!(mcd_ids, mcd_id)
                # @show di, node, mcd_id, mcd, pop, mcds.mcd_pops[mcd_id], mcds.splits[di, mcd_id]
            end
        else
            union!(mcd_ids, get_mcd_splits!(mcds, partition, di, node_set[node], 
                                            level+1))
        end
    end
    return mcd_ids
end

function MCDSplits(
    partition::MultiLevelPartition; 
    ideal_pop::Union{Real, Nothing}=nothing
)::MCDSplits
    mcd_names = Vector{String}(undef, 0)
    mcd_pops = Vector{Float64}(undef, 0)
    mcd_ids = Dict{String, Int}()
    graph = partition.graph
    coarse_graph = graph.graphs_by_level[1]
    node_attributes = coarse_graph.node_attributes
    mcd_col = coarse_graph.mcd_col
    @assert typeof(mcd_col) == String

    mcd_id_count = 1
    for node_id = 1:coarse_graph.num_nodes
        mcds = node_attributes[node_id][mcd_col]
        for (mcd, pop) in mcds
            if haskey(mcd_ids, mcd)
                mcd_pops[mcd_ids[mcd]] += pop
            else
                push!(mcd_pops, pop)
                push!(mcd_names, mcd)
                mcd_ids[mcd] = mcd_id_count
                mcd_id_count += 1
            end
        end
    end
    splits = zeros(partition.num_dists, length(mcd_pops))
    if ideal_pop == nothing
        ideal_pop = coarse_graph.total_pop/partition.num_dists
    end
    return MCDSplits(mcd_names, mcd_pops, mcd_ids, splits, ideal_pop)
end