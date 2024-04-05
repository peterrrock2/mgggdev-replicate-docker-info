""""""
function run_metropolis_hastings!(
    partition::MultiLevelPartition,
    proposal::Union{Function,Vector{Tuple{T, Function}}},
    measure::Measure,
    steps::Union{Int,Tuple{Int,Int}},
    rng::AbstractRNG;
    writer::Union{Writer, Nothing}=nothing,
    output_freq::Int=250
) where T <: Real
    precompute_node_tree_counts!(partition)
    check_proposals_weights(proposal)
    
    initial_step, final_step = set_step_bounds(steps)
    if initial_step == 0 || initial_step == 1
        output(partition, measure, initial_step, 0, writer)
    end

    for step = initial_step:final_step
        proposal!, proposal_index = get_random_proposal(proposal, rng)
        p, update = proposal!(partition, measure, rng)
        if p == 0
            if mod(step, output_freq) == 0 && step != initial_step
                output(partition, measure, step, 0, writer)
            end
            continue
        end
        p *= get_delta_energy(partition, measure, update)

        if rand(rng) < p
            update_partition!(partition, update)
            if haskey(partition.extensions, del_dists::EXTENSIONS) 
                for cd in update[1]
                    partition.extensions[del_dists::EXTENSIONS][cd] += 1
                end
            end
        end
        if mod(step, output_freq) == 0 && step != initial_step
            output(partition, measure, step, 0, writer)
        end
    end
end


""""""
function output(
    partition::MultiLevelPartition,
    measure::Measure,
    step::Integer,
    count::Int,
    writer::Union{Writer, Nothing}
)
    if writer == nothing
        return
    end

    for (desc, f) in writer.map_output_data
        writer.map_param[desc] = f(partition)
    end
    if haskey(partition.extensions, replica_id::EXTENSIONS)
        writer.map_param["replica_id"] = partition.extensions[replica_id::EXTENSIONS]
        writer.map_param["bath_swaps"] = partition.extensions[bath_swaps::EXTENSIONS]
        writer.map_param["del_dists"] = partition.extensions[del_dists::EXTENSIONS]
    elseif haskey(partition.extensions, measure_id::EXTENSIONS)
        writer.map_param["measure_id"] = partition.extensions[measure_id::EXTENSIONS]
        writer.map_param["gamma"] = measure.gamma
        writer.map_param["energies"] = measure.descriptions
        writer.map_param["energy weights"] = measure.weights
    end
    if haskey(partition.extensions, rejection_counter::EXTENSIONS)
        for (k, v) in partition.extensions[rejection_counter::EXTENSIONS]
            writer.map_param[k] = v
        end
    end
##########
    if !writer.output_districting
        d = Dict{Tuple{Vararg{String}}, Int}()
        map = Map{MapParam}("step"*string(step-count), d, 1, writer.map_param)
    else
        map = Map{MapParam}("step"*string(step-count), 
                            partition.node_to_district, 1, writer.map_param)
    end
##########

    try
        addMap(writer.atlas.io, map)
    catch e
        @show writer.map_param
        @show count
        # println("partition.node_to_district: ", partition.node_to_district)
        println("Could not add map to atlas")
        @assert false
    end
    if haskey(partition.extensions, rejection_counter::EXTENSIONS)
        partition.extensions[rejection_counter::EXTENSIONS]["acceptance_wait"] = 0
    end
end


""""""
function update_partition!(
    partition::MultiLevelPartition,
    update::Tuple
)
    # if length(update) == 1 && typeof(update[1])==MultiLevelPartition
    #     partition.district_to_nodes = update[1].district_to_nodes
    #     partition.node_to_district = update[1].node_to_district
    #     partition.dist_populations = update[1].dist_populations
    #     partition.cross_district_edges = update[1].cross_district_edges
    # elseif length(update) == 3
    changed_districts, node_sets_w_pops, edge = update
    for (ii,cd) in enumerate(changed_districts)
        partition.district_to_nodes[cd] = node_sets_w_pops[ii][1][()]
        partition.subgraphs[cd] = MultiLevelSubGraph(partition.graph,
                                                    node_sets_w_pops[ii][1][()])
        partition.dist_populations[cd] = node_sets_w_pops[ii][2]
    end
    partition.node_to_district =
                             construct_node_map(partition.district_to_nodes)
    set_cross_district_edges!(partition, changed_districts)
    # end

    if rejection_counter::EXTENSIONS in keys(partition.extensions)
        for key in keys(partition.extensions[rejection_counter::EXTENSIONS])
            partition.extensions[rejection_counter::EXTENSIONS][key] = 0
        end
    end
end


""""""
@inline function set_step_bounds(steps)
    if typeof(steps)<:Tuple
        return steps
    else
        return 1, steps
    end
end


""""""
@inline function get_random_proposal(
    proposal::Union{Function, Vector{Tuple{T, Function}}},
    rng::AbstractRNG
) where T <: Real
    if !(typeof(proposal) <: Vector)
        return proposal, 1
    end
    proposal_weights = [proposal[i][1] for i = 1:length(proposal)]
    index = findfirst(cumsum(proposal_weights) .> rand(rng))
    return proposal[index][2], index
end


""""""
@inline function check_proposals_weights(
    proposal::Union{Function, Vector{Tuple{T, Function}}}
) where T<:Real
    if !(typeof(proposal) <: Vector)
        return
    end
    weight_sum = sum(proposal[i][1] for i = 1:length(proposal))
    if weight_sum != 1
        throw(
            ArgumentError(
                "Chance of choosing a proposal must sum to 1",
            ),
        )
    end
end
