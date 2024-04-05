import Pkg
push!(LOAD_PATH, "..")

using RandomNumbers
using MultiScaleMapSampler

function run_multiscale2(
    ; pctGraphPath::String,
    subregion_name::String,
    region_name::String,
    population_col::String,
    output_path::Union{Some{String}, Nothing}=nothing,
    num_dists=1,
    rng_seed=42,
    pop_dev=0.2,
    gamma=0,
    steps=1000,
    edge_weights="connections",
)

    nodeData = Set([subregion_name, region_name, population_col])
    base_graph = BaseGraph(
        pctGraphPath,
        population_col,
        inc_node_data=nodeData,
        edge_weights=edge_weights
    )
    graph = MultiLevelGraph(base_graph, [region_name, subregion_name])

    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(graph, num_dists, 0.01))
    add_constraint!(constraints, ConstrainDiscontinuousTraversals(graph))
    add_constraint!(constraints, MaxCoarseNodeSplits(num_dists+1))

    rng = PCG.PCGStateOneseq(UInt64, rng_seed)
    partition = MultiLevelPartition(graph, constraints, num_dists; rng=rng)

    proposal = build_forest_recom2(constraints)
    measure = Measure(gamma)

    if output_path == nothing
        writer = Writer(measure, constraints, partition, stdout)
    else
        output_file_io = smartOpen(output_path.value, "w")
        writer = Writer(measure, constraints, partition, output_file_io)
    end

    run_metropolis_hastings!(
        partition,
        proposal,
        measure,
        steps,
        rng,
        writer=writer,
        output_freq=1
    )
    return nothing
end