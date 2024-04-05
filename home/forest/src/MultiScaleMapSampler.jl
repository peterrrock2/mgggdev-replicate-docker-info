module MultiScaleMapSampler
using JSON
using SimpleWeightedGraphs
using Graphs
using RandomNumbers
using LinearAlgebra
using Hungarian
import Combinatorics

export AbstractGraph,
    BaseGraph,
    MultiLevelGraph,
    MultiLevelSubGraph,
    MultiLevelPartition,
    edge_weight,

    # proposals
    build_single_node_flip,
    build_forest_recom2,
    build_forest_recom3,

    # constraints
    initialize_constraints,
    add_constraint!,
    AbstractConstraint,
    PopulationConstraint,
    ConstrainDiscontinuousTraversals,
    PackNodeConstraint,
    MaxCoarseNodeSplits,
    MaxSharedCoarseNodes,
    MaxSharedNodes,
    AllowedExcessDistsInCoarseNodes,
    MaxHammingDistance,
    
    MultiScaleCuttableTree,

    # Writer
    Writer,
    close_writer,
    push_writer!,

    # mcmc
    run_metropolis_hastings!,
    Measure,
    push_measure!,
    get_isoperimetric_score,
    get_isoperimetric_scores,
    get_polsby_popper_scores,
    get_log_spanning_trees,
    get_log_spanning_forests,
    get_log_linking_edges,
    get_cut_edge_count,
    get_area,
    get_perimeter,
    get_minofracs,
    get_node_counts,
    build_mcd_score,
    build_get_vra_score,

    # abstract chain
    Chain,
    run_chain!,

    # parallel tempering
    parallel_tempering!,
    parse_base_samples,
    parse_base_measure,

    get_vra_score,
    get_vra_scores,

    # cluster graph
    cluster_base_graph

include("./AtlasIO.jl")
include("./SimpleWeightedGraphs_BugFixes.jl")
include("./multi_level_graph.jl")
include("./subgraph.jl")
include("./node_set.jl")
include("./multi_level_subgraph.jl")
include("./multi_level_partition.jl")
include("./constraint_types.jl")
include("./construct_multi_level_partition.jl")
include("./balance_multi_level_graph.jl")
include("./polsby_popper.jl")
include("./node_counts.jl")
include("./mcd_ousted_population.jl")
include("./measure.jl")
include("./constraints.jl")
include("./forest_recom2.jl")
include("./forest_recom3.jl")
include("./single_node_flip.jl")
include("./writer.jl")
include("./mcmc.jl")
include("./chain.jl")
include("./parallel_tempering.jl")
include("./vap_frac.jl")
# include("./parallel_tempering_multiprocessing.jl")
include("./tree.jl")

end # module
