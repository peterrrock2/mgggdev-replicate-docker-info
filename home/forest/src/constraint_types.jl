abstract type AbstractConstraint end

# @enum CONSTRAINT pop pack_node contiguous_traversal max_coarse_splt max_shared_coarse excess_dists_in_coarse max_hamming_distance max_shared_nodes

struct PopulationConstraint <: AbstractConstraint
    min_pop::Real
    max_pop::Real
end

struct PackNodeConstraint <: AbstractConstraint
    nodes::Union{Nothing, Dict{Tuple{Vararg{String}}, Int}}
    ideal_pop::Real
end

struct ConstrainDiscontinuousTraversals <: AbstractConstraint
    max_line_segments::Int
end

struct MaxCoarseNodeSplits <: AbstractConstraint
    max_coarse_node_splits::Int # how many coarse node splits are allowed?
end

struct MaxSharedCoarseNodes <: AbstractConstraint
    max_shared_coarse_nodes::Int # how many coarse node splits are allowed?
end

struct MaxSharedNodes <: AbstractConstraint
    max_shared_nodes::Int # how many quotient node splits are allowed?
end

struct AllowedExcessDistsInCoarseNodes <: AbstractConstraint
    excess_splitting::Int
    ideal_pop::Real
end

struct MaxHammingDistance <: AbstractConstraint
    partition::AbstractPartition # initial partition
    max_distance::Real # maximum fraction of deviation in a single distrct
    norm_type::String # 1 or infinity
end
