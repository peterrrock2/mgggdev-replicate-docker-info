#=import Pkg
codepath = joinpath("/Users/gabrielchuang/Documents/Gabriel/Duke/MultiScaleMapSampler") 
Pkg.activate("mergeSplit"; shared=true)
using Revise

Pkg.develop(path=codepath) =#

using MultiScaleMapSampler
using Test
using RandomNumbers

const testdir = dirname(@__FILE__)

function is_close(a,b)
    if a > 0.01
        0.9 <= a/b && a/b <= 1.1
    else 
        0.6 <= a/b && a/b <= 1.4
    end 
end

function get_observed_districts(
    graph::MultiLevelGraph, 
    constraints::Dict,
    num_districts::Int, 
    measure::Measure, 
    n::Int
)::Dict
    rng = PCG.PCGStateOneseq(UInt64, 1241909)
    partition = MultiLevelPartition(graph, constraints, num_districts; rng=rng);
    proposal = build_forest_recom2(constraints)
    observed_districts = Dict()
    for ii = 1:n
        run_metropolis_hastings!(partition, proposal, measure, 1, rng);
        district_to_nodes = partition.district_to_nodes
        for d2n in district_to_nodes
            if d2n in keys(observed_districts)
                observed_districts[d2n] += 1
            else 
                observed_districts[d2n] = 1 
            end 
        end
    end

    return observed_districts
end 

# given observed_districts, mapping districts to frequencies, count_small_square_districts returns the number of 
# each flavor of districts with population split in (7,9) - 8/8, 7/9 with nub inside, 7/9 with nub outside, respectively. 
# 
# . . . .   . . . .   . . . . 
# . . . .   . . X .   . . . X
# X X X X   X X X X   X X X X 
# X X X X   X X X X   X X X X 
function count_small_square_districts(observed_districts::Dict)::Tuple{Int,Int,Int}
    c1, c2, c3 = 0, 0, 0
    middle_square_values = [("0,0", "1,1"), ("0,1", "1,2"), ("1,0","2,1"),("1,1","2,2")]

    for i in keys(observed_districts) 
        if all([x === nothing for x ∈ values(i)])
            # 8/8 population split 
            c1 += observed_districts[i]
        elseif length(i) == 3 && any([x !== nothing && !(isempty(intersect(middle_square_values, keys(x)))) for x in values(i)])
            # 9/7 population split, with the nub on the inside 
            # i.e., if the 9-population district contains one of the middle squares as the only precinct in its county 
            c2 += observed_districts[i]
        elseif length(i) == 3
            # 9/7 population split, with the nub on the outside 
            c3 += observed_districts[i]
        end
    end
    return c1, c2, c3
end 

small_square_json = joinpath("test_graphs", "4x4pct_2x2cnty.json")
small_square_node_data = Set(["county", "pct", "pop", "area", "border_length"])
small_square_base_graph = BaseGraph(small_square_json, "pop", inc_node_data=small_square_node_data,
                                    area_col="area",node_border_col="border_length", 
                                    edge_perimeter_col="length")
small_square_graph = MultiLevelGraph(small_square_base_graph, ["county", "pct"])
small_square_dists = 2

three_level_json = joinpath("test_graphs", "2x6pct_3level.json")
three_level_node_data = Set(["county", "subcounty", "pct", "pop", "area", "border_length"])
three_level_base_graph = BaseGraph(three_level_json, "pop", inc_node_data=three_level_node_data,
                                    area_col="area",node_border_col="border_length", 
                                    edge_perimeter_col="length")
three_level_graph = MultiLevelGraph(three_level_base_graph, ["county", "subcounty", "pct"])
three_level_dists = 2

tests = [
    "small_square_g1_p88", #small square graph, gamma=1, population=(8,8)
    "small_square_g0_p79", #small square graph, gamma=0, population=(7,9)
    "small_square_g1_p79", #small square graph, gamma=1, population=(7,9)
    "small_square_g1_p79_compactness", #small square graph, gamma=1, population=(7,9), compactness_weight=0.5 
    "small_square_g1_p79_SNF", #small square graph, gamma=1, population=(7,9), single node flip
    "three_level_g1_p66", #3-level graph, gamma=1, population=(6,6)
    "three_level_g0_p57", #3-level graph, gamma=0, population=(5,7)
    "three_level_g1_p57", #3-level graph, gamma=1, population=(5,7)
    "three_level_g1_p57_compactness", #3-level graph, gamma=1, population=(5,7), compactness_weight=0.5 
    "three_level_g1_p57_SNF" # single node flip (very basic check only)
    ]

for t in tests
    tp = joinpath(testdir, "test_cases/$(t).jl")
    include(tp)
end