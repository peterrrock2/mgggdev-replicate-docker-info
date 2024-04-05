@testset "small_square test graph, 2 districts, pop=(7,9), gamma=1, single node flip" begin
    # this test case is VERY BASIC: it just checks that single node flip does in fact sample every possible partition
    # not very sure what the best way to go about checking that it's approximately uniform-ish. 
        constraints = initialize_constraints()
        add_constraint!(constraints, PopulationConstraint(7,9))
    
        measure = Measure(1)
        n = 10000
    
        rng = PCG.PCGStateOneseq(UInt64, 1241909)
        partition = MultiLevelPartition(small_square_graph, constraints, small_square_dists; rng=rng);
        proposal = build_single_node_flip(constraints)
        observed_districts = Dict()
        @time for ii = 1:n
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
    
        @test length(observed_districts) == 140 + 136 + 136
        # there are 140 8/8 splits, 136 7/9 splits. 
    end 