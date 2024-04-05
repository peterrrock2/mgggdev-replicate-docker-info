@testset "three_level test graph, 2 districts, pop=(5,7), gamma=0" begin
    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(5,7))
    measure = Measure(0)
    n = 10000

    observed_districts = get_observed_districts(three_level_graph, constraints, three_level_dists, measure, n)

    #all partitions have 8 possible hierarchhical trees + linking edges 
    @test all([is_close(x, 1/12*2*n) for x in values(observed_districts)])
end