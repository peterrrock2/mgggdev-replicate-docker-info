@testset "three_level test graph, 2 districts, pop=(6,6), gamma=1" begin
    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(6,6))
    measure = Measure(1)
    n = 10000

    observed_districts = get_observed_districts(three_level_graph, constraints, three_level_dists, measure, n)

    # all partitions are equally likely 
    @test all([is_close(x, 1/4*2*n) for x in values(observed_districts)])
end