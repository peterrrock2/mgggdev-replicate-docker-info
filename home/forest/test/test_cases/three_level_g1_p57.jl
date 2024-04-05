@testset "three_level test graph, 2 districts, pop=(5,7), gamma=1" begin
    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(5,7))
    measure = Measure(1)
    n = 10000

    observed_districts = get_observed_districts(three_level_graph, constraints, three_level_dists, measure, n)

    @test all([is_close(x, 1/12*2*n) for x in values(observed_districts)])
end