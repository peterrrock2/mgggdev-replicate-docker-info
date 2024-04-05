@testset "small_square test graph, 2 districts, pop=(7,9), gamma=1" begin

    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(7,9))
    n = 10000
    measure = Measure(1)

    observed_districts = get_observed_districts(small_square_graph, constraints, small_square_dists, measure, n);

    c1, c2, c3 = count_small_square_districts(observed_districts)

    # these fail (which is incorrect; they should pass)
    @test is_close(c1/n/2, 2/18) 
    @test is_close(c2/n, 8/18) 
    @test is_close(c3/n, 8/18) 
end 