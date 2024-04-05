@testset "small_square test graph, 2 districts, pop=(7,9), gamma=0" begin
    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(7,9))
    measure = Measure(0)
    n = 10000

    observed_districts = get_observed_districts(small_square_graph, constraints, small_square_dists, measure, n)
    
    c1, c2, c3 = count_small_square_districts(observed_districts)

    # the three instances should have approximately a 4:2:1 observation ratio. (c1 double counts)
    @test is_close(c1/n/2, 4/7) 
    @test is_close(c2/n, 1/7) 
    @test is_close(c3/n, 2/7) 
end 