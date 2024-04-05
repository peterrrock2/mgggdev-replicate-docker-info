@testset "small_square test graph, 2 districts, pop=8, gamma=1" begin
    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(8,8))
    n = 1000
    measure = Measure(1)

    observed_districts = get_observed_districts(small_square_graph, constraints, small_square_dists, measure, n)

    # sanity check that the total districts sampled is correct 
    @test sum(values(observed_districts)) == small_square_dists * n

    # for small_square with gamma=1, pop=8, the distribution should be uniform on the four possibilities. 
    expected_samples_per_district = n * small_square_dists / length(observed_districts)
    @test all([ is_close(expected_samples_per_district, x) for x in values(observed_districts)])
end