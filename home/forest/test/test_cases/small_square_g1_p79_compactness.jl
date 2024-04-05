@testset "small_square test graph, 2 districts, pop=(7,9), gamma=1, compactness_weight=0.5" begin

    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(7,9))

    measure = Measure(1)
    compactness_weight = 0.5
    push_measure!(measure, get_isoperimetric_score, compactness_weight)

    n = 10000

    observed_districts = get_observed_districts(small_square_graph, constraints, small_square_dists, measure, n);

    c1, c2, c3 = count_small_square_districts(observed_districts)
    # e^{ -w * sum_{d in districts} of perimeter(d)^2 / area(d) }
    c1w = exp(-compactness_weight * (12*12/8 + 12*12/8)) 
    c2w = exp(-compactness_weight * (14*14/9 + 14*14/7))
    c3w = exp(-compactness_weight * (14*14/9 + 12*12/7))

    # with gamma=1, we expect totals proportional to the number of each type of districting. 
    # there are 2 8/8 splits and 8 of the other two. 
    total = 2 * c1w + 8*c2w + 8*c3w
    
    @test is_close(c1/n/2, 2*c1w/total) 
    @test is_close(c2/n, 8*c2w/total) 
    @test is_close(c3/n, 8*c3w/total) 
    
end 