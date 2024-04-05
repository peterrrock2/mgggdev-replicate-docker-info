@testset "three_level test graph, 2 districts, pop=(5,7), gamma=1, compactness_weight=0.5" begin
    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(5,7))
    measure = Measure(1)
    compactness_weight = 0.5
    push_measure!(measure, get_isoperimetric_score, compactness_weight)

    n = 50000

    observed_districts = get_observed_districts(three_level_graph, constraints, three_level_dists, measure, n)

    c1 = 0
    c2 = 0

    for i in keys(observed_districts) 
        if all([x == nothing || all([y === nothing for y in values(x)]) for x ∈ values(i)])
            # 6/6 population split 
            c1 += observed_districts[i]
        else
            # 5/7 split 
            c2 += observed_districts[i]
        end 
    end
    c1w = exp(-compactness_weight * (12*12/6 + 12*12/6))
    c2w = exp(-compactness_weight * (12*12/7 + 10*10/5))
    tot = 4 * c1w + 8*c2w
    @show c1/n/2, 4*c1w/tot
    @test is_close(c1/n/2, 4*c1w/tot) 
    @test is_close(c2/n/2, 8*c2w/tot) 
end