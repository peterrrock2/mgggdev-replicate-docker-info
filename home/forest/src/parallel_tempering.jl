import Base.Threads.@spawn

function parallel_tempering!(
	replicas::Vector{MultiLevelPartition},
	chains::Vector{Chain},
	steps::Int,
	swap_interval::Int,
	base_sampler::Union{Nothing, Tuple{Measure, Vector}}=nothing
)
	threads = Threads.nthreads()
	@assert length(replicas) == length(chains)
	@assert length(replicas) == threads
	initial_step, final_step = set_step_bounds(steps)

	swaps = floor((final_step-initial_step)/swap_interval)

# #######
	# extend partitions
	for ii = 1:length(replicas)
		replicas[ii].extensions[replica_id::EXTENSIONS] = ii
		replicas[ii].extensions[bath_swaps::EXTENSIONS] = 0
		num_dists = replicas[ii].num_dists
		replicas[ii].extensions[del_dists::EXTENSIONS]=zeros(num_dists)
	end
# #######
	error_msg = zeros(Threads.nthreads(), 2)

	pair_start, no_tasks = 0, 0
	tasks = Vector(undef, length(replicas))
	# @show length(tasks)

	for swap = 1:swaps
		tasks = []
		for ii = 1:length(replicas)
			inc = (Int((swap-1)*swap_interval), Int(swap*swap_interval))
			# tasks[ii] = Base.Threads.@spawn run_chain!(replicas[ii], chains[ii], inc)
			push!(tasks, @spawn run_chain!(replicas[ii], chains[ii], inc))
			# msg = run_chain!(replicas[ii], chains[ii], inc)
		 	# error_msg[ii, 1] = msg[1]
		 	# error_msg[ii, 2] = msg[2]
		end
		for ii = 1:length(replicas)
		 	msg = fetch(tasks[ii])
		 	error_msg[ii, 1] = msg[1]
		 	error_msg[ii, 2] = msg[2]
		end
		if maximum(error_msg[:, 1]) > 0
			println("erroring out")
			@show error_msg
			break
		end

		println("swap ", swap, " out of ", swaps); flush(stdout)
		pair_start += 1
		task_ind = 1
		tasks = []
		if base_sampler != nothing && pair_start == 2
			# tasks[task_ind] = Base.Threads.@spawn try_swap_replicas!(replicas, chains, 1,
			# 	                                        base_sampler)
			push!(tasks, @spawn try_swap_replicas!(replicas, chains, 1, Int(swap), 
				                                   base_sampler))
			task_ind += 1
		end
		for ii = pair_start:2:length(replicas)-1
			# tasks[task_ind] = Base.Threads.@spawn try_swap_replicas!(replicas, chains, ii)
			push!(tasks, @spawn try_swap_replicas!(replicas, chains, ii))
			task_ind += 1
		end
		for ii = 1:length(tasks)
		 	wait(tasks[ii])
		end

		pair_start = mod(pair_start, 2)
	end
end

""""""
function try_swap_replicas!(
	replicas::Vector{MultiLevelPartition},
	chains::Vector{Chain},
	index::Int
)
	@assert index < length(chains)
	@assert length(replicas) == length(chains)

	log_p_ii = log_measure(replicas[index], chains[index].measure)
	log_p_jj = log_measure(replicas[index+1], chains[index+1].measure)
	log_p_ij = log_measure(replicas[index], chains[index+1].measure)
	log_p_ji = log_measure(replicas[index+1], chains[index].measure)
	accept_prob = exp(log_p_ji + log_p_ij - log_p_ii - log_p_jj)

	println("try_swap_replicas: ", index, " ", accept_prob, " ", [log_p_ii, log_p_jj, log_p_ij, log_p_ji])
	if rand(chains[index].rng) < accept_prob
		tmp = replicas[index]
		replicas[index] = replicas[index+1]
		replicas[index+1] = tmp
	end
	# re-establish measures on replicas and chains, wherever they are now
	get_log_energy(replicas[index], chains[index].measure)
	get_log_energy(replicas[index+1], chains[index+1].measure)
end

function try_swap_replicas!(
	replicas::Vector{MultiLevelPartition},
	chains::Vector{Chain},
	index::Int,
	swap_ind::Int,
	base_sampler::Tuple{Measure, Vector}
)
	# @assert index < length(chains)
	# @assert length(replicas) == length(chains)

	base_measure, base_samples = base_sampler
	# base_ind = Int(floor(rand(chains[index].rng)*length(base_samples))) + 1
	base_map = base_samples[swap_ind]
	base_partition = MultiLevelPartition(replicas[index].graph, base_map)

	log_p_ii = log_measure(replicas[index], chains[index].measure)
	log_p_jj = log_measure(base_partition, base_measure)
	log_p_ij = log_measure(replicas[index], base_measure)
	log_p_ji = log_measure(base_partition, chains[index].measure)
	accept_prob = exp(log_p_ji + log_p_ij - log_p_ii - log_p_jj)

	println("try_swap_replicas: ", index, " ", accept_prob, " ", [log_p_ii, log_p_jj, log_p_ij, log_p_ji])
	if rand(chains[index].rng) < accept_prob
		base_partition.extensions[replica_id::EXTENSIONS] = 
		                      replicas[index].extensions[replica_id::EXTENSIONS]
		base_partition.extensions[bath_swaps::EXTENSIONS] = 
		                    replicas[index].extensions[bath_swaps::EXTENSIONS]+1
		base_partition.extensions[del_dists::EXTENSIONS] =
		                       replicas[index].extensions[del_dists::EXTENSIONS]
		replicas[index] = base_partition
	end
	# re-establish measures on replicas and chains, wherever they are now
	get_log_energy(replicas[index], chains[index].measure)
end


function parse_base_measure(
	base_sample_path::String,
	graph::MultiLevelGraph;
	ideal_pop::Union{Nothing,Real}=nothing
)
	io=smartOpen(base_sample_path, "r")
    atlas=openAtlas(io);
    gamma = atlas.atlasParam["gamma"]
    energy_weights = atlas.atlasParam["energy weights"]
    energies = atlas.atlasParam["energies"]
    m = nextMap(atlas)
    close(atlas)

    base_measure = Measure(gamma)
    for (e, w) in zip(energies, energy_weights)
    	if e == "get_isoperimetric_score"
    		push_measure!(base_measure, get_isoperimetric_score, w)
    	elseif e == "get_mcd_score"
    		eg_partition = MultiLevelPartition(graph, m.districting)
    		mcd_score = build_mcd_score(eg_partition, 
    			                        ideal_pop=ideal_pop)
    		push_measure!(base_measure, mcd_score, w, desc="get_mcd_score")
    	else
    		println("Unknown measure ", e)
    		@assert false
    	end
    end

    return base_measure
end

function parse_base_samples(base_sample_path::String, burn_in::Int,
	                        samples::Int, rng::AbstractRNG;
	                        manual_max_plans::Int=-1)
	# tot_plans = 0#countlines(base_sample_path)-10-burn_in
	io=smartOpen(base_sample_path, "r")
    atlas=openAtlas(io);
    eachLine=eachline(atlas.io)
    tot_plans = 0
    for l in eachLine 
    	tot_plans += 1
    end
    tot_plans -= 10 + burn_in
    @show tot_plans
    close(io)

	# if manual_max_plans > 0
	# 	tot_plans = min(tot_plans, manual_max_plans-10-burn_in)
	# end
	# @show tot_plans
	base_samples = Vector{Dict{Tuple{Vararg{String}}, Int}}(undef, samples)
	
	plan_inds = [Int(floor(rand(rng)*tot_plans))+burn_in+1 for jj=1:samples]
	ordered_plan_perm = sortperm(plan_inds)

	io=smartOpen(base_sample_path, "r")
    atlas=openAtlas(io);

    skipMap(atlas::Atlas; numSkip=burn_in)
    cur_plan_ind, cur_perm_ind = burn_in, 1
    next_plan_ind = plan_inds[ordered_plan_perm[cur_perm_ind]]
    gathered_samples = 0
    while gathered_samples < samples
    	# @show tot_plans, cur_plan_ind, next_plan_ind, gathered_samples, samples
	    m = nextMap(atlas)
	    cur_plan_ind += 1
	    while cur_plan_ind == next_plan_ind
		    gathered_samples += 1
	    	base_samples[ordered_plan_perm[cur_perm_ind]] = m.districting
	    	tmp = ordered_plan_perm[cur_perm_ind]
		    # @show gathered_samples, cur_perm_ind, tmp, plan_inds[tmp], cur_plan_ind 
	    	cur_perm_ind += 1
	    	if cur_perm_ind > length(ordered_plan_perm)
	    		break
	    	end
		    next_plan_ind = plan_inds[ordered_plan_perm[cur_perm_ind]]
		end
	end

	close(atlas)
    return base_samples
end