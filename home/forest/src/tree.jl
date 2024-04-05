"""
    self_avoiding_walk(g, s, niter, distmx=weights(g); seed=-1, rng=GLOBAL_RNG)

Perform a [loop-erased random walk](https://en.wikipedia.org/wiki/Loop-erased_random_walk)
on graph `g` starting at vertex `s` and continuing for a maximum of `niter` steps.
Return a vector of vertices visited in order.
"""
function loop_erased_randomwalk(
    g::AG, s::Integer, 
    niter::Integer=10000*nv(g)^2,
    distmx::AbstractMatrix{T}=weights(g); 
    f::BitArray, 
    walk_buff::Vector{Int}, # seed::Int=-1,
    rng::AbstractRNG=GLOBAL_RNG,
) where {T <: Real, U, AG <: Graphs.AbstractGraph{U}}
    (s > 0 && s <= nv(g)) || throw(BoundsError())
    # @assert length(walk_buff) == nv(g)

    walk_view = view(walk_buff, 1:1)
    walk_view[1] = s
    i = 1
    cur_pos = 1
    while i <= niter
        cur = walk_view[cur_pos]
        if f[cur]
            break
        end
        nbrs = neighbors(g, cur)
        # length(nbrs) == 0 && throw(ArgumentError())
        wght = [distmx[cur, n] for n in nbrs]
        v = nbrs[findfirst(cumsum(wght) .> rand(rng)*sum(wght))]
        if v in walk_view
            cur_pos = indexin(v, walk_view)[1]
            walk_view = view(walk_buff, 1:cur_pos)
        else
            cur_pos += 1
            walk_buff[cur_pos] = v
            walk_view = view(walk_buff, 1:cur_pos)
        end
        i += 1
    end

    if !f[walk_view[cur_pos]]  
        #check connectivity of graph here
        if typeof(g) == SimpleWeightedGraph
            connected = is_connected_bf(g)
        else
            connected = is_connected(g)
        end
        if connected
            throw(ErrorException("termiating set was not reached; graph *is* connected.  Try increasing niter"))
        else
            throw(ErrorException("termiating set was not reached; graph is *NOT* connected"))
        end
    end
    return walk_view
end

"""
    wilson_rst(g...)

Randomly sample a spanning tree uniformly from a graph g, using rng as input to
Wilson's algorithm. Edge weights may be given using edge_weights and vertex neighbors
are sampled proportional to these weights; otherwise the neighbors are sampled uniformly
at random.
"""
# function wilson_rst end
# @traitfn 
function wilson_rst(g::SimpleWeightedGraph,
    rng::AbstractRNG=GLOBAL_RNG,
    distmx::AbstractMatrix{T}=SimpleWeightedGraphs.weights(g)
) where T <: Real#, U, AG <: SimpleWeightedGraphs.AbstractGraph{U}}
    visited_vertices = BitArray([0 for ii = 1:nv(g)])
    start1 = rand(rng, 1:nv(g))
    start2 = rand(rng, 1:nv(g))
    visited_vertices[start2] = 1
    walk_buff = zeros(Int, nv(g))
    
    walk = loop_erased_randomwalk(g, start1, f=visited_vertices, 
                                  walk_buff=walk_buff, rng=rng)
    edges = Vector{Edge}(undef, nv(g)-1)
    cur_edge = 1
    for i = 1:length(walk)-1
        edges[cur_edge] = Edge(walk[i], walk[i+1])
        cur_edge += 1
    end

    for vw in walk
        visited_vertices[vw] = 1
    end
    cur_visited = length(walk)

    while cur_visited < nv(g)
        v = rand(rng, 1:(nv(g)-cur_visited))
        v = findnext(iszero, visited_vertices, v)
        walk = loop_erased_randomwalk(g, v, f=visited_vertices, 
                                      walk_buff=walk_buff, rng=rng)
        for vw in walk
            visited_vertices[vw] = 1
        end
        for i = 1:length(walk)-1
            edges[cur_edge] = Edge(walk[i], walk[i+1])
            cur_edge += 1
        end
        cur_visited += length(walk)-1
    end
    return edges
end


""""""
function log_nspanning(
    g::AG
) where {U, AG <: Graphs.AbstractGraph{U}}
    # no support for logdet of sparse matrices yet, so must cast to full Matrix
    return logdet(Matrix(view(laplacian_matrix(g), 2:nv(g), 2:nv(g))))
end
