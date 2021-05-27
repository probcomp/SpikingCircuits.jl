function flatten!(arr, map, par, c)
    push!(arr, c)
    map[par] = length(arr)
end

function flatten!(arr, map, par, c::CompositeComponent)
    for (key, subcomp) in pairs(c.subcomponents)
        flatten!(arr, map, (par..., key), subcomp)
    end
end

function flatten_to_subcomponents(c::CompositeComponent)
    arr = []
    map = Dict{Tuple, Int}()
    for (key, subcomp) in pairs(c.subcomponents)
        flatten!(arr, map, (key, ), subcomp)
    end
    return (arr, map)
end

function to_indexed_vals(val::CompositeValue)
    keys = keys_deep(val)
    vals = [val[key] for key in keys]
    old_key_to_idx = Dict(key => i for (i, key) in enumerate(keys))
    return (vals, old_key_to_idx)
end

unnest(p) = (p, )
unnest(p::Pair) = (p.first, unnest(p.second)...)

new_name(in::Input, _, old_to_new, _) = Input(old_to_new[in.id])
new_name(out::Output, _, _, old_to_new) = Output(old_to_new[out.id])
new_name(ci::CompIn, internal_idxs, old_to_new, _) = begin
    fst = ci.comp_name
    snd = ci.in_name
    addr = unnest(fst => snd)
    ind = internal_idxs[addr]
    CompIn(ind)
end
new_name(co::CompOut, internal_idxs, _, old_to_new) = begin
    fst = ci.comp_name
    snd = ci.out_name
    addr = unnest(fst => snd)
    ind = internal_idxs[addr]
    CompOut(ind)
end

function flatten(c::CompositeComponent)
    flattened_subcomponents, old_name_to_new_idxs = flatten_to_subcomponents(c)
    indexed_inputs, old_input_to_new_input = to_indexed_vals(inputs(c))
    indexed_outputs, old_output_to_new_output = to_indexed_vals(outputs(c))
    _new_name(nodename) = new_name(nodename, old_name_to_new_idxs, old_input_to_new_input, old_output_to_new_output)
    new_edges = (_new_name(src) => _new_name(dst) for (src, dst) in get_edges(c))
    inval = CompositeValue(indexed_inputs)
    outval = CompositeValue(indexed_outputs)
    return CompositeComponent(inval, 
                              outval,
                              flattened_subcomponents, 
                              new_edges)
end

function flatten!(arr, map, par, c::CompositeTrajectory)
    for (key, sub) in pairs(c.subtrajectories)
        flatten!(arr, map, (par..., key), sub)
    end
end

function flatten(c::CompositeTrajectory)
    arr = []
    map = Dict{Int, Tuple}()
    for (key, sub) in pairs(c.subtrajectories)
        flatten!(arr, map, (key, ), sub)
    end
    return (CompositeTrajectory(arr, 
                                c.trajectory_length, 
                                c.has_next_spike, 
                                c.next_spike_name), map)
end

function flatten!(arr, map, par, c::CompositeState)
    for (key, sub) in pairs(c.substates)
        flatten!(arr, map, (par..., key), sub)
    end
end

function flatten(c::CompositeState)
    arr = []
    map = Dict{Int, Tuple}()
    for (key, sub) in pairs(c.substates)
        flatten!(arr, map, (key, ), sub)
    end
    return (CompositeState(arr), map)
end
