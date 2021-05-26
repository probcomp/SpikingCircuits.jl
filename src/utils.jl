function flatten!(arr, map, par, c)
    push!(arr, c)
    map[length(arr)] = par
end

function flatten!(arr, map, par, c::CompositeComponent)
    for (key, subcomp) in pairs(c.subcomponents)
        flatten!(arr, map, (par..., key), subcomp)
    end
end

function flatten(c::CompositeComponent)
    arr = []
    map = Dict{Int, Tuple}()
    for (key, subcomp) in pairs(c.subcomponents)
        flatten!(arr, map, (key, ), subcomp)
    end
    return (arr, map)
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
