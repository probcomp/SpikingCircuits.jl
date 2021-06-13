# The flattening algorithm:
#
# Each CompositeComponent has some set of Inputs and Outputs.
# As well as some set of "inner" CompIns and CompOuts.
# 
# The algorithm is recursive. 
#
# For each Input `I` of the composite component, there is some set of edges from `I` to CompIn (or directly to Output).
# If edge E is an edge from `I` to a CompIn `C_in`
# We must associate new edges from `I` to sub-components which are
# routed to by `C_in`.
#
# Let's consider a "1-level" flattening where we have:
#
#     |             _ _ _ _ _ _     |
#     |            |     prim  |    |
# (I) * -> (C_in)  *  -> |> -> * -> * (O)
#     |            |_ _ _ _ _ _|    |
#     |                             |
#
# The new edge for `I` must go to `prim` directly.
# Interface functions used:
# 
# - inputters: for a `NodeName` name (which is the name of a component in a larger component), 
# iterator over all other `NodeName` instances which send an input to name.
# 
# - receivers: for a `NodeName` name, iterator over all other `NodeName` instance which receive
# from name.

function to_indexed_vals(val::CompositeValue)
    keys = keys_deep(val)
    vals = [val[key] for key in keys]
    old_key_to_idx = Dict(key => i for (i, key) in enumerate(keys))
    return (vals, old_key_to_idx)
end

struct FlatBuilder
    inp_indexed_vals
    inp_old_key_to_idx
    out_indexed_vals
    out_old_key_to_idx
    new_input_edges
    new_output_edges
end

function FlatBuilder(c::CompositeComponent)
    FlatBuilder(to_indexed_vals(inputs(c))..., 
                to_indexed_vals(outputs(c))...,
                [], [])
end

function build(fb::FlatBuilder)
end

function handle!(fb::FlatBuilder, c::CompositeComponent, 
        outp::Output, k::CompOut)
    sub = c[k.comp_name]
    recs = inputters(sub, Input(Circuits.valname(k)))
    append!(fb.new_input_edges, map(recs) do r
                r => outp
            end)
end

function handle!(fb::FlatBuilder, c::CompositeComponent, 
        inp::Input, k::CompIn)
    sub = c[k.comp_name]
    sends = receivers(sub, Input(Circuits.valname(k)))
    append!(fb.new_input_edges, map(sends) do s
                inp => s
            end)
end

function create_new_edges!(fb::FlatBuilder, c::CompositeComponent, inp::Input)
    for k in receivers(c, inp)
        handle!(fb, c, inp, k)
    end
end

function create_new_edges!(fb::FlatBuilder, c::CompositeComponent, outp::Output)
    for k in inputters(c, outp)
        handle!(fb, c, outp, k)
    end
end

function flatten!(fb::FlatBuilder, c::CompositeComponent)
    for inp in filter(x -> x isa Input, c.idx_to_node)
        create_new_edges!(fb, c, inp)
    end
    for outp in filter(x -> x isa Output, c.idx_to_node)
        create_new_edges!(fb, c, outp)
    end
end

# Toplevel.
function flatten(c::CompositeComponent)
    inp_vals, inp_old_keys_to_idx = to_indexed_vals(inputs(c))
    out_vals, out_old_keys_to_idx = to_indexed_vals(outputs(c))
end

# Depth first flatten.
function dfs_flatten(c::CompositeComponent)
end

second(t::Tuple{T, K}) where {T, K} = t[2 : end]

# George.
function flatten(c::CompositeComponent)
    flt_subs_and_mappings = map(flatten, c.subcomponents)
    flt_subs = map(first, flt_subs_and_mappings)
    subcomp_mappings = map(second, flt_subs_and_mappings)
    # each element of subcomp_mappings is `input_map, output_map` giving the name
    # maps for a subcomponent
    idx_to_in = inputs(c) |> keys_deep |> collect
    idx_to_out = outputs(c) |> keys_deep |> collect
    in_to_idx = Dict(v => k for (k, v) in enumerate(idx_to_in))
    out_to_idx = Dict(v => k for (k, v) in enumerate(idx_to_out))
    idx_to_subsubcomp_name = []
    for (name, subcomp) in pairs(c.subcomponents)
        for (innername, _) in pairs(subcomp.subcomponents)
            push!(idx_to_subsubcomp_name, (name, innername))
        end
    end
    subsubcomp_name_to_idx = Dict(v => k for (k, v) in enumerate(idx_to_subsubcomp_name))
    new_edges = []
    for (src, dst) in get_edges(c)
        srcnames = new_names(src, c, subsubcomp_name_to_idx, in_to_idx, out_to_idx, subcomp_mappings)
        dstnames = new_names(src, c, subsubcomp_name_to_idx, in_to_idx, out_to_idx, subcomp_mappings)
        for (s, d) in Iterators.product(srcnames, dstnames)
            push!(new_edges, s => d)
        end
    end
    new_comp = CompositeComponent(
        IndexedValue([inputs(c)[name] for name in idx_to_in]),
        IndexedVlaue([outputs(c)[name] for name in idx_to_out]),
        flt_subs,
        new_edges
    )
    return (new_comp, in_to_idx, out_to_idx)
end
