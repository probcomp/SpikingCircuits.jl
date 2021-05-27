"""
    SpikingSimulator

A event-based simulator for the `Spiking` target.

At any time during simulation, each component has a `State`
and a `Trajectory`.  A `t::Trajectory` describes how the state will evolve
over the next `trajectory_length(t)` milliseconds if no input spikes are received in this time.
A trajectory may specify the next output spike which will be emmitted from this component if no inputs are received.
`next_spike(::Component, t::Trajectory)` returns `nothing` if the trajectory does not extend to a time
when the component emits an output spike, or `output_name` if it does.
If the trajectory extends to a spike emission, `trajectory_length(t)` is the time of this spike emission.

The function `extend_trajectory!(::Component, ::State, t::Trajectory)` returns the a new `trajectory`
such that `trajectory_length(trajectory) > trajectory_length(t)`.  [TODO: do we want `>` or `≥`?]

The function `advance_time_by!(::Component, ::State, t::Trajectory, ΔT)` returns
a tuple `(new_state, new_traj, output_spike_names)` giving a new state and trajectory
for the component obtained by advancing the time by `ΔT`, and an iterator over the names
of the outputs from this component which spiked during this time.  It is required
that `ΔT ≤ trajectory_length(t)`, so we are only advancing time to a point we have already
determined in the trajectory.  Note that this means that any spikes emitted by the component
are emitted at exactly time `ΔT`, since we require that any output spikes the trajectory `t`
includes occur at `trajectory_length(t)`.

The function `receive_input_spike!(::Component, ::State, ::Trajectory, inname)`
returns a tuple `(new_state, new_traj, spiking_output_names)` of a new state and (possibly empty)
trajectory reflecting knowledge of the inputted spike.  `spiking_output_names` is an iterator
over the names of the component outputs which spike immediatly upon receiving this input spike.
This will often return an empty trajectory, since the old trajectory is rendered invalid by the knowledge
that a new spike arrived.

Currently the simulator operates on immutable objects; it may be more performant to move to mutable objects.

Note: I think that the code probably requires some condition about the order in which input spikes are received
and output spikes are emitted in one moment of time not mattering.  It may be worthwhile to articulate this at some point.
It may also be the right call to move to a simulator where nothing can happen at the same instant, and when things in theory happen
at the same instant, we instead sequence them with some tiny ϵ-milisecond delay instead.
"""
module SpikingSimulator

using Printf: @sprintf
using DataStructures: OrderedDict, Queue, enqueue!, dequeue!
import Circuits: Component, PrimitiveComponent, CompositeComponent, CompositeValue
import Circuits: Output, Input, CompOut, CompIn, NodeName, receivers, does_output, keys_deep, inputs, outputs, get_edges
import ..Spiking

const PrimComp = PrimitiveComponent{>:Spiking}

# TODO: could we do this better/in a more specific way?  This currently probably enables _some_ type checking but not performance improvements.
const Name = Union{Integer, Symbol, Pair}

#############
# Interface #
#############

"""
    State

A spiking component state.
"""
abstract type State end

"""
    Trajectory

A spiking component trajectory.

A `t::Trajectory` describes how the state will evolve
over the next `trajectory_length(t)` milliseconds, assuming that no input spikes are received in this time.

Invariant: if no input spikes arrive in the next `trajectory_length(t)` milliseconds, no output spikes are emitted
until `trajectory_length(t)` milliseconds have passed.  (Ie. the soonest an output spike could occur is in exactly `trajectory_length(t)` ms.)
"""
abstract type Trajectory end

"""
    initial_state(::Component)

The initial `State` for a spiking component.
"""
initial_state(::Component) = error("Not implemented.")

"""
    empty_trajectory(::Component)

A `t::Trajectory` for the spiking component such that `trajectory_length(t) == 0`.
"""
empty_trajectory(::Component) = EmptyTrajectory()

"""
    next_spike(::Component, t::Trajectory)

If the trajectory does not extend to a time when the component outputs a spike, returns `nothing`.
If the trajectory does extend to a time when the component outputs a spike, returns the name of the component
output which spikes at `trajectory_length(t)`.
"""
# TODO: do we want `next_spike`? or something like `next_spikes` in case multiple spikes happen at once?
next_spike(::Component, ::Trajectory) = error("Not implemented.")

"""
    extended_trajectory = extend_trajectory!(::Component, ::State, old_trajectory::Trajectory)

Return a new trajectory for the component such that `trajectory_length(extended_trajectory) > trajectory_length(old_trajectory)`.
May mutate `old_trajectory`.

Input condition: `next_spike(::Component, old_trajectory)` must be `nothing`.  Ie. we cannot extend a trajectory
past the next output spike.
"""
extend_trajectory!(::Component, ::State, ::Trajectory) = error("Not implemented.")

function advancing_too_far_error(t, ΔT)
    @assert ΔT > trajectory_length(t) "This error should only be called when `ΔT > trajectory_length(t)`!  (But `ΔT = $ΔT` and `trajectory_length(t)` = $(trajectory_length(t))"
    error("`advance_by_time` called to advance a trajectory past its length.  (Asked to extend a trajectory of length $(trajectory_length(t)) by $ΔT seconds.")
end

"""
    (new_state, new_traj, spiking_output_names) = advance_time_by!(::Component, ::State, t::Trajectory, ΔT)

Returns a new state and trajectory obtained by advancing along the current trajectory by `ΔT`.  Also returns
the iterator `spiking_output_names` over the names of the component's outputs which spike within those `ΔT` seconds.
May mutate the old state and trajectory.

Input condition: `trajectory_length(t) ≥ ΔT`.

Note that this condition, and the invariant that a trajectory length cannot extend past the next output spike,
implies that if `spiking_output_names` is not empty, then all spikes in `spiking_output_names`
occurred at exactly `ΔT` seconds.
(Ie. if `Ts` is the time of the spike output, the trajectory length invariant tells us
`Ts ≥ trajectory_length(t)`, and thus the input condition tells us `Ts ≥ ΔT`.
Thus, if `Ts ≤ ΔT` so that `spiking_output_names` is nonempty, we must have `Ts = ΔT`.)
"""
advance_time_by!(::Component, ::State, ::Trajectory, ΔT) = error("Not implemented.")

"""
    (new_state, new_traj, spiking_output_names) = receive_input_spike!(::Component, ::State, ::Trajectory, inname)

Returns a new state and trajectory consistent with the component receiving a spike in the input with name `inname`.
Also returns an iterator `spiking_output_names` over names of component outputs which spike immediately when
this input is received.
May mutate the new state and trajectory.

`new_traj` may be any valid trajectory for the component consistent with this spike being received.
(One common choice is to have `new_traj` simply be an empty trajectory.)
"""
receive_input_spike!(::Component, ::State, ::Trajectory, inname) = error("Not implemented.")

### General State/Trajectory Types ###

"""
    EmptyState <: State

Empty state containing no information.
"""
struct EmptyState <: State end

"""
    OnOffState <: State

State denoting that a component is either on or off.
`state.on` is a boolean returning whether it is on.
"""
struct OnOffState <: State
    on::Bool
end
Base.:(==)(a::OnOffState, b::OnOffState) = a.on == b.on

"""
    EmptyTrajectory <: Trajectory

A trajectory of length 0.
"""
struct EmptyTrajectory <: Trajectory end
Base.show(io::IO, t::EmptyTrajectory) = print(io, "EmptyTrajectory()")

"""
    NextSpikeTrajectory <: Trajectory

A trajectory containing the time of the next spike (given by
`traj.time_to_next_spike`).
"""
struct NextSpikeTrajectory <: Trajectory
    time_to_next_spike::Float64
end
Base.show(io::IO, t::NextSpikeTrajectory) = print(io, "NextSpikeTrajectory($(t.time_to_next_spike))")

next_spike(::Component, ::EmptyTrajectory) = nothing # no next spike is within this trajectory

trajectory_length(::EmptyTrajectory) = 0.
trajectory_length(t::NextSpikeTrajectory) = t.time_to_next_spike

advance_without_statechange(s::State, t::NextSpikeTrajectory, ΔT) =
    let remaining_time = t.time_to_next_spike - ΔT
        if remaining_time == 0
            (s, EmptyTrajectory(), (:out,))
        elseif remaining_time > 0
            (s, NextSpikeTrajectory(remaining_time), ())
        else
            advancing_too_far_error(t, ΔT)
        end
    end

################################
# High-Level Simulator Methods #
################################

"""
    abstract type Event end

An event which can occur during a spiking simulation we may want to observe.

An `Event` occurs at some time for some component (but the time and component
are not part of the `Event` object).
"""
abstract type Event end

abstract type Spike <: Event end
struct InputSpike <: Spike
    name::Name
end
struct OutputSpike <: Spike
    name::Name
end

# """
#     Spike <: Event
#     Spike(outputname)

# The event that a spike was output by a component from the output with name `outputname`.
# """
# struct Spike <: Event
#     outputname::Name
# end

# """
#     StateChange <: Event
#     StateChange(new_state)

# The event that a component's state changed to `new_state`.
# """
# struct StateChange <: Event
#     new_state::State
# end

"""
    (new_state, new_traj, spiking_output_names) = advance_time_by!(c::Component, ::State, ::Trajectory, ΔT, f::Function)

Same as `advance_time_by!`, but calls `f(itr, dt)` at each `dt` at which an `Event` occurs
in the next `ΔT` seconds (including `ΔT`).  `itr` will be an iterator over pairs
`(compname, event)`, specifying each `event::Event` occurring for `c` or a subcomponent of `c`
specified by `compname` at time `dt`.  `compname` is a subcomponent name (possibly nested, or `nothing` to refer to `c`);
see documentation for `CompositeComponent`.

By default, for any component other than a `CompositeComponent`, it is assumed that there are no events occuring for internal components,
and so `f` is only called for events occurring at the top level to `c` (ie. output spikes from `c` and state changes for `c`).
"""
function advance_time_by!(c::Component, s::State, t::Trajectory, ΔT, f::Function)
    (newstate, t, son) = advance_time_by!(c, s, t, ΔT)
    f(((nothing, OutputSpike(name)) for name in son), ΔT)
    # f(
    #     Iterators.flatten((
    #         ((nothing, OutputSpike(name)) for name in son),
    #         newstate == s ? () : ((nothing, StateChange(newstate)),)
    #     )),
    #     ΔT
    # )

    return (newstate, t, son)
end

# TODO: don't repeat so much from the `advance_time_by` above
"""
    (new_state, new_traj, spiking_output_names) = receive_input_spike!(c::Component, s::State, t::Trajectory, inname, f::Function)

Same as `receive_input_spike!`, but calls `f(itr)` before returning, where `itr` is an iterator
over values `(compname, event)` specifying each `event::Event` occurring in `c` or a subcomponent
(specified by `compname`) when this input spike is received.  `compname` is a subcomponent name (possibly nested, or `nothing` to refer to `c`);
see documentation for `CompositeComponent`.

By default, for any component other than a `CompositeComponent`, it is assumed that there are no events occuring for internal components,
and so `f` is only called for events occurring at the top level to `c` (ie. output spikes from `c` and state changes for `c`).
"""
function receive_input_spike!(c::Component, s::State, t::Trajectory, inname, f::Function)
    (newstate, t, son) = receive_input_spike!(c, s, t, inname)
    f(
        Iterators.flatten((
            ((nothing, InputSpike(inname)),),
            ((nothing, OutputSpike(name)) for name in son) #,
            # newstate == s ? () : ((nothing, StateChange(newstate)),)
        ))
    )

    return (newstate, t, son)
end

# TODO: more general method with a way to give input spikes during the simulation.
# The main desideratum at this point is the ability to see the output spikes
# before deciding what the next input spikes are going to be.
# (This could probably be accomplished in a somewhat hacky way using the `inputs` iterator,
# but it would be better to have a well-designed mechanism.)

"""
    simulate_for_time(callback, c::Component, ΔT, s::State=initial_state(c), t::Trajectory=empty_trajectory(c);
        inputs=(), event_filter=((compname, event)->true)
    )

Simulates the component for `ΔT` milliseconds.
    
To send inputs into the circuit, a kwarg `inputs` may be provided.
This should be an iterator over `(time, itr)` pairs ordered by time,
where `itr` is an iterator over input names where a spike should be delivered at time `time` into the simulation.
(If an input name appears N times, N spikes will be delivered to that input.)

If all inputs should be delivered at time=0, instead of `inputs`, one can use kwarg `initial_inputs`
to provide an iterator over the input names to deliver spikes to at t=0.

Calls `callback(itr, time)` at every time since the start of the simulation
when any `Event`s occur occur for `c` or one of its subcomponents.
`itr` will be an iterator over pairs `(component_name, event::Event)`
giving each event to occur for `c` or one of its subcomponents at `time`.
Any `event` for a component with name `compname` such that `!event_filter(compname, event)` will be skipped.
Invariant: if `callback(itr1, t1)` is called before `callback(itr2, t2)`, then `t1 ≤ t2`.

`component_name` will be a subcomponent name (possibly nested; `nothing` to refer to `c`); see documentation
for `CompositeComponent`.

See also: `simulate_for_time_and_get_events`.

Note: this function assumes no additional inputs enter `c` after the `initial_inputs`. (So, eg., there must be no recurrent connections
back into `c`'s inputs.  If recurrent connections are needed, nest the recurrent component within an outer `CompositeComponent`.)
"""
function simulate_for_time(
    callback, c::Component, ΔT, s::State=initial_state(c), t::Trajectory=empty_trajectory(c);
    initial_inputs=(), inputs=[(0., initial_inputs)], event_filter=((compname, event)->true)
)
    function filtered_callback(itr, dt)
        if dt <= ΔT # only track events before the time limit
            callback(
                Iterators.filter(args -> event_filter(args...), itr),
                dt
            )
        end
    end

    next_inputs, inputs = isempty(inputs) ? ((Inf, ()), ()) : Iterators.peel(inputs)
    time_passed = 0.
    while time_passed < ΔT
        t = extend_trajectory!(c, s, t)
        extending_by = trajectory_length(t)
        # println("Extended by $extending_by")

        if trajectory_length(t) == Inf
            break;
        end
        time_to_next_input = next_inputs[1] - time_passed
        @assert time_to_next_input ≥ 0.
        if trajectory_length(t) ≥ time_to_next_input
            # println("Feeding in input[s] $(next_inputs[2])")
            extending_by = time_to_next_input
            (s, t, _) = advance_time_by!(c, s, t, time_to_next_input, (itr, t) -> filtered_callback(itr, t + time_passed))
            for input in next_inputs[2]
                (s, t, _) =
                    try
                        receive_input_spike!(c, s, t, input, itr -> filtered_callback(itr, next_inputs[1]))
                    catch e
                        @error("Error encountered when sending in input $input at time $(next_inputs[1]).", exception=(e, catch_backtrace()))
                        throw(e)
                    end
            end

            next_inputs, inputs = isempty(inputs) ? ((Inf, ()), ()) : Iterators.peel(inputs)
            println("Just delivered inputs at time $(time_passed + time_to_next_input).  Next inputs delivered at time $(next_inputs[1]).")
        else
            (s, t, _) = advance_time_by!(c, s, t, trajectory_length(t), (itr, t) -> filtered_callback(itr, t + time_passed))
        end
        time_passed += extending_by
    end

    return nothing
end
# simulate_for_time(
#     callback, c::Component, ΔT, s::State=initial_state(c), t::Trajectory=empty_trajectory(c);
#     initial_inputs, event_filter=((compname, event)->true)
# ) = simulate_for_time(callback, c, ΔT, s, t;
#     inputs = [(0, initial_inputs)], event_filter
# )

default_log_str(time, compname, event) =
    let timestr = @sprintf("%.4f", time)
        "$timestr : $compname $event"
    end

"""
    simulate_for_time_and_get_events(c::Component, ΔT, s::State=initial_state(c), t::Trajectory=empty_trajectory(c);
        log_filter=((_,_,_)->true), log=false, log_str=default_log_str, [inputs=() | initial_inputs=()], event_filter=((compname, evt)->true)
    )

Same as `simulate_for_time`, but returns a vector of triples `(time, component_name, event)`
giving the events which occurred for `c` and subcomponents during the simulation,
instead of using a callback function.
    
If `log` is true, will log every `(time, compname, event)` triple which passes `log_filter`.
The logging reports the string returned by `log_str(time, compname, event)`.

See also: `simulate_for_time`.
"""
function simulate_for_time_and_get_events(args...;
    log_filter=((_,_,_) -> true), log=false,
    log_str=default_log_str,
    kwargs...
)
    events = Tuple{Float64, Union{Nothing, Name}, Event}[]

    # function to add events to the array
    callback = (
        if !log
            (itr, time) -> begin
                for (compname, event) in itr
                    push!(events, (time, compname, event))
                end
            end
        else
            (itr, time) -> begin
                for (compname, event) in itr
                    push!(events, (time, compname, event))
                    
                    if log_filter(time, compname, event)
                        @info log_str(time, compname, event)
                    end    
                end
            end
        end
    )

    try
        simulate_for_time(callback, args...; kwargs...)
    catch e
        @error("Simulation terminated due to exception.", exception=(e, catch_backtrace()))
    end

    return events
end

# simulate_for_time_and_get_spikes_and_primitive_statechanges(c, args...; kwargs...) =
#     simulate_for_time_and_get_events(c, args..., kwargs..., event_filter= (compname, e) -> (
#         (e isa Spike) || (e isa StateChange && c[compname] isa PrimitiveComponent{>:Spiking})
#     ))

#############
# Composite #
#############

"""
    CompositeState <: State

State for a `CompositeComponent.`
"""
struct CompositeState{S} <: State
    substates::S
    CompositeState(substates::S) where {S <: Union{<:Vector, <:Dict}} = new{S}(substates)
end
Base.:(==)(a::CompositeState, b::CompositeState) = a.substates == b.substates

"""
    CompositeTrajectory <: Trajectory

Trajectory for a `CompositeComponent.`
"""
struct CompositeTrajectory{T} <: Trajectory
    subtrajectories::T
    trajectory_length::Float64
    has_next_spike::Bool
    next_spike_name::Union{Nothing, Name}
    CompositeTrajectory(st::T, args...) where {T <: Union{<:Vector, <:Dict}} = new{T}(st, args...)
end
function Base.show(io::IO, t::CompositeTrajectory)
    print(io, "CompositeTrajectory((")
    for (i, st) in enumerate(t.subtrajectories)
        print(io, st)
        if i != length(t.subtrajectories)
            print(io, ", ")
        end
    end
    print(io, "), $(t.trajectory_length), $(t.has_next_spike), $(t.next_spike_name))")
end

to_mutable_version_map(f, t::Tuple, T) = T[f(x) for x in t]
to_mutable_version_map(f, nt::NamedTuple, T) = Dict{Symbol, T}(key => f(val) for (key, val) in pairs(nt))

# state_type(...)
initial_state(c::CompositeComponent) = CompositeState(to_mutable_version_map(initial_state, c.subcomponents, State))
empty_trajectory(c::CompositeComponent) = CompositeTrajectory(to_mutable_version_map(empty_trajectory, c.subcomponents, Trajectory), 0.0, false, nothing)
trajectory_length(t::CompositeTrajectory) = t.trajectory_length
next_spike(::CompositeComponent, t::CompositeTrajectory) = t.has_next_spike ? t.next_spike_name : nothing

function first_pair_with_nonnothing_value(itr) # helper function
    for (name, v) in itr
        if v !== nothing
            return (name, v)
        end
    end
    return nothing
end
function extend_trajectory!(c::CompositeComponent, s::CompositeState, t::CompositeTrajectory)
    # If there are no subcomponents, this will never emit an output spike (unless just passing through an input),
    # so we can extend to `Inf` (until we receive an input spike)
    if isempty(c.subcomponents)
        @assert trajectory_length(t) == 0 || trajectory_length(t) == Inf
        return CompositeTrajectory(t.subtrajectories, Inf, false, nothing)
    end

    trajectories = t.subtrajectories
    mintime = Inf
    at_min_time = [] # TODO: is there a more performant structure than an array?

    # extend trajectories & note which 
    for (key, subcomp) in pairs(c.subcomponents)
        new_traj = try
            extend_trajectory!(subcomp, s.substates[key], trajectories[key])
        catch e
            @error "Error when extending traj for  $key" exception=(e, catch_backtrace())
        end
        trajectories[key] = new_traj
        tl = trajectory_length(new_traj)
        if tl < mintime
            mintime = tl
            at_min_time = [key]
        elseif tl == mintime
            push!(at_min_time, key)
        end
    end

    outputted_spikes = Iterators.filter(
        ((compname, outname),) -> !isnothing(outname) && does_output(c, CompOut(compname, outname)),
        ((name, next_spike(c.subcomponents[name], trajectories[name])) for name in at_min_time)
    )
    has_output_spike = !isempty(outputted_spikes)

    outputname = if has_output_spike
        first(r for s in outputted_spikes for r in receivers(c, CompOut(s...)) if r isa Output).id
    else
        nothing
    end

    return CompositeTrajectory(trajectories, mintime, has_output_spike, outputname)
end

function nest_callback(f, nest_at)
    function nested(itr, dt...)
        nested_itr = (
            (   isnothing(name) ? nest_at : (nest_at => name),
                event
            )
            for (name, event) in itr
        )
        f(nested_itr, dt...)
    end
end

# returns a collection of the same top-level type mapping `name => name`
names(t::Tuple) = Tuple(1:length(t))
names(n::NamedTuple) = (;(k=>k for k in keys(n))...)

# by the invariants, this:
# (1) does not extend time past where we have extended the trajectories to, and
# (2) does not extend time past the first spike which occurs in this component
function advance_time_by!(c::CompositeComponent, s::CompositeState, t::CompositeTrajectory, ΔT, f::Function)
    initial_trajectory_length = trajectory_length(t)
    had_next_spike = t.has_next_spike
    initial_next_spike_name = t.next_spike_name
    @assert (initial_trajectory_length >= ΔT) "Should not advance a time past trajectory length! (Trajlen was $initial_trajectory_length; ΔT = $ΔT."

    states = s.substates
    trajectories = t.subtrajectories
    spikes_to_process = []
    for (key, subcomponent) in pairs(c.subcomponents)
        states[key], trajectories[key], outspikes = advance_time_by!(subcomponent, states[key], trajectories[key], ΔT, nest_callback(f, key))
        for outname in outspikes
            push!(spikes_to_process, CompOut(key, outname))
        end
    end

    # send internal spikes to correct components; note the output spikes
    handled_spikes = !isempty(spikes_to_process)
    (outspikes, _) = process_internal_spiking!(c, states, trajectories, spikes_to_process, itr -> f(itr, ΔT))

    # callback
    # f(Iterators.flatten((
    #     ((nothing, OutputSpike(name)) for name in outspikes),
    #     s == new_state ? () : ((nothing, StateChange(new_state)),)
    # )), ΔT)
    f(((nothing, OutputSpike(name)) for name in outspikes), ΔT)

    return (
        CompositeState(states),
        CompositeTrajectory(
            trajectories,
            initial_trajectory_length - ΔT,
            handled_spikes && had_next_spike,
            handled_spikes ? nothing : initial_next_spike_name,
        ),
        outspikes
    )
end
advance_time_by!(c::CompositeComponent, s::CompositeState, t::CompositeTrajectory, ΔT) = advance_time_by!(c, s, t, ΔT, (_,_)->nothing)

function receive_input_spike!(c::CompositeComponent, s::CompositeState, t::CompositeTrajectory, inname, f::Function)
    states, trajectories = s.substates, t.subtrajectories
    (outspikes, state_changed) = process_internal_spiking!(c, states, trajectories, (Input(inname),), f)

    # callback to note new events
    f(Iterators.flatten((
        ((nothing, InputSpike(inname)),),
        ((nothing, OutputSpike(name)) for name in outspikes) #,
        # state_changed ? ((nothing, StateChange(new_state)),) : ()
    )))

    return (
        CompositeState(states),
        CompositeTrajectory(trajectories, 0.0, false, nothing),
        outspikes
    )
end
receive_input_spike!(c::CompositeComponent, s::CompositeState, t::CompositeTrajectory, inname) =
    receive_input_spike!(c, s, t, inname, (_, _) -> nothing)

# TODO: document the interface for these functions
function process_internal_spiking!(
    c::CompositeComponent, s, t,
    initial_spikes, f
)
    spike_queue = Queue{NodeName}()
    for spike in initial_spikes
        enqueue!(spike_queue, spike)
    end

    outspikes = []
    state_changed = false

    while !isempty(spike_queue)
        spike = dequeue!(spike_queue)
        for receiver in receivers(c, spike)
            enqueue!(spike_queue, receiver)
        end

        sc = handle_spike!(c, s, t, spike, outspikes, spike_queue, f)
        state_changed = state_changed || sc
    end

    return (outspikes, state_changed)
end

function handle_spike!(_, _, _, r::Output, outspikes, _, _)
    push!(outspikes, r.id)
    return false
end
function handle_spike!(c, s, t, receiver::CompIn, outspikes, spike_queue, f)
    cn = receiver.comp_name
    oldstate = s[cn]
    s[cn], t[cn], out_spike_names = receive_input_spike!(
        c.subcomponents[cn], s[cn], t[cn], receiver.in_name,
        nest_callback(f, cn)
    )
    for outname in out_spike_names
        enqueue!(spike_queue, CompOut(cn, outname))
    end

    return oldstate == s[cn]
end
handle_spike!(_, _, _, r::Union{Input, CompOut}, _, _, _) = false

include("utils.jl")
export flatten, to_indexed_vals

end # module
