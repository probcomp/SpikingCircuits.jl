"""
    InputFunctionPoisson <: PrimitiveComponent{Spiking}
    InputFunctionPoisson(
        input_functions::Vector{Function},
        memories::Vector{Float64},
        rate_fn::Function
    )

A Poisson neuron whose potential is
``u(t) = ∑_i{F_i(N_i(t))}``,
where ``F_i`` denotes `input_functions[i]`,
and `N_i(t)` is the number of spikes received by the `i`th input
in the past `memories[i]` time units.
The rate at time `t` is `rate_fn(u(t))`.
"""
struct InputFunctionPoisson <: PrimitiveComponent{Spiking}
    input_functions::Tuple{Vararg{Function}}
    memories::Tuple{Vararg{Float64}}
    rate_fn::Function
end
Circuits.inputs(p::InputFunctionPoisson) = IndexedValues(SpikeWire() for _ in p.input_functions)
Circuits.outputs(::InputFunctionPoisson) = CompositeValue((out=SpikeWire(),))

Base.:(==)(a::InputFunctionPoisson, b::InputFunctionPoisson) = a.input_functions == b.input_functions && a.memories == b.memories && a.rate_fn == b.rate_fn
Base.hash(a::InputFunctionPoisson, h::UInt) = hash(a.input_functions, hash(a.memories, hash(a.rate_fn, h)))

### simulation methods ###
"""
    InputTimesState <: SpikingSimulator.State

State storing how long ago each input spike was received for each input.
For `i::InputTimesState`, `i.input_times[j]` is a vector of each value `t`
such that a spike was received for input `j` `t` time units ago,
and this spike is still remembered.
"""
struct InputTimesState <: Sim.State
    input_times::Vector{Vector{Float64}}
end
# how long until we one of the currently remembered spikes stops being remembered?
time_to_removal(i::InputTimesState, memories) =
    minimum(
        collect(memories) - map(x -> maximum(x; init=-Inf), i.input_times); 
        init=Inf
    )
receive_spike(st::InputTimesState, inidx) =
    InputTimesState([
        i == inidx ? push!(copy(times), 0.) : times
        for (i, times) in enumerate(st.input_times)
    ])

advance(st::InputTimesState, ΔT, memories) =
    InputTimesState([
        filter(<(memories[i]), times .+ ΔT)
        for (i, times) in enumerate(st.input_times)
    ])

rate(p::InputFunctionPoisson, st::InputTimesState) =
    p.rate_fn(sum(
        in_fn(length(times))
        for (in_fn, times) in zip(p.input_functions, st.input_times)
    ))

Sim.initial_state(p::InputFunctionPoisson) = InputTimesState(
    [Float64[] for _ in p.input_functions]
)
Sim.next_spike(::InputFunctionPoisson, ::Sim.NextSpikeTrajectory) = :out

_sample_spike(p::InputFunctionPoisson, st::InputTimesState) =
    let λ = rate(p, st)
        λ == Inf ? 0. : exponential(λ)
    end

function Sim.extend_trajectory(p::InputFunctionPoisson, st::InputTimesState, ::Sim.EmptyTrajectory)
    time_to_outspike = _sample_spike(p, st)
    elapsed_time = 0.
    while time_to_removal(st, p.memories) < time_to_outspike + elapsed_time
        elapsed_time += time_to_removal(st, p.memories)
        st = advance(st, time_to_removal(st, p.memories), p.memories)
        time_to_outspike = _sample_spike(p, st)
    end
    Sim.NextSpikeTrajectory(time_to_outspike + elapsed_time)
end

Sim.extend_trajectory(::InputFunctionPoisson, ::InputTimesState, t::Sim.NextSpikeTrajectory) = t

Sim.advance_time_by(p::InputFunctionPoisson, st::InputTimesState, t::Sim.NextSpikeTrajectory, ΔT) =
    (
        advance(st, ΔT, p.memories),
        (let remaining_time = t.time_to_next_spike - ΔT
            if remaining_time > 0.
                (Sim.NextSpikeTrajectory(remaining_time), ())
            elseif remaining_time == 0
                (Sim.EmptyTrajectory(), (:out,))
            else
                Sim.advancing_too_far_error()
            end
        end)...
    )

Sim.receive_input_spike(p::InputFunctionPoisson, st::InputTimesState, ::Sim.Trajectory, inidx) =
    let new_st = receive_spike(st, inidx)
        if rate(p, new_st) == Inf
            prev_rate_msg = rate(p, st) == Inf ? " (In this case, the rate was Inf before this spike was received as well.) " : ""
            @warn("InputFunctionPoisson just recieved an input spike; after this, the rate is `Inf`. $prev_rate_msg This can cause an infinite number of spikes in a single instant.")
            (new_st, EmptyTrajectory(), (:out,))
        else
            (new_st, EmptyTrajectory(), ())
        end
    end
