# TODO:
# 1. Test case: simple one (no recurrence)
# 2. Test case: simple recurrence.
# 3. Test case: explicit passthrough in subcomponent.
# 4. Test case: passthrough which then recurs.
#         - (goes back in for a finite number of cycles).
#
# TODO: also output some data structure which lets us lookup
# the nested name for a neuron before flattening.

struct PlaceholderNeuron <: PrimitiveComponent{Spiking} end
Circuits.inputs(::PlaceholderNeuron) = IndexedValues(SpikeWire() for _ in 1 : 1)
Circuits.outputs(::PlaceholderNeuron) = CompositeValue((out = SpikeWire(),),)

@testset "flatten - unroll inner component" begin
    inner_comp = CompositeComponent(
        CompositeValue((SpikeWire(), SpikeWire(), SpikeWire())),
        CompositeValue((SpikeWire(), SpikeWire(), SpikeWire())),
            (neuron=PlaceholderNeuron(),),
                (
                Input(1) => Output(1),
                Input(2) => Output(2),
                Input(3) => CompIn(:neuron, 1),
                CompOut(:neuron, :out) => Output(3)
                )
        )
    outer_comp = CompositeComponent(
        CompositeValue((SpikeWire(),)), CompositeValue((SpikeWire(),)),
        (subcomp=inner_comp,),
        (
            Input(1) => CompIn(:subcomp, 1),
            CompOut(:subcomp, 1) => CompIn(:subcomp, 2),
            CompOut(:subcomp, 2) => CompIn(:subcomp, 3),
            CompOut(:subcomp, 3) => Output(1)
        )
    )

    flattened = Circuits.flatten(outer_comp)
    println(Circuits.inputs(flattened))
    println(Circuits.outputs(flattened))
end
