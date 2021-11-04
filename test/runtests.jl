module TestSpikingCircuits

using Test
using SpikingCircuits
import Circuits
using Circuits: inline, PrimitiveComponent, 
                CompositeValue, CompIn, CompOut,
                Input, Output, CompositeComponent,
                PrimitiveValue, IndexedValues

include("inline.jl")

end # module
