module TestSpikingCircuits

using Test
using SpikingCircuits
import Circuits
using Circuits: flatten, PrimitiveComponent, 
                CompositeValue, CompIn, CompOut,
                Input, Output, CompositeComponent,
                PrimitiveValue, IndexedValues

include("flatten.jl")

end # module
