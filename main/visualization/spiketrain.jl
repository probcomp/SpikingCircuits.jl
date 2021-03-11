module SpiketrainViz
using GLMakie, AbstractPlotting, Colors
AbstractPlotting.inline!(false)

export draw_spiketrain_figure, get_spiketrain_figure

"""
    draw_spiketrain_figure(
        spiketrains; # Vector of `Vector{Float64}`s
        names=String[], # Vector of neuron names for the first `length(names)` spiketrains
        colors=Color[], # Vector of colors for the first `length(colors)` spiketrains
        resolution=(1280, 720),
        figure_title="Spiketrain",
        time=0.
    )
"""
function draw_spiketrain_figure(args...; kwargs...)
    (f, colors) = get_spiketrain_figure(args...; kwargs...)
    display(f)

    return (f, colors)
end

function get_spiketrain_figure(
    spiketrains; # Vector of `Vector{Float64}`s
    names=String[], # Vector of neuron names for the first `length(names)` spiketrains
    colors=Color[], # Vector of colors for the first `length(colors)` spiketrains
    resolution=(1280, 720),
    figure_title="Spiketrain",
    time=0.
)
    f = Figure(;resolution)
    ax = f[1, 1] = Axis(f; title = figure_title)

    draw_spiketrain!(ax, spiketrains, names, colors, time)

    return (f, colors)
end

function draw_spiketrain!(ax, spiketrains, names, colors, time)
    hideydecorations!(ax, ticklabels=false)

    # set neuron names on axis label
    ypositions = 1:length(spiketrains)
    trainheight = 1

    colors = collect(Iterators.flatten((colors, Iterators.repeated(RGB(0, 0, 0), length(names) - length(colors)))))
    @assert length(names) == length(colors)

    for (spiketrain, pos, color) in zip(spiketrains, ypositions, colors)        
        draw_single_spiketrain!(ax, spiketrain, pos, trainheight, time, color)
    end

    ylims!(ax, (first(ypositions) - 1, last(ypositions) + 1))
    println("Set ylims to $((first(ypositions) - 1, last(ypositions) + 1))")
    ax.yticks = (ypositions[1:length(names)], names)
    ax.yticklabelcolor[] = colors

    return colors
end

function draw_single_spiketrain!(ax, spiketimes::Vector{Float64}, ypos, height, current_time, color=RGB(0, 0, 0))
    y1 = ypos - height/2; y2 = ypos + height/2
    times = vcat([
        [Point2f0(t - current_time, y1), Point2f0(t - current_time, y2)]
        for t in spiketimes
    ]...)

    if !isempty(times)
        linesegments!(ax, times; color)
    end
end

end