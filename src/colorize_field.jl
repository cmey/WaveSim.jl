using Makie.Colors


function colorize_field(images)
    num_colormap_entries = 100
    cmap = colormap("RdBu", num_colormap_entries)

    min, max = extrema(images)
    images_scaled_for_indexing = clamp.((images .- min) / (max - min) * num_colormap_entries .+ 1, 1, num_colormap_entries)
    images_indexed = Int.(round.(images_scaled_for_indexing))

    return cmap[images_indexed]
end
