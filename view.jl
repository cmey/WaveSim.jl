include("bilog.jl")
using ImageView
using Colors

function colorize_field(images)
    num_colormap_entries = 100
    cmap = colormap("RdBu", num_colormap_entries)

    images_bilog = bilog(images)
    min, max = extrema(images_bilog)
    images_scaled_for_indexing = clamp.((images_bilog - min) / (max - min) * num_colormap_entries + 1, 1, num_colormap_entries)
    images_indexed = Int.(round.(images_scaled_for_indexing))

    return cmap[images_indexed]
end

function imshow4d(images)
    imshow(colorize_field(images), axes=(2, 3))
end

