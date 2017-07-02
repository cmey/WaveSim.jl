include("bilog.jl")
using Colors, Gtk.ShortNames, ImageView


function colorize_field(images)
    num_colormap_entries = 100
    cmap = colormap("RdBu", num_colormap_entries)

    min, max = extrema(images)
    images_scaled_for_indexing = clamp.((images - min) / (max - min) * num_colormap_entries + 1, 1, num_colormap_entries)
    images_indexed = Int.(round.(images_scaled_for_indexing))

    return cmap[images_indexed]
end


# Get the beam profile spatial map and transmit time of beam energy, where each pixel indicates the maximum energy that was seen at that place, and at what time.
function beam_energy_map_and_transmit_time_map(images)
    maxval, linindices = findmax(images, 3)

    beam_energy_map = squeeze(maxval, 3)

    transmit_time_map = similar(beam_energy_map)

    for linind in eachindex(linindices)
        x, y, t = ind2sub(images, linindices[linind])
        transmit_time_map[linind] = t
    end

    return beam_energy_map, transmit_time_map
end


function imshowall(images)
    beam_energy_map, transmit_time_map = beam_energy_map_and_transmit_time_map(images)

    imshow(colorize_field(bilog(images)))
    imshow(bilog(beam_energy_map))
    imshow(transmit_time_map)
end
