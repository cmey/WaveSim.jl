include("bilog.jl")
include("colorize_field.jl")
using ImageView


function imshowall(sim_params, images, beam_energy_map, transmit_time_map, dbrange=40)
    imshow(colorize_field(bilog(images, dbrange)))
    imshow(bilog(beam_energy_map, dbrange))
    imshow(transmit_time_map)
end
