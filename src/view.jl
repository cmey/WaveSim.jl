include("bilog.jl")
include("colorize_field.jl")
using ImageView
using Parameters


function imshowall(images, beam_energy_map, transmit_time_map, sim_params)
    @unpack dbrange = sim_params
    imshow(colorize_field(bilog(images, dbrange)))
    imshow(bilog(beam_energy_map, dbrange))
    imshow(transmit_time_map)

    return  # nothing
end
