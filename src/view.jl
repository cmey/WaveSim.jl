include("bilog.jl")
include("colorize_field.jl")
using ImageView
using Parameters

function imshowall(images, beam_energy_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params)
    @unpack dbrange, orientation = sim_params

    if orientation == :vertical
        images = mapslices(rotr90, images; dims=[1, 2])
        beam_energy_map = rotr90(beam_energy_map)
        transmit_time_map = rotr90(transmit_time_map)
        peak_to_peak_time_delta_map = rotr90(peak_to_peak_time_delta_map)
    end

    ImageView.imshow(colorize_field(bilog(images, dbrange)))
    ImageView.imshow(bilog(beam_energy_map, dbrange))
    ImageView.imshow(transmit_time_map)
    ImageView.imshow(peak_to_peak_time_delta_map)

    return  # nothing
end
