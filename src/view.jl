include("bilog.jl")
include("colorize_field.jl")
using ImageView
using Parameters

function imshowall(images, windowed_energy_map, integrated_energy_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params)
    @unpack dbrange, orientation = sim_params

    if orientation == :vertical
        images = mapslices(rotr90, images; dims=[1, 2])
        windowed_energy_map = rotr90(windowed_energy_map)
        integrated_energy_map = rotr90(integrated_energy_map)
        transmit_time_map = rotr90(transmit_time_map)
        peak_to_peak_time_delta_map = rotr90(peak_to_peak_time_delta_map)
    end

    ImageView.imshow(colorize_field(bilog(images, dbrange/2)))
    ImageView.imshow(bilog(windowed_energy_map, dbrange))
    ImageView.imshow(bilog(integrated_energy_map, dbrange))
    ImageView.imshow(transmit_time_map)
    ImageView.imshow(peak_to_peak_time_delta_map)

    return  # nothing
end
