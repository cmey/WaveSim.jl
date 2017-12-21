using PyPlot
include("bilog.jl")
include("colorize_field.jl")

const wave_propagation_filename = "wave_propagation.png"
const beam_energy_map_filename = "beam_energy_map.png"
const transmit_time_map_filename = "transmit_time_map.png"


function saveall(sim_params, images, beam_energy_map, transmit_time_map, output_path, dbrange=40)
    fov = sim_params["fov"]
    extent=[0, fov[2], -1/2 * fov[1], 1/2 * fov[1]]

    #= colorized_field = colorize_field(bilog(images, dbrange)) =#
    #= figure() =#
    #= for i_time in size(colorized_field)[3] =#
    #=     imshow(colorized_field[:, :, i_time], extent=extent) =#
    #=     colorbar() =#
    #=     title("wave amplitude [dB]") =#
    #=     xlabel("Depth [m]") =#
    #=     ylabel("Azimuth [m]") =#

    #=     img_filename = split(wave_propagation_filename, '.') =#
    #=     insert!(img_filename, 2, ".") =#
    #=     insert!(img_filename, 2, string(i_time)) =#
    #=     img_filename = join(img_filename) =#

    #=     savefig(joinpath(output_path, img_filename)) =#
    #= end =#
    #= close() =#

    figure()
    imshow(bilog(beam_energy_map, dbrange), extent=extent)
    colorbar()
    title("beam energy map [dB]")
    xlabel("Depth [m]")
    ylabel("Azimuth [m]")
    savefig(joinpath(output_path, beam_energy_map_filename))
    close()

    figure()
    imshow(transmit_time_map .* 1e6, extent=extent)
    colorbar()
    title("transmit time map [µs]")
    xlabel("Depth [m]")
    ylabel("Azimuth [m]")
    savefig(joinpath(output_path, transmit_time_map_filename))
    close()
end
