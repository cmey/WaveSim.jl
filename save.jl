using PyPlot
include("bilog.jl")

const beam_energy_map_filename = "beam_energy_map.pdf"
const transmit_time_map_filename = "transmit_time_map.pdf"


function saveall(sim_params, images, beam_energy_map, transmit_time_map, output_path)
    #= colorized_field = colorize_field(bilog(images)) =#
    #= figure() =#
    #= for i_time in size(colorize_field)[3] =#
    #=     imshow(colorized_field[:, :, i_time]) =#
    #=     savefig() =#
    #= end =#
    fov = sim_params["fov"]
    extent=[0, fov[2], -1/2 * fov[1], 1/2 * fov[1]]

    figure()
    imshow(bilog(beam_energy_map), extent=extent)
    colorbar()
    title("beam energy map [dB]")
    xlabel("Depth [m]")
    ylabel("Azimuth [m]")
    savefig(output_path * beam_energy_map_filename)
    close()

    figure()
    imshow(transmit_time_map .* 1e6, extent=extent)
    colorbar()
    title("transmit time map [µs]")
    xlabel("Depth [m]")
    ylabel("Azimuth [m]")
    savefig(output_path * transmit_time_map_filename)
    close()
end