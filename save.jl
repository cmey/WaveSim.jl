using PyPlot
include("bilog.jl")

const beam_energy_map_filename = "beam_energy_map.pdf"
const transmit_time_map_filename = "transmit_time_map.pdf"


function saveall(images, beam_energy_map, transmit_time_map, output_path)
    #= colorized_field = colorize_field(bilog(images)) =#
    #= figure() =#
    #= for i_time in size(colorize_field)[3] =#
    #=     imshow(colorized_field[:, :, i_time]) =#
    #=     savefig() =#
    #= end =#
    figure()
    imshow(bilog(beam_energy_map))
    colorbar()
    title("beam energy map [dB]")
    xlabel("Depth [m]")
    ylabel("Azimuth [m]")
    savefig(output_path * beam_energy_map_filename)
    close()

    figure()
    imshow(transmit_time_map .* 1e6)
    colorbar()
    title("transmit time map [µs]")
    xlabel("Depth [m]")
    ylabel("Azimuth [m]")
    savefig(output_path * transmit_time_map_filename)
    close()
end
