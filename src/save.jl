using PyCall
using PyPlot
include("bilog.jl")
include("colorize_field.jl")

@pyimport matplotlib.animation as anim
const wave_propagation_filename = "wave_propagation"
const beam_energy_map_filename = "beam_energy_map.png"
const transmit_time_map_filename = "transmit_time_map.png"

function saveall(sim_params, images, beam_energy_map, transmit_time_map, output_path="images", dbrange=40)
    fov = sim_params["fov"]
    extent=[0, fov[2], -1/2 * fov[1], 1/2 * fov[1]]

    wave_field = bilog(images, dbrange)
    vmin, vmax = extrema(wave_field)

    fig = figure()
    ims = []
    for i_time in 1:size(wave_field, 3)
        im = imshow(wave_field[:, :, i_time], extent=extent, cmap="viridis", vmin=vmin, vmax=vmax)
        if 1 == i_time
            colorbar()
        end
        title("Wave amplitude [dB]")
        xlabel("Depth [m]")
        ylabel("Azimuth [m]")
        push!(ims, PyCall.PyObject[im])
    end
    #= close() =#

    # If matplotlib complains, ensure that
    # a) ffmpeg is installed with libx264 support, and
    # b) matplotlib is built with ffmpeg support enabled
    ani = anim.ArtistAnimation(fig, ims, interval=30, blit=true, repeat=false)
    anim_filename = wave_propagation_filename * ".mp4"
    ani[:save](joinpath(output_path, anim_filename),
               extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
    #= ani[:save]("ani.gif"); =#

    figure()
    imshow(bilog(beam_energy_map, dbrange), extent=extent)
    colorbar()
    title("Beam energy map [dB]")
    xlabel("Depth [m]")
    ylabel("Azimuth [m]")
    savefig(joinpath(output_path, beam_energy_map_filename))
    close()

    figure()
    imshow(transmit_time_map .* 1e6, extent=extent)
    colorbar()
    title("Transmit time map [µs]")
    xlabel("Depth [m]")
    ylabel("Azimuth [m]")
    savefig(joinpath(output_path, transmit_time_map_filename))
    close()
end
