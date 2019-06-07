using Parameters
using PyCall
using PyPlot
include("bilog.jl")
include("colorize_field.jl")

anim = pyimport("matplotlib.animation")  # After using PyPlot, so matplotlib installed via conda.
const wave_propagation_filename = "wave_propagation.mp4"
const beam_energy_map_filename = "beam_energy_map.png"
const transmit_time_map_filename = "transmit_time_map.png"

function saveall(images, beam_energy_map, transmit_time_map, sim_params, output_path="images")
    @unpack fov, dbrange, orientation = sim_params

    extent = [0, fov[2], -1/2 * fov[1], 1/2 * fov[1]]
    xlabel = "Depth [m]"
    ylabel = "Azimuth [m]"

    vertical_enabled = orientation == :vertical
    if vertical_enabled
        extent = [-1/2 * fov[1], 1/2 * fov[1], fov[2], 0]
        xlabel, ylabel = ylabel, xlabel

        images = mapslices(rotr90, images; dims=[1, 2])
        beam_energy_map = rotr90(beam_energy_map)
        transmit_time_map = rotr90(transmit_time_map)
    end

    wave_field = bilog(images, dbrange)
    vmin, vmax = extrema(wave_field)

    fig = figure()
    ims = []
    for i_time in 1:size(wave_field, 3)
        im = PyPlot.imshow(wave_field[:, :, i_time], extent=extent, cmap="viridis", vmin=vmin, vmax=vmax)
        if 1 == i_time
            PyPlot.colorbar()
        end
        PyPlot.title("Wave amplitude [dB]")
        PyPlot.xlabel(xlabel)
        PyPlot.ylabel(ylabel)
        push!(ims, PyCall.PyObject[im])
    end
    #= close() =#

    # If matplotlib complains, ensure that
    # a) ffmpeg is installed with libx264 support, and
    # b) matplotlib is built with ffmpeg support enabled
    ani = anim.ArtistAnimation(fig, ims, interval=30, blit=true, repeat=false)
    anim_filename = wave_propagation_filename
    ffmpeg_writer = anim.FFMpegWriter()
    ani.save(joinpath(output_path, anim_filename), ffmpeg_writer)
    #= ani.save("ani.gif"); =#

    PyPlot.figure()
    PyPlot.imshow(bilog(beam_energy_map, dbrange), extent=extent)
    PyPlot.colorbar()
    PyPlot.title("Beam energy map [dB]")
    PyPlot.xlabel(xlabel)
    PyPlot.ylabel(ylabel)
    PyPlot.savefig(joinpath(output_path, beam_energy_map_filename))
    PyPlot.close()

    PyPlot.figure()
    PyPlot.imshow(transmit_time_map .* 1e6, extent=extent)
    PyPlot.colorbar()
    PyPlot.title("Transmit time map [µs]")
    PyPlot.xlabel(xlabel)
    PyPlot.ylabel(ylabel)
    PyPlot.savefig(joinpath(output_path, transmit_time_map_filename))
    PyPlot.close()

    return  # nothing
end
