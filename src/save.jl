using CairoMakie
using Makie.Colors
using Parameters
include("bilog.jl")

const wave_propagation_filename = "wave_propagation.gif"
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

    # Wave propagation movie
    images = permutedims(images, (2, 1, 3));
    wave_field_images = bilog(images, dbrange)
    vmin, vmax = extrema(wave_field_images)
    fig = Figure(size=size(wave_field_images[:, :, 1]))
    ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel=xlabel, ylabel=ylabel, title="Wave amplitude [dB]")
    centers_x = range(extent[1], extent[2], length=size(beam_energy_map)[1])
    centers_y = range(extent[3], extent[4], length=size(beam_energy_map)[2])
    framerate = 30  # [fps]
    duration = 3  # [s]
    time_indices = Int.(round.(range(1, size(wave_field_images)[3], length=framerate*duration)))
    record(fig, joinpath(output_path, wave_propagation_filename), time_indices; framerate=framerate) do time_index
        println("Saving simulation frame ", time_index, "/", size(wave_field_images)[3], " ")
        hm = heatmap!(ax, centers_x, centers_y, wave_field_images[:, :, time_index]; colormap=:RdBu_11, colorrange=(vmin, vmax))
        if time_index == 1
            Colorbar(fig[:, end+1], hm)
        end
    end

    # Beam energy map
    beam_energy_map = beam_energy_map';
    fig = Figure(size = size(beam_energy_map))
    ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel=xlabel, ylabel=ylabel, title="Beam energy map [dB]")
    centers_x = range(extent[1], extent[2], length=size(beam_energy_map)[1])
    centers_y = range(extent[3], extent[4], length=size(beam_energy_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, bilog(beam_energy_map, dbrange))
    Colorbar(fig[:, end+1], hm)
    save(joinpath(output_path, beam_energy_map_filename), fig)

    # Transmit time map
    transmit_time_map = transmit_time_map';
    fig = Figure(size = size(transmit_time_map))
    ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel=xlabel, ylabel=ylabel, title="Transmit time map [µs]")
    centers_x = range(extent[1], extent[2], length=size(beam_energy_map)[1])
    centers_y = range(extent[3], extent[4], length=size(beam_energy_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, transmit_time_map .* 1e6)
    Colorbar(fig[:, end+1], hm)
    save(joinpath(output_path, transmit_time_map_filename), fig)

    return  # nothing
end
