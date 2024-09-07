using CairoMakie
using Makie.Colors
using Parameters
include("bilog.jl")

const wave_propagation_filename = "wave_propagation.gif"
const beam_energy_map_filename = "beam_energy_map.png"
const transmit_time_map_filename = "transmit_time_map.png"
const peak_to_peak_time_delta_map_filename = "peak_to_peak_time_delta_map.png"
const line_scan_filename = "line_scan.png"

function saveall(images, beam_energy_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params, output_path="images")
    mkpath(output_path)

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
        peak_to_peak_time_delta_map = rotr90(peak_to_peak_time_delta_map)
    end

    # Wave propagation movie
    images = permutedims(images, (2, 1, 3));
    wave_field_images = bilog(images, dbrange)
    vmin, vmax = extrema(wave_field_images)
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(wave_field_images)[1], height=size(wave_field_images)[2], xlabel=xlabel, ylabel=ylabel, title="Wave amplitude [dB]")
    centers_x = range(extent[1], extent[2], length=size(wave_field_images)[1])
    centers_y = range(extent[3], extent[4], length=size(wave_field_images)[2])
    hm = heatmap!(ax, centers_x, centers_y, wave_field_images[:, :, 1]; colormap=:RdBu_11, colorrange=(vmin, vmax))
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    framerate = 30  # [fps]
    duration = 3  # [s]
    time_indices = Int.(round.(range(1, size(wave_field_images)[3], length=framerate*duration)))
    record(fig, joinpath(output_path, wave_propagation_filename), time_indices; framerate=framerate) do time_index
        heatmap!(ax, centers_x, centers_y, wave_field_images[:, :, time_index]; colormap=:RdBu_11, colorrange=(vmin, vmax))
        println("Saving simulation frame ", time_index, "/", size(wave_field_images)[3], " ")
    end

    # Beam energy map
    beam_energy_map = beam_energy_map';
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(beam_energy_map)[1], height=size(beam_energy_map)[2], xlabel=xlabel, ylabel=ylabel, title="Beam energy map [dB]")
    centers_x = range(extent[1], extent[2], length=size(beam_energy_map)[1])
    centers_y = range(extent[3], extent[4], length=size(beam_energy_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, bilog(beam_energy_map, dbrange))
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, beam_energy_map_filename), fig)

    # Transmit time map
    transmit_time_map = transmit_time_map';
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(transmit_time_map)[1], height=size(transmit_time_map)[2], xlabel=xlabel, ylabel=ylabel, title="Transmit time map [µs]")
    centers_x = range(extent[1], extent[2], length=size(beam_energy_map)[1])
    centers_y = range(extent[3], extent[4], length=size(beam_energy_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, transmit_time_map .* 1e6)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, transmit_time_map_filename), fig)

    # Peak-to-peak time-delta map
    peak_to_peak_time_delta_map = peak_to_peak_time_delta_map';
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(peak_to_peak_time_delta_map)[1], height=size(peak_to_peak_time_delta_map)[2], xlabel=xlabel, ylabel=ylabel, title="Peak-to-peak time-delta map [µs]")
    centers_x = range(extent[1], extent[2], length=size(peak_to_peak_time_delta_map)[1])
    centers_y = range(extent[3], extent[4], length=size(peak_to_peak_time_delta_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, peak_to_peak_time_delta_map .* 1e6)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, peak_to_peak_time_delta_map_filename), fig)

    # Line scan at end depth
    line_scan = beam_energy_map'[:, end]
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Azimuth [m]", ylabel="Peak-to-peak amplitude [AU]", title="Linear scan - azimuth at $(fov[2]) [m] depth")
    centers_x = range(extent[3], extent[4], length=length(line_scan))
    scatterlines!(ax, centers_x, line_scan)
    resize_to_layout!(fig)
    save(joinpath(output_path, line_scan_filename), fig)

    return  # nothing
end
