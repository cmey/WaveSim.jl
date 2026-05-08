using CairoMakie
using Interpolations
# using JLD
using Makie.Colors
using Parameters

function make_filename(original_filename, conditions_string)
    if isempty(conditions_string)
        return original_filename
    else
        return "$(conditions_string)_$(original_filename)"
    end
end

function slice_lateral_label(beamplot_axes::Symbol)
  # Human-readable label for the lateral axis of the selected slice.
  return beamplot_axes == :elevation_depth ? "Elevation" : "Azimuth"
end

function orient_beamplot_polygons(element_polygons, vertical_enabled::Bool)
    if !vertical_enabled
        return element_polygons
    end

    return [
        ntuple(i -> SVector{2, Float32}(polygon[i][2], polygon[i][1]), length(polygon))
        for polygon in element_polygons
    ]
end

function saveall(images, windowed_energy_map, integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params, output_path="images", conditions_string="")
    println("Saving simulation results to: ", output_path)
    mkpath(output_path)

    # save(joinpath(output_path, make_filename("saved_data.jld", conditions_string)), "images", images; compress=true)

    @unpack tx_frequency, pulse_cycles, fov, dbrange, orientation, beamplot_axes = sim_params
    lateral_label = slice_lateral_label(beamplot_axes)
    slice_name = string(beamplot_axes)

    extent = [0, fov[2], -1/2 * fov[1], 1/2 * fov[1]]
    xlabel = "Depth [m]"
    ylabel = "Azimuth [m]"

    vertical_enabled = orientation == :vertical
    if vertical_enabled
        extent = [-1/2 * fov[1], 1/2 * fov[1], fov[2], 0]
        xlabel, ylabel = ylabel, xlabel

        images = mapslices(rotr90, images; dims=[1, 2])
        windowed_energy_map = rotr90(windowed_energy_map)
        integrated_energy_map = rotr90(integrated_energy_map)
        peak_to_peak_map = rotr90(peak_to_peak_map)
        transmit_time_map = rotr90(transmit_time_map)
        peak_to_peak_time_delta_map = rotr90(peak_to_peak_time_delta_map)
    end

    # Windowed energy map
    windowed_energy_map = windowed_energy_map';
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(windowed_energy_map)[1], height=size(windowed_energy_map)[2], xlabel=xlabel, ylabel=ylabel, title="Windowed energy map [dB]\n$(conditions_string)")
    centers_x = range(extent[1], extent[2], length=size(windowed_energy_map)[1])
    centers_y = range(extent[3], extent[4], length=size(windowed_energy_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, WaveSim.normlog(windowed_energy_map, dbrange), colormap=:jet1)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, make_filename("$(slice_name)_windowed_energy_map.png", conditions_string)), fig)

    # Integrated energy map
    integrated_energy_map = integrated_energy_map';
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(integrated_energy_map)[1], height=size(integrated_energy_map)[2], xlabel=xlabel, ylabel=ylabel, title="Time-integrated energy map [dB]\n$(conditions_string)")
    centers_x = range(extent[1], extent[2], length=size(integrated_energy_map)[1])
    centers_y = range(extent[3], extent[4], length=size(integrated_energy_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, WaveSim.normlog(integrated_energy_map, dbrange), colormap=:jet1)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, make_filename("$(slice_name)_integrated_energy_map.png", conditions_string)), fig)

    # Peak-to-peak map
    peak_to_peak_map = peak_to_peak_map';
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(peak_to_peak_map)[1], height=size(peak_to_peak_map)[2], xlabel=xlabel, ylabel=ylabel, title="Peak-to-peak map [dB]\n$(conditions_string)")
    centers_x = range(extent[1], extent[2], length=size(peak_to_peak_map)[1])
    centers_y = range(extent[3], extent[4], length=size(peak_to_peak_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, WaveSim.normlog(peak_to_peak_map, dbrange), colormap=:jet1)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, make_filename("$(slice_name)_peak_to_peak_map.png", conditions_string)), fig)

    # Transmit time map
    transmit_time_map = transmit_time_map';
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(transmit_time_map)[1], height=size(transmit_time_map)[2], xlabel=xlabel, ylabel=ylabel, title="Transmit time map [µs]\n$(conditions_string)")
    centers_x = range(extent[1], extent[2], length=size(peak_to_peak_map)[1])
    centers_y = range(extent[3], extent[4], length=size(peak_to_peak_map)[2])
    hm = heatmap!(ax, centers_x, centers_y, transmit_time_map .* 1e6, colormap=:jet1)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, make_filename("$(slice_name)_transmit_time_map.png", conditions_string)), fig)

    # Peak-to-peak time-delta map
    peak_to_peak_time_delta_map = peak_to_peak_time_delta_map';
    max_display_peak_to_peak_duration = 1.0f0 / tx_frequency * pulse_cycles
    peak_to_peak_time_delta_map = clamp.(peak_to_peak_time_delta_map, -max_display_peak_to_peak_duration, max_display_peak_to_peak_duration)
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(peak_to_peak_time_delta_map)[1], height=size(peak_to_peak_time_delta_map)[2], xlabel=xlabel, ylabel=ylabel, title="Peak-to-peak time-delta map [µs]\n$(conditions_string)")
    centers_x = range(extent[1], extent[2], length=size(peak_to_peak_time_delta_map)[1])
    centers_y = range(extent[3], extent[4], length=size(peak_to_peak_time_delta_map)[2])
    # Define custom red–green–red colormap
    red_green_red = cgrad([:red, :green, :red], 256, categorical = false)
    hm = heatmap!(ax, centers_x, centers_y, peak_to_peak_time_delta_map .* 1e6, colormap=red_green_red)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, make_filename("$(slice_name)_peak_to_peak_time_delta_map.png", conditions_string)), fig)

    # Line scan at end depth
    line_scan = peak_to_peak_map[end, :]
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="$(lateral_label) [m]", ylabel="Peak-to-peak amplitude [AU]", title="Linear scan - $(lateral_label) at $(fov[2]) [m] depth\n$(conditions_string)")
    centers_x = range(extent[3], extent[4], length=length(line_scan))
    scatterlines!(ax, centers_x, line_scan)
    resize_to_layout!(fig)
    save(joinpath(output_path, make_filename("$(slice_name)_line_scan.png", conditions_string)), fig)

    # Arc scan passing through end depth
    center = (1, size(peak_to_peak_map)[2] / 2)  # (x_center, y_center)
    radius = min(size(peak_to_peak_map)[1], size(peak_to_peak_map)[2]) / 2 - 1
    num_points = 200  # number of points along the curve
    angles = range(-π/2, π/2, length=num_points)
    x_arc = floor.([center[1] + radius * cos(θ) for θ in angles])
    y_arc = floor.([center[2] + radius * sin(θ) for θ in angles])
    itp = interpolate(peak_to_peak_map, BSpline(Linear()))
    arc_values = [itp(x, y) for (x, y) in zip(x_arc, y_arc)]
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Angle [degs]", ylabel="Peak-to-peak amplitude [AU]", title="Arc scan - $(lateral_label) at $(fov[2]) [m] depth\n$(conditions_string)")
    scatterlines!(ax, angles / π * 180, arc_values)
    resize_to_layout!(fig)
    save(joinpath(output_path, make_filename("$(slice_name)_arc_scan.png", conditions_string)), fig)

    # Peak-to-peak map marked with arc
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(peak_to_peak_map)[1], height=size(peak_to_peak_map)[2], xlabel=xlabel, ylabel=ylabel, title="Peak-to-peak map [dB]\n$(conditions_string)")
    centers_x = range(extent[1], extent[2], length=size(peak_to_peak_map)[1])
    centers_y = range(extent[3], extent[4], length=size(peak_to_peak_map)[2])
    image = WaveSim.normlog(peak_to_peak_map, dbrange)
    for (x, y) in zip(x_arc, y_arc)
        image[round(Int, x), round(Int, y)] = dbrange  # mark the arc
    end
    hm = heatmap!(ax, centers_x, centers_y, image)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, make_filename("$(slice_name)_peak_to_peak_map_with_arc.png", conditions_string)), fig)

    # Wave propagation movie
    images = permutedims(images, (2, 1, 3));
    wave_field_images = WaveSim.bilog(images, dbrange/2)
    vmin, vmax = extrema(wave_field_images)
    element_polygons = orient_beamplot_polygons(WaveSim.beamplot_element_polygons(sim_params), vertical_enabled)
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(wave_field_images)[1], height=size(wave_field_images)[2], xlabel=xlabel, ylabel=ylabel, title="Wave amplitude [dB]\n$(conditions_string)")
    centers_x = range(extent[1], extent[2], length=size(wave_field_images)[1])
    centers_y = range(extent[3], extent[4], length=size(wave_field_images)[2])

    # Use an Observable for the data to update the plot efficiently
    data_obs = Observable(wave_field_images[:, :, 1])
    hm = heatmap!(ax, centers_x, centers_y, data_obs; colormap=:RdBu_11, colorrange=(vmin, vmax))
    for polygon in element_polygons
        poly!(ax, Point2f[(corner[1], corner[2]) for corner in polygon], color=(:white, 0.25), strokecolor=:black, strokewidth=1.0)
    end
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    framerate = 30  # [fps]
    duration = 3  # [s]
    time_indices = round.(Int, range(1, size(wave_field_images)[3], length=framerate*duration))
    record(fig, joinpath(output_path, make_filename("$(slice_name)_wave_propagation.gif", conditions_string)), time_indices; framerate=framerate) do time_index
        data_obs[] = wave_field_images[:, :, time_index]
        # Skip printing every frame for 1000+ frames to avoid console overhead
        if time_index % 10 == 0 || time_index == time_indices[end]
            println("Saving simulation frame ", time_index, "/", size(wave_field_images)[3])
        end
    end

    return  # nothing
end
