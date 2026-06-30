using CairoMakie
using Interpolations
# using JLD
using Makie.Colors
using Parameters
using Base.Threads

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

function save_sim_heatmap(data, extent, xlabel, ylabel, title, colormap, dbrange, output_path, filename; normlog_enabled=true, red_green_red=false)
    centers_x = range(extent[1], extent[2], length=size(data, 1))
    centers_y = range(extent[3], extent[4], length=size(data, 2))
    
    fig = Figure()
    ax = Axis(fig[1, 1], width=size(data, 1), height=size(data, 2), xlabel=xlabel, ylabel=ylabel, title=title)
    
    plot_data = normlog_enabled ? normlog(data, dbrange) : data
    
    cmap = colormap
    if red_green_red
        cmap = cgrad([:red, :green, :red], 256, categorical = false)
    end
    
    hm = heatmap!(ax, centers_x, centers_y, plot_data, colormap=cmap)
    Colorbar(fig[:, end+1], hm)
    resize_to_layout!(fig)
    save(joinpath(output_path, filename), fig)
end

function saveall(images, windowed_energy_map, integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params, output_path="images", conditions_string="")
    println("Saving simulation results to: ", output_path)
    mkpath(output_path)

    # save(joinpath(output_path, make_filename("saved_data.jld", conditions_string)), "images", images; compress=true)

    @unpack tx_frequency, pulse_cycles, fov, dbrange, orientation, beamplot_axes = sim_params
    lateral_label = slice_lateral_label(beamplot_axes)
    slice_name = string(beamplot_axes)

    vertical_enabled = orientation == :vertical
    extent = [0, fov[2], -1/2 * fov[1], 1/2 * fov[1]]
    xlabel = "Depth [m]"
    ylabel = "Azimuth [m]"

    if vertical_enabled
        extent = [-1/2 * fov[1], 1/2 * fov[1], fov[2], 0]
        xlabel, ylabel = ylabel, xlabel
    end

    # Define transformation helper
    transform_map(M) = vertical_enabled ? rotr90(M)' : M'

    @sync begin
        # Windowed energy map
        Threads.@spawn save_sim_heatmap(transform_map(windowed_energy_map), extent, xlabel, ylabel, "Windowed energy map [dB]\n$(conditions_string)", :jet1, dbrange, output_path, make_filename("$(slice_name)_windowed_energy_map.png", conditions_string))

        # Integrated energy map
        Threads.@spawn save_sim_heatmap(transform_map(integrated_energy_map), extent, xlabel, ylabel, "Time-integrated energy map [dB]\n$(conditions_string)", :jet1, dbrange, output_path, make_filename("$(slice_name)_integrated_energy_map.png", conditions_string))

        # Peak-to-peak map
        Threads.@spawn save_sim_heatmap(transform_map(peak_to_peak_map), extent, xlabel, ylabel, "Peak-to-peak map [dB]\n$(conditions_string)", :jet1, dbrange, output_path, make_filename("$(slice_name)_peak_to_peak_map.png", conditions_string))

        # Transmit time map
        Threads.@spawn save_sim_heatmap(transform_map(transmit_time_map) .* 1e6, extent, xlabel, ylabel, "Transmit time map [µs]\n$(conditions_string)", :jet1, dbrange, output_path, make_filename("$(slice_name)_transmit_time_map.png", conditions_string); normlog_enabled=false)

        # Peak-to-peak time-delta map
        Threads.@spawn begin
            max_display_peak_to_peak_duration = 1.0f0 / tx_frequency * pulse_cycles
            data_delta = clamp.(transform_map(peak_to_peak_time_delta_map), -max_display_peak_to_peak_duration, max_display_peak_to_peak_duration)
            save_sim_heatmap(data_delta .* 1e6, extent, xlabel, ylabel, "Peak-to-peak time-delta map [µs]\n$(conditions_string)", :jet1, dbrange, output_path, make_filename("$(slice_name)_peak_to_peak_time_delta_map.png", conditions_string); normlog_enabled=false, red_green_red=true)
        end

        # Line scan at end depth
        Threads.@spawn begin
            p2p_transformed = transform_map(peak_to_peak_map)
            line_scan = p2p_transformed[end, :]
            fig = Figure()
            ax = Axis(fig[1, 1], xlabel="$(lateral_label) [m]", ylabel="Peak-to-peak amplitude [AU]", title="Linear scan - $(lateral_label) at $(fov[2]) [m] depth\n$(conditions_string)")
            centers_lat = range(extent[3], extent[4], length=length(line_scan))
            scatterlines!(ax, centers_lat, line_scan)
            resize_to_layout!(fig)
            save(joinpath(output_path, make_filename("$(slice_name)_line_scan.png", conditions_string)), fig)
        end

        # Arc scan passing through end depth
        Threads.@spawn begin
            p2p_transformed = transform_map(peak_to_peak_map)
            center = (1, size(p2p_transformed, 2) / 2)
            radius = min(size(p2p_transformed, 1), size(p2p_transformed, 2)) / 2 - 1
            num_points = 200
            angles = range(-π/2, π/2, length=num_points)
            x_arc = floor.([center[1] + radius * cos(θ) for θ in angles])
            y_arc = floor.([center[2] + radius * sin(θ) for θ in angles])
            itp = interpolate(p2p_transformed, BSpline(Linear()))
            arc_values = [itp(x, y) for (x, y) in zip(x_arc, y_arc)]
            fig = Figure()
            ax = Axis(fig[1, 1], xlabel="Angle [degs]", ylabel="Peak-to-peak amplitude [AU]", title="Arc scan - $(lateral_label) at $(fov[2]) [m] depth\n$(conditions_string)")
            scatterlines!(ax, angles / π * 180, arc_values)
            resize_to_layout!(fig)
            save(joinpath(output_path, make_filename("$(slice_name)_arc_scan.png", conditions_string)), fig)
        end

        # Peak-to-peak map marked with arc
        Threads.@spawn begin
            p2p_transformed = transform_map(peak_to_peak_map)
            center = (1, size(p2p_transformed, 2) / 2)
            radius = min(size(p2p_transformed, 1), size(p2p_transformed, 2)) / 2 - 1
            num_points = 200
            angles = range(-π/2, π/2, length=num_points)
            x_arc = floor.([center[1] + radius * cos(θ) for θ in angles])
            y_arc = floor.([center[2] + radius * sin(θ) for θ in angles])

            fig = Figure()
            ax = Axis(fig[1, 1], width=size(p2p_transformed, 1), height=size(p2p_transformed, 2), xlabel=xlabel, ylabel=ylabel, title="Peak-to-peak map [dB]\n$(conditions_string)")
            centers_x = range(extent[1], extent[2], length=size(p2p_transformed, 1))
            centers_y = range(extent[3], extent[4], length=size(p2p_transformed, 2))
            img = normlog(p2p_transformed, dbrange)
            for (x, y) in zip(x_arc, y_arc)
                img[round(Int, x), round(Int, y)] = dbrange
            end
            hm = heatmap!(ax, centers_x, centers_y, img)
            Colorbar(fig[:, end+1], hm)
            resize_to_layout!(fig)
            save(joinpath(output_path, make_filename("$(slice_name)_peak_to_peak_map_with_arc.png", conditions_string)), fig)
        end

        # Wave propagation movie (now parallelized with the rest)
        Threads.@spawn begin
            println("Generating wave propagation movie...")
            max_val = maximum(abs, images)
            element_polygons = orient_beamplot_polygons(WaveSim.beamplot_element_polygons(sim_params), vertical_enabled)
            
            # We'll use a single frame to get vmin/vmax for colorscale
            sample_frame = transform_map(images[:, :, end÷2])
            vmin, vmax = extrema(bilog(sample_frame, max_val, dbrange/2))

            fig_movie = Figure()
            size_x, size_y = vertical_enabled ? (size(images, 1), size(images, 2)) : (size(images, 2), size(images, 1))
            ax_movie = Axis(fig_movie[1, 1], width=size_x, height=size_y, xlabel=xlabel, ylabel=ylabel, title="Wave amplitude [dB]\n$(conditions_string)")
            centers_mx = range(extent[1], extent[2], length=size_x)
            centers_my = range(extent[3], extent[4], length=size_y)

            data_obs = Observable(bilog(transform_map(images[:, :, 1]), max_val, dbrange/2))
            hm_movie = heatmap!(ax_movie, centers_mx, centers_my, data_obs; colormap=:RdBu_11, colorrange=(vmin, vmax))
            for polygon in element_polygons
                poly!(ax_movie, Point2f[(corner[1], corner[2]) for corner in polygon], color=(:white, 0.25), strokecolor=:black, strokewidth=1.0)
            end
            Colorbar(fig_movie[:, end+1], hm_movie)
            resize_to_layout!(fig_movie)
            
            framerate = 30
            duration = 3
            time_indices = round.(Int, range(1, size(images, 3), length=framerate*duration))
            record(fig_movie, joinpath(output_path, make_filename("$(slice_name)_wave_propagation.gif", conditions_string)), time_indices; framerate=framerate) do time_index
                data_obs[] = bilog(transform_map(images[:, :, time_index]), max_val, dbrange/2)
                if time_index % 10 == 0 || time_index == time_indices[end]
                    println("Saving simulation frame ", time_index, "/", size(images, 3))
                end
            end
        end
    end

    println("All simulation results saved to ", output_path)
    return
end
