# WaveSim.jl, wave propagation simulator.
# Christophe Meyer, 2016-2026
module WaveSim

using Parameters
using ProgressMeter
using StaticArrays
using Statistics

include("bilog.jl")
include("colorize_field.jl")
include("view.jl")
include("save.jl")

export WaveSimParameters
export imshowall
export saveall

@enum ApodizationShape Rect Hann

# Configuration
@with_kw mutable struct WaveSimParameters{pulse_shape_func_T, directivity_func_T}
  # Transmit center frequency.
  tx_frequency::Float32 = 3_000_000  # [Hz]
  # Length of the transmit pulse.
  pulse_cycles::Float32 = 2  # [cycles]
  # Speed of sound.
  c::Float32 = 1540  # [m/s]
  # Time span to simulate.
  end_simulation_time::Float32 = 56.0e-6  # [s] Starts at 0 s.
  # Resolution at which to divide the simulation time span.
  temporal_res::Float32 = 0.1e-6  # [s]
  # Field of view, x then z. x=0 centered on aperture center, z=0 at aperture plane.
  fov::SVector{2, Float32} = @SVector [4e-2, 8e-2]  # [m]
  # Spatial resolution of the simulation, x then z.
  spatial_res::SVector{2, Int} = @SVector [256, 512]  # [pixels]
  # Spacing between physical elements of the transducer array in azimuth and elevation.
  transducer_pitch::Float32 = 208e-6  # [m]
  transducer_pitch_elevation::Float32 = 208e-6  # [m]
  # Active aperture size of the transducer array in azimuth and elevation. Will be centered.
  aperture_size::Float32 = 0.02  # [m] Will simulate at least 1 element.
  aperture_size_elevation::Float32 = 0.0  # [m] Will simulate at least 1 element.
  # Radius of curvature of the aperture in azimuth and elevation. Inf means a flat aperture.
  aperture_radius::Float32 = Inf  # [m]
  aperture_radius_elevation::Float32 = Inf  # [m]
  # Distance from transducer to focus point, Inf and -Inf mean plane wave. 0 is not allowed.
  # Positive values mean focusing, negative values mean diverging.
  focus_depth::Float32 = 0.03  # [m]
  focus_depth_elevation::Float32 = 0.03  # [m]
  # Steering angle in azimuth.
  steer_angle::Float32 = 10.0  # [deg]
  steer_angle_elevation::Float32 = 0.0  # [deg]
  # Which 2D slice is being simulated and rendered, :azimuth_depth or :elevation_depth.
  beamplot_axes::Symbol = :azimuth_depth
  # Shape of the transmit pulse.
  pulse_shape_func::pulse_shape_func_T = phase_divided_by_pi -> cospi(phase_divided_by_pi)
  # Shape of the transmit aperture apodization.
  apodization_shape::ApodizationShape = Rect
  # Display dynamic range.
  dbrange::Float32 = 40  # [dB]
  # Plot orientation: beam is horizontal or vertical.
  orientation::Symbol = :horizontal
  # Attenuation coefficient in dB/cm/MHz.
  attenuation_coefficient::Float32 = 0.0  # [dB/cm/MHz]
  # Directivity function of angle and frequency.
  directivity_func::directivity_func_T = default_directivity
end

@inline function besselj1_approx(x::Float32)::Float32
  x2 = x * x
  return x * (0.5f0 - x2 * (1f0 / 16f0 - x2 * (1f0 / 384f0 - x2 * (1f0 / 18432f0 - x2 * (1f0 / 1474560f0)))))
end

# Concrete default directivity function with fully-typed signature to help inference.
# Uses a GPU-friendly polynomial approximation for J1 that is accurate over the angle range used in beam plots.
function default_directivity(θ::Float32, tx_frequency::Float32, c::Float32, transducer_pitch::Float32)::Float32
  # Account for mechanical crosstalk between elements
  crosstalk_factor = 1.2f0  # Assume 20% crosstalk factor
  # From Umchid 2009 Directivity Pattern Measurement of Ultrasound Transducers:
  # https://www.thaiscience.info/Journals/Article/IABE/10892457.pdf
  element_surface_diameter = transducer_pitch * crosstalk_factor
  a = element_surface_diameter / 2.0f0
  λ = c / tx_frequency
  k = 2.0f0 * π / λ
  ka_sin_θ = k * a * sin(θ)
  if abs(ka_sin_θ) < 1f-5
    return 1.0f0
  end
  return convert(Float32, 2.0f0 * besselj1_approx(ka_sin_θ) / ka_sin_θ)
end

@inline function centered_axis_coordinates(num_elements::Int, pitch::Float32)
  # Return element centers measured from the array midpoint.
  if num_elements <= 1
    return Float32[0.0]
  end
  center_offset = (num_elements + 1) / 2
  return Float32[((ielem - center_offset) * pitch) for ielem in 1:num_elements]
end

@inline function curvature_sagitta(coord::Float32, radius::Float32)
  # Map a lateral coordinate onto a curved surface height.
  if !isfinite(radius)
    return 0.0f0
  end
  return radius - sqrt(max(radius * radius - coord * coord, 0.0f0))
end

@inline function curvature_slope(coord::Float32, radius::Float32)
  # Compute the local surface slope used to build the element normal.
  if !isfinite(radius)
    return 0.0f0
  end
  return coord / sqrt(max(radius * radius - coord * coord, 1.0f-12))
end

@inline function element_normal(x::Float32, y::Float32, radius_x::Float32, radius_y::Float32)
  # Build a unit normal for an element on the curved 2D aperture surface.
  dzdx = curvature_slope(x, radius_x)
  dzdy = curvature_slope(y, radius_y)
  nx = -dzdx
  ny = -dzdy
  nz = 1.0f0
  inv_norm = inv(sqrt(nx * nx + ny * ny + nz * nz))
  return SVector{3,Float32}(nx * inv_norm, ny * inv_norm, nz * inv_norm)
end

function transducer_geometry(sim_params)
  @unpack transducer_pitch, transducer_pitch_elevation, aperture_size, aperture_size_elevation, aperture_radius, aperture_radius_elevation = sim_params

  n_transducers_azimuth = max(1, round(Int, aperture_size / transducer_pitch))
  n_transducers_elevation = max(1, round(Int, aperture_size_elevation / transducer_pitch_elevation))

  x_transducers = centered_axis_coordinates(n_transducers_azimuth, transducer_pitch)
  y_transducers = centered_axis_coordinates(n_transducers_elevation, transducer_pitch_elevation)

  if isfinite(aperture_radius)
    @assert maximum(abs.(x_transducers)) <= aperture_radius + 0.00001f0 "Aperture radius too small for given azimuth aperture size"
  end
  if isfinite(aperture_radius_elevation)
    @assert maximum(abs.(y_transducers)) <= aperture_radius_elevation + 0.00001f0 "Aperture radius too small for given elevation aperture size"
  end

  x_matrix = repeat(reshape(x_transducers, :, 1), 1, n_transducers_elevation)
  y_matrix = repeat(reshape(y_transducers, 1, :), n_transducers_azimuth, 1)
  z_matrix = curvature_sagitta.(x_matrix, aperture_radius) .+ curvature_sagitta.(y_matrix, aperture_radius_elevation)
  z_matrix .-= minimum(z_matrix)

  transducers = [SVector{3,Float32}(x_matrix[i, j], y_matrix[i, j], z_matrix[i, j]) for j in 1:n_transducers_elevation, i in 1:n_transducers_azimuth]
  transducer_normals = [element_normal(x_matrix[i, j], y_matrix[i, j], aperture_radius, aperture_radius_elevation) for j in 1:n_transducers_elevation, i in 1:n_transducers_azimuth]
  return transducers, transducer_normals
end

function beamplot_element_polygons(sim_params)
  @unpack beamplot_axes, transducer_pitch, transducer_pitch_elevation = sim_params

  transducers, transducer_normals = transducer_geometry(sim_params)
  lateral_index = beamplot_axes == :elevation_depth ? 2 : 1
  offaxis_index = beamplot_axes == :elevation_depth ? 1 : 2
  selected_pitch = slice_directivity_pitch(beamplot_axes, transducer_pitch, transducer_pitch_elevation)
  min_offaxis = minimum(abs(transducer[offaxis_index]) for transducer in transducers)
  offaxis_tolerance = max(selected_pitch * 1.0f-3, 1.0f-9)

  polygons = Vector{NTuple{4, SVector{2, Float32}}}()
  half_pitch = selected_pitch / 2.0f0
  for (transducer, normal) in zip(transducers, transducer_normals)
    if abs(abs(transducer[offaxis_index]) - min_offaxis) > offaxis_tolerance
      continue
    end

    center = SVector{2,Float32}(transducer[3], transducer[lateral_index])
    normal_in_slice = SVector{2,Float32}(normal[3], normal[lateral_index])
    inv_norm = inv(max(sqrt(sum(abs2, normal_in_slice)), 1.0f-6))
    normal_in_slice *= inv_norm
    tangent_in_slice = SVector{2,Float32}(-normal_in_slice[2], normal_in_slice[1])
    display_center = center + half_pitch * normal_in_slice

    push!(polygons, (
      display_center - half_pitch * normal_in_slice - half_pitch * tangent_in_slice,
      display_center - half_pitch * normal_in_slice + half_pitch * tangent_in_slice,
      display_center + half_pitch * normal_in_slice + half_pitch * tangent_in_slice,
      display_center + half_pitch * normal_in_slice - half_pitch * tangent_in_slice,
    ))
  end

  return polygons
end

# compute dependent parameters given global configuration
function init(trans_delays, sim_params)
  @unpack tx_frequency, transducer_pitch, transducer_pitch_elevation, aperture_size, aperture_size_elevation, aperture_radius, aperture_radius_elevation, spatial_res, temporal_res, end_simulation_time, fov, pulse_cycles, apodization_shape, beamplot_axes, focus_depth, focus_depth_elevation, steer_angle, steer_angle_elevation = sim_params

  n_transducers_azimuth = max(1, round(Int, aperture_size / transducer_pitch))
  n_transducers_elevation = max(1, round(Int, aperture_size_elevation / transducer_pitch_elevation))

  transducers, transducer_normals = transducer_geometry(sim_params)
  time_vec = collect(0.0f0:temporal_res:end_simulation_time)
  pulse_length = pulse_cycles / tx_frequency
  if apodization_shape == Rect
      apodization_matrix = ones(Float32, n_transducers_azimuth, n_transducers_elevation)
  else
      apodization_azimuth = n_transducers_azimuth == 1 ? ones(Float32, 1) : Float32[(sin((i - 1) * pi / (n_transducers_azimuth - 1)))^2 for i in 1:n_transducers_azimuth]
      apodization_elevation = n_transducers_elevation == 1 ? ones(Float32, 1) : Float32[(sin((i - 1) * pi / (n_transducers_elevation - 1)))^2 for i in 1:n_transducers_elevation]
      apodization_matrix = apodization_azimuth .* reshape(apodization_elevation, 1, :)
  end
  apodization_vec = vec(apodization_matrix)
  # Precompute pixel coordinates for the image grid to avoid recomputing every time step
  dx = Float32(fov[1]) / Float32(spatial_res[1])
  dz = Float32(fov[2]) / Float32(spatial_res[2])
  pix_coord_xs = Float32[((i - (spatial_res[1] + 1) / 2) * dx) for i in 1:spatial_res[1]]
  pix_coord_ys = Float32[j * dz for j in 1:spatial_res[2]]

  # Indices of transducers that will ever fire (non-negative delays)
  transducers_that_are_firing = findall(trans_delays .>= 0)

  return transducers, transducer_normals, time_vec, pulse_length, apodization_vec, pix_coord_xs, pix_coord_ys, transducers_that_are_firing
end

@inline function axis_delay_law(axis_coords::AbstractVector{Float32}, focus_depth::Float32, steer_angle::Float32, c::Float32)
  # Compute a 1D delay law for one aperture axis, including steering.
  delays = similar(axis_coords)
  if isinf(focus_depth)
    for i in eachindex(axis_coords)
      delays[i] = axis_coords[i] * sind(steer_angle) / c
    end
  elseif focus_depth > 0.0f0
    focus_distance = focus_depth
    focus_coord = focus_depth * tand(steer_angle)
    for i in eachindex(axis_coords)
      coord = axis_coords[i]
      delays[i] = -sqrt((coord - focus_coord)^2 + focus_distance^2) / c
    end
  elseif focus_depth < 0.0f0
    # Negative focus depth reuses the old virtual-source construction for diverging waves.
    virtual_source_coord = -sind(steer_angle) * abs(focus_depth)
    virtual_source_depth = -cosd(steer_angle) * abs(focus_depth)
    for i in eachindex(axis_coords)
      coord = axis_coords[i]
      delays[i] = sqrt((coord - virtual_source_coord)^2 + virtual_source_depth^2) / c
    end
  else
    error("Invalid focus depth value 0")
  end
  delays .-= minimum(delays)
  return delays
end

# Compute transmit time delays for transducer elements given aperture_size [m] (impacts how many elements are firing) and focus depth [m] and steer angle [deg] (both impact amount of delay).
# Note: Steer angle is about center of active aperture.
# Since current implementation places the active aperture about the center of the physical aperture, steer angle is also about center of physical aperture.
# Focus equation adapted from: Tumsys, 2014, http://dx.doi.org/10.5755/j01.eee.20.3.3638
# Steer equation adapted from: Ramm, 1983, http://dx.doi.org/10.1109/TBME.1983.325149
function delays_from_focus_and_steer(sim_params)
  @unpack focus_depth, focus_depth_elevation, steer_angle, steer_angle_elevation, aperture_size, aperture_size_elevation, c, transducer_pitch, transducer_pitch_elevation = sim_params

  n_transducers_azimuth = max(1, round(Int, aperture_size / transducer_pitch))
  n_transducers_elevation = max(1, round(Int, aperture_size_elevation / transducer_pitch_elevation))

  x_transducers = centered_axis_coordinates(n_transducers_azimuth, transducer_pitch)
  y_transducers = centered_axis_coordinates(n_transducers_elevation, transducer_pitch_elevation)

  azimuth_delays = axis_delay_law(x_transducers, focus_depth, steer_angle, c)
  elevation_delays = axis_delay_law(y_transducers, focus_depth_elevation, steer_angle_elevation, c)

  trans_delays = Float32[azimuth_delays[i] + elevation_delays[j] for j in 1:n_transducers_elevation, i in 1:n_transducers_azimuth]
  trans_delays = vec(trans_delays)
  trans_delays .-= minimum(trans_delays)

  return trans_delays
end

# Compute temporal and spatial resolutions fine enough to support the pulse,
# and end of simulation time just long enough to reach the corner of the FOV.
function autores(sim_params, trans_delays; multiplier=1.0f0, min_spatial_res=(256, 512))
  @unpack tx_frequency, fov, aperture_size, aperture_size_elevation, c, pulse_cycles = sim_params
  # For our implementation, temporal resolution is way more important than spatial resolution
  # we're not interested in the detailed look of the pulse cycle, but rather
  # at each pixel, we need good temporal sampling for correct interference buildup.
  temporal_res = 1/tx_frequency / 8 / multiplier  # 8 time points per cycle
  wavelength = c / tx_frequency
  spatial_res_v = Int.(round.(fov / wavelength)) * 4 * multiplier # 4 samples per wavelength
  spatial_res_v = max.(spatial_res_v, min_spatial_res)  # but at least min_spatial_res, for human visualization purposes.
  spatial_res = SVector(spatial_res_v[1], spatial_res_v[2])
  # End simulation when pulse reaches outside of FOV (for "worst" case i.e. longest path).
  # since setup is symmetric, computing one side is enough.
  max_aperture_size = max(aperture_size, aperture_size_elevation)
  time_top_trans_to_bot_corner = sqrt((fov[1]/2 + max_aperture_size/2)^2 + fov[2]^2) / c
  end_simulation_time = time_top_trans_to_bot_corner + maximum(trans_delays) + pulse_cycles * 1/tx_frequency
  WaveSimParameters(sim_params; temporal_res=temporal_res, spatial_res=spatial_res, end_simulation_time=end_simulation_time)
end

@inline function slice_pixel_coordinates(beamplot_axes::Symbol, lateral::Float32, depth::Float32)
  # Convert the 2D image grid into the requested physical slice coordinates.
  if beamplot_axes == :elevation_depth
    return 0.0f0, lateral, depth
  end
  return lateral, 0.0f0, depth
end

@inline function slice_directivity_pitch(beamplot_axes::Symbol, transducer_pitch::Float32, transducer_pitch_elevation::Float32)
  # Use the element width that corresponds to the simulated beamplot plane.
  return beamplot_axes == :elevation_depth ? transducer_pitch_elevation : transducer_pitch
end

function cuda_backend_available()
  ext = Base.get_extension(@__MODULE__, :WaveSimCUDAExt)
  return ext !== nothing && ext.cuda_backend_available_impl()
end

function wavesim_cuda(trans_delays, sim_params)
  ext = Base.get_extension(@__MODULE__, :WaveSimCUDAExt)
  if ext === nothing
    error("CUDA backend is unavailable. Install CUDA.jl, import CUDA, and load WaveSim in an environment with a functional NVIDIA GPU.")
  end
  return ext.wavesim_cuda_impl(trans_delays, sim_params)
end

# simulate one time step of wave propagation
function simulate_one_time_step!(image, t, trans_delays, pulse_length, tx_frequency, c, spatial_res, pulse_shape_func, apodization_vec, directivity_func, transducer_pitch, attenuation_coefficient, pix_coord_xs, pix_coord_ys, transducers, transducer_normals, beamplot_axes, transducers_that_are_firing)
  one_over_c = 1.0f0 / c  # For computation speed improvement [s/m]

  # Attenuation [dB] = α [dB/cm/MHz] * distance [cm] * frequency [MHz]
  # For linear amplitude reduction (A = A0 * 10^(-dB/20)):
  # amplitude_factor = 10 ^ (- (α * dist_cm * freq_MHz) / 20)
  α_factor = - (attenuation_coefficient * tx_frequency / 1f6) / 20.0f0

  @inbounds for y in 1:spatial_res[2]
    pix_coord_y = pix_coord_ys[y]
    @inbounds for x in 1:spatial_res[1]
      pix_coord_x = pix_coord_xs[x]
      pix_x, pix_y, pix_z = slice_pixel_coordinates(beamplot_axes, pix_coord_x, pix_coord_y)
      amp = 0.0f0
      @simd for i_trans in transducers_that_are_firing
        t_coord = transducers[i_trans]
        n_coord = transducer_normals[i_trans]
        trans_delay = trans_delays[i_trans]
        dx = pix_x - t_coord[1]
        dy = pix_y - t_coord[2]
        dz = pix_z - t_coord[3]
        dist_to_transducer_not_zero = max(sqrt(dx*dx + dy*dy + dz*dz), 1.0f-6)
        wave_spreading_factor = 1.0f0 / dist_to_transducer_not_zero
        time_to_transducer = dist_to_transducer_not_zero * one_over_c
        time_to_reach = time_to_transducer + trans_delay
        if time_to_reach <= t <= time_to_reach + pulse_length
          # directivity
          cos_theta = clamp((dx * n_coord[1] + dy * n_coord[2] + dz * n_coord[3]) / dist_to_transducer_not_zero, -1.0f0, 1.0f0)
          directivity_factor = 0.0f0
          if cos_theta >= 0.0f0  # are we in front of the transducer? if not, directivity is 0 and we can skip the expensive acos and directivity function evaluation
            θ = acos(cos_theta)
            directivity_factor = directivity_func(θ, tx_frequency, c, transducer_pitch)
          end

          dist_cm = dist_to_transducer_not_zero * 100.0f0
          attenuation_factor = 10.0f0 ^ (α_factor * dist_cm)

          # The phase is in radians (the π factor in the more common notation "2*π*f" is inside pulse_shape_func's cospi)
          amp += apodization_vec[i_trans] * pulse_shape_func((t - time_to_reach) * tx_frequency * 2.0f0) * wave_spreading_factor * directivity_factor * attenuation_factor
        end
      end
      image[x,y] = amp
    end
  end
end

# run the simulation time steps
function _wavesim_cpu(trans_delays, sim_params)
  @unpack tx_frequency, pulse_cycles, spatial_res, c, fov, pulse_shape_func, directivity_func, transducer_pitch, transducer_pitch_elevation, attenuation_coefficient, beamplot_axes = sim_params
  transducers, transducer_normals, tvec, pulse_length, apodization_vec, pix_coord_xs, pix_coord_ys, transducers_that_are_firing = init(trans_delays, sim_params)
  images = zeros(Float32, (spatial_res[1], spatial_res[2], length(tvec)))
  effective_transducer_pitch = slice_directivity_pitch(beamplot_axes, transducer_pitch, transducer_pitch_elevation)

  @showprogress Threads.@threads for i_time in 1:length(tvec)
    t = tvec[i_time]
    image = view(images, :, :, i_time)
    simulate_one_time_step!(image, t, trans_delays, pulse_length, tx_frequency, c, spatial_res, pulse_shape_func, apodization_vec, directivity_func, effective_transducer_pitch, attenuation_coefficient, pix_coord_xs, pix_coord_ys, transducers, transducer_normals, beamplot_axes, transducers_that_are_firing)
  end

  return images
end

function wavesim(trans_delays, sim_params; backend::Symbol = :cpu)
  if backend === :cpu
    return _wavesim_cpu(trans_delays, sim_params)
  elseif backend === :cuda
    return wavesim_cuda(trans_delays, sim_params)
  else
    error("Unknown backend: $backend. Expected :cpu or :cuda.")
  end
end

# Get the beam profile spatial map and transmit time of beam energy, where each pixel indicates the maximum energy that was seen at that place, and at what time.
function beam_energy_map_and_transmit_time_map(images, sim_params)
    @unpack temporal_res = sim_params

    integrated_energy_map = dropdims(sum(images.^2, dims=3), dims=3)  # SOS integration
    # integrated_energy_map = dropdims(sqrt.(mean(images.^2, dims=3)), dims=3)  # RMS
    # integrated_energy_map = dropdims(mean(abs.(images), dims=3), dims=3)  # MeanAbs

    # peak to peak squared amplitude
    maxval, maxlinindices = findmax(images, dims=3)
    minval, minlinindices = findmin(images, dims=3)
    peak_to_peak_map = dropdims(maxval .- minval, dims=3).^2

    transmit_time_map = similar(peak_to_peak_map)
    for linind in eachindex(maxlinindices)
        x, y, t = Tuple(CartesianIndices(images)[maxlinindices[linind]])
        transmit_time_map[linind] = t * temporal_res
    end

    peak_to_peak_time_delta_map = similar(peak_to_peak_map)
    for linind in eachindex(maxlinindices)
        x, y, tmax = Tuple(CartesianIndices(images)[maxlinindices[linind]])
        x, y, tmin = Tuple(CartesianIndices(images)[minlinindices[linind]])
        peak_to_peak_time_delta_map[linind] = (tmax - tmin) * temporal_res
    end

    # maximum windowed energy map. window is size of pulse length, and we take the maximum over time of the windowed energy.
    window_size = round(Int, sim_params.pulse_cycles / sim_params.tx_frequency / sim_params.temporal_res)
    windowed_energy_map = similar(peak_to_peak_map)
    for linind in eachindex(maxlinindices)
        x, y, t = Tuple(CartesianIndices(images)[linind])
        # energy_over_time = [sqrt.(mean(images[x, y, t:t+window_size-1].^2)) for t in 1:(size(images, 3)-window_size+1)]
        energy_over_time = [sum(images[x, y, t:t+window_size-1].^2) for t in 1:(size(images, 3)-window_size+1)]
        windowed_energy_map[x, y] = maximum(energy_over_time)
    end

    return windowed_energy_map, integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map
end

end  # module WaveSim
