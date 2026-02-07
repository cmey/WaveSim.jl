# WaveSim.jl, wave propagation simulator.
# Christophe Meyer, 2016-2026
module WaveSim

using Parameters
using ProgressMeter
using SpecialFunctions
using StaticArrays
using Statistics

export WaveSimParameters

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
  end_simulation_time::Float32 = 56.0e-6  # [s] starts at 0 s
  # Resolution at which to divide the simulation time span.
  temporal_res::Float32 = 0.1e-6  # [s]
  # Field of view, x then z. x=0 centered on aperture center, z=0 at aperture plane.
  fov::SVector{2, Float32} = @SVector [4e-2, 8e-2]  # [m]
  # Spatial resolution of the simulation, x then z.
  spatial_res::SVector{2, Int} = @SVector [256, 512]  # [pixels]
  # Spacing between physical elements of transducer array.
  transducer_pitch::Float32 = 208e-6  # [m]
  # Active aperture size of transducer (centered)
  aperture_size::Float32 = 0.02  # [m]
  # Radius of curvature of the aperture. Inf means a flat aperture.
  aperture_radius = Inf  # [m]
  # Distance from transducer to focus point, Inf and -Inf mean plane wave.
  # Positive values mean focusing, negative values mean diverging.
  focus_depth::Float32 = 0.03  # [m]
  # Steering angle in azimuth.
  steer_angle::Float32 = 10.0  # [deg]
  # Shape of the transmit pulse.
  pulse_shape_func::pulse_shape_func_T = phase_divided_by_pi -> cospi(phase_divided_by_pi)
  # Shape of the transmit aperture apodization.
  apodization_shape::ApodizationShape = Rect
  # Display dynamic range.
  dbrange::Float32 = 40  # [dB]
  # Plot orientation: beam is horizontal or vertical.
  orientation::Symbol = :horizontal
  # Directivity function of angle and frequency.
  directivity_func::directivity_func_T = default_directivity
end

# Concrete default directivity function with fully-typed signature to help inference
function default_directivity(sinθ::Float32, tx_frequency::Float32, c::Float32, transducer_pitch::Float32)::Float32
  # Account for mechanical crosstalk between elements
  crosstalk_factor = 1.2f0  # Assume 20% crosstalk factor
  # From Umchid 2009 Directivity Pattern Measurement of Ultrasound Transducers:
  # https://www.thaiscience.info/Journals/Article/IABE/10892457.pdf
  element_surface_diameter = transducer_pitch * crosstalk_factor
  a = element_surface_diameter / 2.0f0
  λ = c / tx_frequency
  k = 2.0f0 * π / λ
  ka_sin_θ = k * a * sinθ
  if abs(ka_sin_θ) < 1f-5
    return 1.0f0
  end
  return convert(Float32, 2.0f0 * besselj1_approx(ka_sin_θ) / ka_sin_θ)
  # return convert(Float32, 2.0f0 * besselj1(ka_sin_θ) / ka_sin_θ)
end

# CUDA-friendly polynomial approximation for besselj1(x) for Float32 inputs
@inline function besselj1_approx(x::Float32)
  xk = x
  x2 = xk * xk
  return xk * (0.5f0 - x2 * (1/16f0 - x2 * (1/384f0)))
end

# Compute temporal and spatial resolutions fine enough to support the pulse,
# and end of simulation time just long enough to reach the corner of the FOV.
function autores(sim_params, trans_delays; multiplier=1.0f0, min_spatial_res=(256, 512))
  @unpack tx_frequency, fov, aperture_size, c, pulse_cycles = sim_params
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
  time_top_trans_to_bot_corner = sqrt((fov[1]/2 + aperture_size/2)^2 + fov[2]^2) / c
  end_simulation_time = time_top_trans_to_bot_corner + maximum(trans_delays) + pulse_cycles * 1/tx_frequency
  WaveSimParameters(sim_params; temporal_res=temporal_res, spatial_res=spatial_res, end_simulation_time=end_simulation_time)
end


# compute dependent parameters given global configuration
function init(trans_delays, sim_params)
  @unpack tx_frequency, transducer_pitch, aperture_radius, spatial_res, temporal_res, end_simulation_time, fov, pulse_cycles, apodization_shape = sim_params
  n_transducers = length(trans_delays)
  x_transducers = [transducer_pitch * itrans for itrans in 1:n_transducers]  # [m]
  x_transducers .+= fov[1] / 2 - mean(extrema(x_transducers))  # center around FOV in x
  y_transducers = zeros(Float32, n_transducers)
  if isfinite(aperture_radius)
    x_transducers_centered = x_transducers .- mean(extrema(x_transducers))
    x_half_width = maximum(abs.(x_transducers_centered))
    @assert aperture_radius >= x_half_width + 0.00001f0 "Aperture radius too small for given aperture size"
    y_transducers .= sqrt.( aperture_radius.^2 .- x_transducers_centered.^2)
    y_transducers .-= minimum(y_transducers)
  end
  transducers = [SVector{2,Float32}(x_transducers[i], y_transducers[i]) for i in 1:length(x_transducers)]
  time_vec = collect(0.0f0:temporal_res:end_simulation_time)
  pulse_length = pulse_cycles / tx_frequency
  if apodization_shape == Rect
      apodization_vec = ones(Float32, n_transducers)
  else
      apodization_vec = Float32[(sin(i*pi/(n_transducers-1)))^2 for i in 1:n_transducers]
  end
  # Precompute pixel coordinates for the image grid to avoid recomputing every time step
  dx = Float32(fov[1]) / Float32(spatial_res[1])
  dz = Float32(fov[2]) / Float32(spatial_res[2])
  pix_coord_xs = Float32[i * dx for i in 1:spatial_res[1]]
  pix_coord_ys = Float32[j * dz for j in 1:spatial_res[2]]

  # Indices of transducers that will ever fire (non-negative delays)
  transducers_that_are_firing = findall(trans_delays .>= 0)

  return transducers, time_vec, pulse_length, apodization_vec, pix_coord_xs, pix_coord_ys, transducers_that_are_firing
end


# simulate one time step of wave propagation
function simulate_one_time_step!(image, t, trans_delays, pulse_length, tx_frequency, c, spatial_res, pulse_shape_func, apodization_vec, directivity_func, transducer_pitch, pix_coord_xs, pix_coord_ys, transducers, transducers_that_are_firing)
  one_over_c = 1.0f0 / c  # For computation speed improvement [s/m]

  @inbounds for y in 1:spatial_res[2]
    pix_coord_y = pix_coord_ys[y]
    @inbounds for x in 1:spatial_res[1]
      pix_coord_x = pix_coord_xs[x]
      amp = 0.0f0
      @simd for i_trans in transducers_that_are_firing
        t_coord = transducers[i_trans]
        trans_delay = trans_delays[i_trans]
        dx = pix_coord_x - t_coord[1]
        dy = pix_coord_y - t_coord[2]
        dist_to_transducer_not_zero = max(sqrt(dx*dx + dy*dy), 1.0f-6)
        wave_spreading_factor = 1.0f0 / dist_to_transducer_not_zero
        time_to_transducer = dist_to_transducer_not_zero * one_over_c
        time_to_reach = time_to_transducer + trans_delay
        if time_to_reach <= t <= time_to_reach + pulse_length
          # compute sin(theta) = opposite / hypotenuse where opposite = horizontal difference
          sinθ = dx / dist_to_transducer_not_zero
          directivity_factor = directivity_func(sinθ, tx_frequency, c, transducer_pitch)
          # The phase is in radians (the π factor in the more common notation "2*π*f" is inside pulse_shape_func's cospi)
          amp += apodization_vec[i_trans] * pulse_shape_func((t - time_to_reach) * tx_frequency * 2.0f0) * wave_spreading_factor * directivity_factor
        end
      end
      image[x,y] = amp
    end
  end
end


# run the simulation time steps
function wavesim(trans_delays, sim_params)
  @unpack tx_frequency, pulse_cycles, spatial_res, c, fov, pulse_shape_func, directivity_func, transducer_pitch = sim_params
  transducers, tvec, pulse_length, apodization_vec, pix_coord_xs, pix_coord_ys, transducers_that_are_firing = init(trans_delays, sim_params)
  images = zeros(Float32, (spatial_res[1], spatial_res[2], length(tvec)))

  @showprogress Threads.@threads for i_time in 1:length(tvec)
    t = tvec[i_time]
    image = view(images, :, :, i_time)
    simulate_one_time_step!(image, t, trans_delays, pulse_length, tx_frequency, c, spatial_res, pulse_shape_func, apodization_vec, directivity_func, transducer_pitch, pix_coord_xs, pix_coord_ys, transducers, transducers_that_are_firing)
  end

  return images
end


# Get the beam profile spatial map and transmit time of beam energy, where each pixel indicates the maximum energy that was seen at that place, and at what time.
function beam_energy_map_and_transmit_time_map(images, sim_params)
    @unpack temporal_res = sim_params

    # integrated_energy_map = dropdims(sqrt.(mean(images.^2, dims=3)), dims=3)  # RMS
    integrated_energy_map = dropdims(mean(abs.(images), dims=3), dims=3)  # MeanAbs

    # peak to peak amplitude
    maxval, maxlinindices = findmax(images, dims=3)
    minval, minlinindices = findmin(images, dims=3)
    peak_to_peak_map = dropdims(maxval .- minval, dims=3)

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

    return integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map
end


# Compute transmit time delays for transducer elements given aperture_size [m] (impacts how many elements are firing) and focus depth [m] and steer angle [deg] (both impact amount of delay).
# Note: Steer angle is about center of active aperture.
# Since current implementation places the active aperture about the center of the physical aperture, steer angle is also about center of physical aperture.
# Focus equation adapted from: Tumsys, 2014, http://dx.doi.org/10.5755/j01.eee.20.3.3638
# Steer equation adapted from: Ramm, 1983, http://dx.doi.org/10.1109/TBME.1983.325149
function delays_from_focus_and_steer(sim_params)
    @unpack focus_depth, steer_angle, aperture_size, c, transducer_pitch, fov = sim_params

    num_elements = round(Int, aperture_size / transducer_pitch)

    x_transducers = [transducer_pitch * ielem for ielem in 1:num_elements]  # [m]
    x_transducers .+= -mean(extrema(x_transducers))  # center around 0 in x

    # Focus
    if focus_depth == Inf || focus_depth == -Inf
        # Focus at infinity is considered planewave.
        first_firing_element = -steer_angle < 0 ? x_transducers[1] : x_transducers[end]
        sin_a = abs(sind(steer_angle))  # sind : sine of degrees
        trans_delays = Float32[abs(element_location - first_firing_element) * sin_a / c for element_location in x_transducers]
    elseif focus_depth > 0
        # Focused beam.
        # Define variables as in the reference paper (see above).
        l0 = focus_depth
        A = aperture_size
        e = transducer_pitch
        n = num_elements
        v1 = c
        # Compute dependent variables.
        l02 = l0^2
        Ae = A - e
        first_part = sqrt(l02 + (Ae / 2) ^ 2)  # longest dist from focus to elem (i.e. to the edge)

        # for transducers that will fire...
        trans_delays = zeros(Float32, num_elements)
        for (i_elem, x_elem) = enumerate(x_transducers)
            k = i_elem
            second_part = sqrt(l02 + (Ae * abs(n - 2k + 1) / (2n - 1)) ^ 2)  # offset_dist
            trans_delays[i_elem] = (first_part - second_part) / v1
        end

        # Separately add steering
        # Define variables as in the reference paper (see above).
        θ = deg2rad(steer_angle)
        d = transducer_pitch
        # n * d = x_elem
        # for transducers that will fire...
        for (i_elem, x_elem) = enumerate(x_transducers)
            trans_delays[i_elem] += x_elem / c * sin(θ)
        end
    elseif focus_depth < 0
        # Diverging beam.
        virtual_source_x = -sind(steer_angle) * abs(focus_depth)
        virtual_source_y = -cosd(steer_angle) * abs(focus_depth)
        trans_delays = Float32[sqrt((element_location - virtual_source_x)^2 + virtual_source_y^2) / c for element_location in x_transducers]
    else
        error("Invalid focus depth value 0")
    end

    # Causal delays
    trans_delays .-= minimum(trans_delays)

    return trans_delays
end


end  # module WaveSim
