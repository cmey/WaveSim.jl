# WaveSim.jl, wave propagation simulator.
# Christophe Meyer, 2016-2024
module WaveSim

using Parameters
using ProgressMeter
using SpecialFunctions
using StaticArrays
using Statistics

export WaveSimParameters

@enum ApodizationShape Rect Hann

# Configuration
@with_kw struct WaveSimParameters{pulse_shape_func_T}
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
  # Distance from transducer to focus point, Inf means plane wave.
  focus_depth::Float32 = 0.03  # [m]
  # Steering angle in azimuth.
  steer_angle::Float32 = 10.0  # [deg]
  # Shape of the transmit pulse.
  pulse_shape_func::pulse_shape_func_T = phase_over_pi -> cospi(phase_over_pi)
  # Shape of the transmit aperture apodization.
  apodization_shape::ApodizationShape = Rect
  # Display dynamic range.
  dbrange::Float32 = 40  # [dB]
  # Plot orientation: beam is horizontal or vertical.
  orientation::Symbol = :horizontal
  # Directivity function of angle and frequency.
  directivity_func::Function = function(θ, tx_frequency, c, transducer_pitch)
    # account for mechanical crosstalk between elements
    crosstalk_factor = 1.2f0  # assume 20% crosstalk factor
    element_surface_diameter = transducer_pitch * crosstalk_factor
    # from Umchid 2009 Directivity Pattern Measurement of Ultrasound Transducers
    # https://www.thaiscience.info/Journals/Article/IABE/10892457.pdf
    a = element_surface_diameter / 2
    λ = c / tx_frequency
    k = 2.0f0 * π / λ
    ka_sin_θ = k * a * sin(θ) + 0.00001f0  # avoid division by zero, assuming θ ∈ [-π/2 to π/2]
    h = convert(Float32, 2 * besselj1(ka_sin_θ) / ka_sin_θ)
    return h
  end
end


# compute temporal and spatial resolutions fine enough to support the pulse,
# and end of simulation time just long enough to reach the corner of the FOV.
function autores(sim_params, trans_delays; multiplier=1.0f0)
  @unpack tx_frequency, fov, aperture_size, c, pulse_cycles = sim_params
  # for our implementation, temporal resolution is way more important than spatial resolution
  # we're not interested in the detailed look of the pulse cycle, but rather
  # at each pixel, we need good temporal sampling for correct interference buildup.
  temporal_res = 1/tx_frequency / 8 / multiplier  # 8 time points per cycle
  wavelength = c / tx_frequency
  spatial_res_v = Int.(round.(fov / wavelength)) * 4 * multiplier # 4 samples per wavelength
  spatial_res_v = max(spatial_res_v, [256, 512])  # but at least 256x512, for human visualization purposes.
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
  if aperture_radius != Inf
    x_transducers_centered = x_transducers .- mean(extrema(x_transducers))
    @assert aperture_radius >= sum(extrema(x_transducers_centered))/2
    y_transducers .= sqrt.( aperture_radius.^2 .- x_transducers_centered.^2)
    y_transducers .-= minimum(y_transducers)
  end
  transducers = collect(zip(x_transducers, y_transducers))
  time_vec = collect(0.0f0:temporal_res:end_simulation_time)
  pulse_length = pulse_cycles / tx_frequency
  if apodization_shape == Rect
      apodization_vec = ones(Float32, n_transducers)
  else
      apodization_vec = Float32[(sin(i*pi/(n_transducers-1)))^2 for i in 1:n_transducers]
  end
  return transducers, time_vec, pulse_length, apodization_vec
end


# simulate one time step of wave propagation
function simulate_one_time_step!(image, t, trans_delays, transducers, pulse_length, tx_frequency, c, spatial_res, fov, pulse_shape_func, apodization_vec, directivity_func, transducer_pitch)
  # println("simulate_one_time_step")
  image_pitch = fov ./ spatial_res
  transducers_that_are_firing = findall(trans_delays .>= 0)
  trans_coord_x::Float32 = 0  # [m] space
  trans_coord_y::Float32 = 0  # [m] space
  pix_coord_x::Float32 = 0  # [m] space
  pix_coord_y::Float32 = 0  # [m] space
  one_over_c = 1.0f0 / c  # For computation speed improvement [s/m]

  for y in 1:spatial_res[2]
    pix_coord_y = y * image_pitch[2]  # [m]
    for x in 1:spatial_res[1]
      pix_coord_x = x * image_pitch[1]  # [m]

      # for transducers that will fire...
      amp = 0.0f0
      @inbounds for i_trans in transducers_that_are_firing
        trans_coord_x, trans_coord_y = transducers[i_trans]  # index
        trans_delay = trans_delays[i_trans]
        # from pixel to transducer
        dist_to_transducer = sqrt((trans_coord_x - pix_coord_x)^2 +
                                  (trans_coord_y - pix_coord_y)^2)
        # wave is spreading spherically in space, energy spreading loss goes with r^2, amp with r
        # see: https://ccrma.stanford.edu/~jos/pasp/Spherical_Waves_Point_Source.html
        wave_spreading_factor = 1.0f0 / dist_to_transducer  # no perf diff if put inside conditional
        time_to_transducer = dist_to_transducer * one_over_c
        time_to_reach = time_to_transducer + trans_delay
        # if the transducer wave has reached this pixel...
        if time_to_reach <= t <= time_to_reach + pulse_length
          θ = atan(pix_coord_x - trans_coord_x, pix_coord_y - trans_coord_y)
          directivity_factor = directivity_func(θ, tx_frequency, c, transducer_pitch)
          # wave interference
          amp += apodization_vec[i_trans] * pulse_shape_func((t - time_to_reach) * tx_frequency * 2.0f0) * wave_spreading_factor * directivity_factor
        end
      end
      image[x,y] = amp  # write back to memory only once per pixel (precompute stack variable amp for all transducers)
    end
  end
end


# run the simulation time steps
function wavesim(trans_delays, sim_params)
  @unpack tx_frequency, pulse_cycles, spatial_res, c, fov, pulse_shape_func, directivity_func, transducer_pitch = sim_params
  transducers, tvec, pulse_length, apodization_vec = init(trans_delays, sim_params)

  images = zeros(Float32, (spatial_res[1], spatial_res[2], length(tvec)))

  @showprogress Threads.@threads for i_time in 1:length(tvec)
    t = tvec[i_time]
    image = view(images, :, :, i_time)
    simulate_one_time_step!(image, t, trans_delays, transducers, pulse_length, tx_frequency, c, spatial_res, fov, pulse_shape_func, apodization_vec, directivity_func, transducer_pitch)
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
    x_transducers .+= fov[1] / 2 - mean(extrema(x_transducers))  # center around FOV in x

    trans_delays = zeros(Float32, num_elements)

    # Focus
    focus_coord = [0.0, focus_depth]
    if Inf == focus_depth
        # Focus at infinity is considered planewave. Zero delay is fine.
    else
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
        for (i_elem, x_elem) = enumerate(x_transducers)
            k = i_elem
            second_part = sqrt(l02 + (Ae * abs(n - 2k + 1) / (2n - 1)) ^ 2)  # offset_dist
            trans_delays[i_elem] = (first_part - second_part) / v1
        end
    end

    # Steer
    # Define variables as in the reference paper (see above).
    θ = deg2rad(steer_angle)
    d = transducer_pitch
    # n * d = x_elem
    # for transducers that will fire...
    for (i_elem, x_elem) = enumerate(x_transducers)
        trans_delays[i_elem] += x_elem / c * sin(θ)
    end

    # Causal delays
    trans_delays .-= minimum(trans_delays)

    return trans_delays
end


end  # module WaveSim
