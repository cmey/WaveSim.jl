# WaveSim.jl, wave propagation simulator.
# Wave propagation computation.
# Christophe MEYER, 2016-2018
module WaveSim

using Parameters

export WaveSimParameters

# Configuration
@with_kw struct WaveSimParameters
  # Transmit center frequency.
  tx_frequency::Float64 = 3_000_000.0  # [Hz]
  # Length of the transmit pulse.
  pulse_cycles::Float64 = 2  # [cycles]
  # Speed of sound.
  c::Float64 = 1540.0  # [m/s]
  # Resolution at which to divide the simulation time span.
  temporal_res::Float64 = 0.1e-6  # [s]
  # Time span to simulate.
  end_simulation_time::Float64 = 56.0e-6  # [s] starts at 0 s
  # Spatial resolution of the simulation.
  spatial_res::Array{Int, 1} = [256, 512]  # [pixels]
  # Field of view, x then z. x=0 centered on aperture center, z=0 at aperture plane.
  fov::Array{Float64, 1} = [4e-2, 8e-2]  # [m]
  # Physical length of the transducer array.
  transducer_array_size::Float64 = 0.03  # [m]
  # Spacing between physical elements of transducer array.
  transducer_pitch::Float64 = 208e-6  # [m]
  # Active aperture size of transducer (centered)
  aperture_size::Float64 = 0.02  # [m]
  # Distance from transducer to focus point, Inf means plane wave.
  focus_depth::Float64 = 0.03  # [m]
  # Steering angle in azimuth.
  steer_angle::Float64 = 10.0  # [deg]
  # Shape of the transmit pulse.
  pulse_shape_func = phase -> cos(phase)
  # Display range.
  dbrange::Float64 = 40
end


# compute dependent parameters given global configuration
function init(trans_delays, sim_params)
  @unpack tx_frequency, transducer_pitch, spatial_res, temporal_res, end_simulation_time, fov, pulse_cycles = sim_params
  n_transducers = length(trans_delays)
  x_transducers = [transducer_pitch * itrans for itrans in 1:n_transducers]  # [m]
  x_transducers += fov[1] / 2 - mean(extrema(x_transducers))  # center around FOV in x
  time_vec = collect(0.0:temporal_res:end_simulation_time)
  pulse_length = pulse_cycles / tx_frequency
  return x_transducers, time_vec, pulse_length
end


# simulate one time step of wave propagation
function simulate_one_time_step!(image, t, trans_delays, x_transducers, pulse_length, sim_params)
  # println("simulate_one_time_step")
  @unpack tx_frequency, c, spatial_res, fov, pulse_shape_func = sim_params
  image_pitch = fov ./ spatial_res
  transducers_that_are_firing = find(trans_delays .>= 0)
  trans_coord = [0.0, 0.0]  # [m] space
  pix_coord = [0.0, 0.0]
  one_over_c = 1.0 / c  # For computation speed improvement [s/m]

  for y in 1:spatial_res[2]
    pix_coord[2] = y * image_pitch[2]  # [m]
    for x in 1:spatial_res[1]
      pix_coord[1] = x * image_pitch[1]  # [m]

      # for transducers that will fire...
      @inbounds for i_trans in transducers_that_are_firing
        xt = x_transducers[i_trans]  # index
        trans_delay = trans_delays[i_trans]
        trans_coord[1] = xt
        trans_coord[2] = 0.0
        # from pixel to transducer
        dist_to_transducer = sqrt((trans_coord[1] - pix_coord[1])^2 +
                                  (trans_coord[2] - pix_coord[2])^2)
        # wave is spreading spherically in space, energy spreading loss goes with r^2, amp with r
        # see: https://ccrma.stanford.edu/~jos/pasp/Spherical_Waves_Point_Source.html
        wave_spreading_factor = 1 / dist_to_transducer
        time_to_transducer = dist_to_transducer * one_over_c
        time_to_reach = time_to_transducer + trans_delay
        # if the transducer wave reached this pixel...
        if time_to_reach <= t <= time_to_reach + pulse_length
          amp = pulse_shape_func((t - time_to_reach) * tx_frequency * 2pi) * wave_spreading_factor
        else
          amp = 0.0
        end
        # wave interference
        image[x,y] += amp
      end
    end
  end
end


# run the simulation time steps
function wavesim(trans_delays, sim_params)
  @unpack tx_frequency, pulse_cycles, spatial_res = sim_params
  x_transducers, tvec, pulse_length = init(trans_delays, sim_params)

  images = zeros(Float32, (spatial_res[1], spatial_res[2], length(tvec)))

  Threads.@threads for i_time in 1:length(tvec)
    t = tvec[i_time]
    image = view(images, :, :, i_time)
    simulate_one_time_step!(image, t, trans_delays, x_transducers, pulse_length, sim_params)
  end

  return images
end


# Get the beam profile spatial map and transmit time of beam energy, where each pixel indicates the maximum energy that was seen at that place, and at what time.
function beam_energy_map_and_transmit_time_map(images, sim_params)
    @unpack temporal_res = sim_params

    maxval, linindices = findmax(images, 3)

    beam_energy_map = squeeze(maxval, 3)

    transmit_time_map = similar(beam_energy_map)

    for linind in eachindex(linindices)
        x, y, t = ind2sub(images, linindices[linind])
        transmit_time_map[linind] = t * temporal_res
    end

    return beam_energy_map, transmit_time_map
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
    x_transducers += fov[1] / 2 - mean(extrema(x_transducers))  # center around FOV in x

    trans_delays = zeros(num_elements)

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
    trans_delays -= minimum(trans_delays)

    return trans_delays
end


end  # module WaveSim
