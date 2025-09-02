using WaveSim

# Run the simulator and display results.
function main()
  # Define simulation parameters (use many default values, see WaveSimParameters).
  sim_params = WaveSimParameters(
    focus_depth = 0.03,  # [m]
    steer_angle = 10.0,  # [deg]
    aperture_size = 0.02,  # [m]
    # tx_frequency = 3_000_000.0,  # [Hz]
    # focus_depth = Inf,  # [m]
    # steer_angle = 0.0,  # [deg]
    # transducer_pitch = 205e-6,  # [m]
    # aperture_size = 0.01312,  # [m]
    dbrange = 50,  # [dB]
  )

  # Compute focusing delays for the elements of the phased array.
  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)

  # Compute optimized "best" spatial and temporal parameters.
  sim_params = WaveSim.autores(sim_params, trans_delays)

  # Run the simulation.
  images = WaveSim.wavesim(trans_delays, sim_params)
  beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params)

  return beam_energy_map, transmit_time_map, images, sim_params
end

beam_energy_map, transmit_time_map, images, sim_params = main()

# Display results.
include("view.jl")
imshowall(images, beam_energy_map, transmit_time_map, sim_params)

# Save results.
include("save.jl")
saveall(images, beam_energy_map, transmit_time_map, sim_params, "images")




# With a lens:

using WaveSim
using StaticArrays
include("view.jl")
include("save.jl")

function main_lens()
  # Simulate a lens in 2 steps.
  # 1. Simulate wave propagation from transducer in lens material.
  # 2. Record the waveforms at the lens exterior surface.
  # 3. Simulate wave propagation from lens exterior surface in medium.

  # Radius of curvature of the lens.
  lens_radius = 1e-2  # [m]
  # Speed of sound.
  lens_c = 900  # [m/s]

  # 1. Simulate propagation from transducer in lens material.
  sim_params = WaveSimParameters(
    focus_depth = Inf,  # [m]
    steer_angle = 0.0,  # [deg]
    tx_frequency = 4_000_000.0,  # [Hz]
    transducer_pitch = 205e-6,  # [m]
    aperture_size = 1.4e-3,  # [m]
    aperture_radius = Inf,  # [m]
    fov = [4e-3, 8e-3],  # [m]
    c = lens_c,  # [m/s]  SPEED IN THE LENS
    dbrange = 50  # [dB]
  )

  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)
  sim_params = WaveSim.autores(sim_params, trans_delays)

  images = WaveSim.wavesim(trans_delays, sim_params)

  # intermediate display (lens prop)
  beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params)
  imshowall(images, beam_energy_map, transmit_time_map, sim_params)

  # 2. Record the waveforms at the lens exterior surface.
  transducers, tvec = WaveSim.init(trans_delays, sim_params)
  lens_radius
  pulse_shape_waveforms
  waveforms = images[x, y, i_time]
  function pulse_shape_func(t)
    return 1
  end

  # 3. Simulate wave propagation from lens exterior surface in medium.
  sim_params = WaveSimParameters(
    focus_depth = Inf,  # [m]
    steer_angle = 0.0,  # [deg]
    tx_frequency = 4_000_000.0,  # [Hz]
    pulse_shape_func = pulse_shape_func,
    # transducer_pitch = 205e-6,  # [m]
    # aperture_size = 1.4e-3,  # [m]
    aperture_radius = lens_radius,  # [m]
    fov = [4e-2, 8e-2],  # [m]
    c = 1540,  # [m/s]  SPEED IN THE IMAGED MEDIUM
    dbrange = 50  # [dB]
  )

  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)
  sim_params = WaveSim.autores(sim_params, trans_delays)

  images = WaveSim.wavesim(trans_delays, sim_params)

  beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params)
  imshowall(images, beam_energy_map, transmit_time_map, sim_params)
  saveall(images, beam_energy_map, transmit_time_map, sim_params, "images")
end

main_lens()
