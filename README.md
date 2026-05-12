# Wave propagation simulator

Simulates the interference and propagation of waves from multiple transmitting elements.

## Package status

![CI status](https://github.com/cmey/WaveSim.jl/actions/workflows/ci.yml/badge.svg)

## Installation

The package is registered, so one can simply:
```
using Pkg
Pkg.add("WaveSim")
```

For development, open the Julia environment from inside this package folder. It will re-precompile an up-to-date version with any local change:
```
julia --project=.
```

## Usage

```julia
using WaveSim

# Run the simulator and display results.
function main()
  # Define simulation parameters (use many default values, see WaveSimParameters).
  sim_params = WaveSimParameters(
    focus_depth = 0.03,  # [m]
    steer_angle = 10.0,  # [deg]
    tx_frequency = 3_000_000.0,  # [Hz]
    transducer_pitch = 208e-6,  # [m]
    aperture_size = 0.02,  # [m]
    dbrange = 50,  # [dB]
  )

  # Compute focusing delays for the elements of the phased array.
  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)

  # Compute optimized "best" spatial and temporal parameters.
  sim_params = WaveSim.autores(sim_params, trans_delays, multiplier=0.5)

  # Run the simulation.
  images = WaveSim.wavesim(trans_delays, sim_params)
  windowed_energy_map, integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params)

  return windowed_energy_map, integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map, images, sim_params
end

windowed_energy_map, integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map, images, sim_params = main()

# Display results.
imshowall(images, windowed_energy_map, integrated_energy_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params)

# Save results.
saveall(images, windowed_energy_map, integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params, "images")

```

Visualize the wave propagating through space, over time:

![wave propagation animation](images/wave_propagation.gif)

Get a spatial heatmap of where the max-SPL-windowed energy went:

![windowed energy map](images/windowed_energy_map.png)

Get a spatial heatmap of where the time-integrated energy went:

![integrated energy map](images/integrated_energy_map.png)

Get a spatial heatmap of the peak-to-peak amplitude:

![peak-to-peak map](images/peak_to_peak_map.png)

Get a spatial heatmap of the transmit time delay to the peak amplitude:

![transmit time map](images/transmit_time_map.png)

Get a spatial heatmap of the peak-to-peak time delta (can indicate where the peak-to-peak estimations are reasonable):

![peak-to-peak time delta map](images/peak_to_peak_time_delta_map.png)


## Tips

### Parallelization

The code supports multi-threading, make use of it by starting Julia with multiple threads: `julia --threads 4 --project=.`

### CUDA backend

If `CUDA.jl` is available in the active environment (`import CUDA`) and a functional
NVIDIA GPU is present, the simulation can be run on the GPU with:
`WaveSim.wavesim(trans_delays, sim_params; backend = :cuda)`.
A comparison script is available at `scripts/benchmark_cuda.jl`.
