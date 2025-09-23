# Wave propagation simulator

Simulates the propagation of waves from multiple transmitting elements.

## Package status

![CI status](https://github.com/cmey/WaveSim.jl/actions/workflows/ci.yml/badge.svg)

## Installation

Until this package gets registered, open the Julia environment from inside the package folder:
```
julia --project=.
```

## Usage

```
using WaveSim

# Define simulation parameters (use many default values, see the WaveSimParameters struct).
sim_params = WaveSimParameters(
    focus_depth = 0.03,  # [m]
    steer_angle = 10.0,  # [deg]
    aperture_size = 0.02,  # [m]
);

# Compute focusing delay for the elements of the phased array.
trans_delays = WaveSim.delays_from_focus_and_steer(sim_params);

# (optional) Compute optimized "best" spatial and temporal parameters.
sim_params = WaveSim.autores(sim_params, trans_delays)

# Run the simulation.
images = WaveSim.wavesim(trans_delays, sim_params);
integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params);

# Display results.
include("src/view.jl")
imshowall(images, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params);

# Save results.
include("src/save.jl")
saveall(images, integrated_energy_map, peak_to_peak_map, transmit_time_map, peak_to_peak_time_delta_map, sim_params, "images")
```

Visualize the wave propagating through space, over time:

![wave propagation animation](images/wave_propagation.gif)

Get a spatial heatmap of where the time-integrated energy went:

![integrated energy map](images/integrated_energy_map.png)

Get a spatial heatmap of the peak-to-peak amplitude:

![peak-to-peak map](images/peak_to_peak_map.png)

Get a spatial heatmap of the transmit time delay to peak amplitude:

![transmit time map](images/transmit_time_map.png)

Get a spatial heatmap of the peak-to-peak time delta (can indicate where the peak-to-peak estimations are reasonable):

![peak-to-peak time delta map](images/peak_to_peak_time_delta_map.png)


## Tips

### Parallelization

The code supports multi-threading, make use of it by starting Julia with multiple threads: `julia --threads 4 --project=.`
