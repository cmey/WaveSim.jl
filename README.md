# Wave propagation simulator

Simulates the propagation of waves from multiple transmitting elements.

## Package status

| macOS | Linux | Windows |
|-------|-------|---------|
|[![Build Status](https://travis-matrix-badges.herokuapp.com/repos/cmey/WaveSim.jl/branches/master/2)](https://travis-ci.org/cmey/WaveSim.jl)|[![Build Status](https://travis-matrix-badges.herokuapp.com/repos/cmey/WaveSim.jl/branches/master/1)](https://travis-ci.org/cmey/WaveSim.jl)|[![Build status](https://ci.appveyor.com/api/projects/status/8pqnoxopn8g8fstv?svg=true)](https://ci.appveyor.com/project/cmey/wavesim-jl)|

## Usage

```
include("src/WaveSim.jl")
using WaveSim

# Define simulation parameters (use many default values, see WaveSimParameters).
sim_params = WaveSimParameters(
    focus_depth = 0.03,  # [m]
    steer_angle = 10.0,  # [deg]
    aperture_size = 0.02,  # [m]
    temporal_res = 0.1e-6,  # [s]
    spatial_res = [128, 256]  # [pixels]
)

# Compute focusing delays for the elements of the phased array.
trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)

# Run the simulation.
images = WaveSim.wavesim(trans_delays, sim_params)
beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params)

# Display results.
include("src/view.jl")
imshowall(images, beam_energy_map, transmit_time_map, sim_params)
```

Visualize the wave propagating through space, over time:

![wave propagation animation](images/wave_propagation.gif)

Get a spatial heatmap of where the transmitted energy is sent:

![beam energy map](images/beam_energy_map.png)

Know the transmit time delay everwhere in space:

![beam energy map](images/transmit_time_map.png)

## Tips

### Parallelization

The code supports multi-threading, make use of it by setting:

`export JULIA_NUM_THREADS=4` (or whatever number of cores your machine has), before starting `julia`

or start `julia` directly with:

    JULIA_NUM_THREADS=`getconf _NPROCESSORS_ONLN` julia
