using Test
using InteractiveUtils
using StaticArrays
using SpecialFunctions

try
	using CUDA
catch
end

using WaveSim

@testset "slice helpers" begin
	@test WaveSim.slice_lateral_label(:azimuth_depth) == "Azimuth"
	@test WaveSim.slice_lateral_label(:elevation_depth) == "Elevation"
	@test WaveSim.slice_pixel_coordinates(:azimuth_depth, 1.0f0, 2.0f0) == (1.0f0, 0.0f0, 2.0f0)
	@test WaveSim.slice_pixel_coordinates(:elevation_depth, 1.0f0, 2.0f0) == (0.0f0, 1.0f0, 2.0f0)
	@test WaveSim.slice_directivity_pitch(:azimuth_depth, 205f-6, 1400f-6) == 205f-6
	@test WaveSim.slice_directivity_pitch(:elevation_depth, 205f-6, 1400f-6) == 1400f-6
	original_polygons = [(
		SVector{2,Float32}(1.0f0, 2.0f0),
		SVector{2,Float32}(3.0f0, 4.0f0),
		SVector{2,Float32}(5.0f0, 6.0f0),
		SVector{2,Float32}(7.0f0, 8.0f0),
	)]
	@test WaveSim.orient_beamplot_polygons(original_polygons, false) == original_polygons
	@test WaveSim.orient_beamplot_polygons(original_polygons, true) == [(
		SVector{2,Float32}(2.0f0, 1.0f0),
		SVector{2,Float32}(4.0f0, 3.0f0),
		SVector{2,Float32}(6.0f0, 5.0f0),
		SVector{2,Float32}(8.0f0, 7.0f0),
	)]
	flat_params = WaveSimParameters(
		aperture_size = 0.004f0,
		aperture_size_elevation = 0.0f0,
		transducer_pitch = 0.001f0,
		transducer_pitch_elevation = 0.0015f0,
		aperture_radius = Inf,
		aperture_radius_elevation = Inf,
		beamplot_axes = :azimuth_depth,
	)
	azimuth_polygons = WaveSim.beamplot_element_polygons(flat_params)
	@test length(azimuth_polygons) == 4
	@test all(minimum(corner[1] for corner in polygon) >= -1f-7 for polygon in azimuth_polygons)
	@test all(isapprox(sum(corner[1] for corner in polygon) / 4, 0.0005f0; atol=1f-7) for polygon in azimuth_polygons)

	elevation_params = WaveSimParameters(
		aperture_size = 0.001f0,
		aperture_size_elevation = 0.003f0,
		transducer_pitch = 0.001f0,
		transducer_pitch_elevation = 0.001f0,
		aperture_radius = Inf,
		aperture_radius_elevation = Inf,
		beamplot_axes = :elevation_depth,
	)
	elevation_polygons = WaveSim.beamplot_element_polygons(elevation_params)
	@test length(elevation_polygons) == 3
end

@testset "default directivity" begin
	x_values = Float32.(range(0f0, 2.6f0, length=11))
	for x in x_values
		expected = Float32(besselj1(x))
		measured = WaveSim.besselj1_approx(x)
		@test isapprox(measured, expected; atol=1f-4, rtol=5f-4)
	end

	tx_frequency = 5.0f6
	c = 1540.0f0
	transducer_pitch = 205f-6
	for theta_deg in (-90.0f0, -60.0f0, -30.0f0, -15.0f0, 0.0f0, 15.0f0, 30.0f0, 60.0f0, 90.0f0)
		theta = theta_deg * Float32(pi) / 180.0f0
		crosstalk_factor = 1.2f0
		a = (transducer_pitch * crosstalk_factor) / 2.0f0
		lambda = c / tx_frequency
		k = 2.0f0 * Float32(pi) / lambda
		x = k * a * sin(theta)
		expected = abs(x) < 1f-5 ? 1.0f0 : Float32(2.0f0 * besselj1(x) / x)
		measured = WaveSim.default_directivity(theta, tx_frequency, c, transducer_pitch)
		@test isapprox(measured, expected; atol=5f-5, rtol=5f-4)
		@test abs(measured) <= 1.0f0
	end
end

@testset "forward hemisphere gating" begin
	image_forward = zeros(Float32, 1, 1)
	WaveSim.simulate_one_time_step!(
		image_forward,
		1.0f-3 / 1540.0f0,
		Float32[0.0f0],
		1.0f-6,
		1.0f6,
		1540.0f0,
		SVector{2,Int}(1, 1),
		_ -> 1.0f0,
		Float32[1.0f0],
		(θ, tx_frequency, c, transducer_pitch) -> 1.0f0,
		205f-6,
		0.0f0,
		Float32[0.0f0],
		Float32[1.0f-3],
		[SVector{3,Float32}(0.0f0, 0.0f0, 0.0f0)],
		[SVector{3,Float32}(0.0f0, 0.0f0, 1.0f0)],
		:azimuth_depth,
		[1],
	)
	@test image_forward[1, 1] > 0.0f0

	image_behind = zeros(Float32, 1, 1)
	WaveSim.simulate_one_time_step!(
		image_behind,
		1.0f-3 / 1540.0f0,
		Float32[0.0f0],
		1.0f-6,
		1.0f6,
		1540.0f0,
		SVector{2,Int}(1, 1),
		_ -> 1.0f0,
		Float32[1.0f0],
		(θ, tx_frequency, c, transducer_pitch) -> 1.0f0,
		205f-6,
		0.0f0,
		Float32[0.0f0],
		Float32[0.0f0],
		[SVector{3,Float32}(0.0f0, 0.0f0, 1.0f-3)],
		[SVector{3,Float32}(0.0f0, 0.0f0, 1.0f0)],
		:azimuth_depth,
		[1],
	)
	@test image_behind[1, 1] == 0.0f0
end

@testset "delay synthesis" begin
	axis_coords = Float32[-0.002f0, 0.0f0, 0.002f0]
	diverging_delays = WaveSim.axis_delay_law(axis_coords, -0.03f0, 45.0f0, 1540.0f0)
	expected_diverging = Float32[
		sqrt((coord + sind(45.0f0) * 0.03f0)^2 + (cosd(45.0f0) * 0.03f0)^2) / 1540.0f0 for coord in axis_coords
	]
	expected_diverging .-= minimum(expected_diverging)
	@test all(isapprox.(diverging_delays, expected_diverging; atol=1f-7, rtol=1f-6))

	focused_params = WaveSimParameters(
		tx_frequency = 3.0f6,
		pulse_cycles = 1.0f0,
		c = 1540.0f0,
		end_simulation_time = 4e-6,
		temporal_res = 1e-6,
		fov = SVector{2,Float32}(0.01f0, 0.01f0),
		spatial_res = SVector{2,Int}(2, 2),
		transducer_pitch = 0.001f0,
		transducer_pitch_elevation = 0.001f0,
		aperture_size = 0.004f0,
		aperture_size_elevation = 0.002f0,
		aperture_radius = Inf,
		aperture_radius_elevation = Inf,
		focus_depth = 0.03f0,
		focus_depth_elevation = 0.02f0,
		steer_angle = 10.0f0,
		steer_angle_elevation = -5.0f0,
		directivity_func = (θ, tx_frequency, c, transducer_pitch) -> 1.0f0
	)

	focused_delays = WaveSim.delays_from_focus_and_steer(focused_params)
	@test length(focused_delays) == 8
	@test minimum(focused_delays) == 0.0f0
	@test all(isfinite, focused_delays)
	@test length(unique(round.(focused_delays; digits=7))) > 1

	plane_params = WaveSimParameters(
		tx_frequency = 3.0f6,
		pulse_cycles = 1.0f0,
		c = 1540.0f0,
		end_simulation_time = 2e-6,
		temporal_res = 1e-6,
		fov = SVector{2,Float32}(0.01f0, 0.01f0),
		spatial_res = SVector{2,Int}(2, 2),
		transducer_pitch = 0.001f0,
		transducer_pitch_elevation = 0.001f0,
		aperture_size = 0.004f0,
		aperture_size_elevation = 0.0f0,
		aperture_radius = Inf,
		aperture_radius_elevation = Inf,
		focus_depth = Inf,
		focus_depth_elevation = Inf,
		steer_angle = 15.0f0,
		steer_angle_elevation = 0.0f0,
		directivity_func = (θ, tx_frequency, c, transducer_pitch) -> 1.0f0
	)

	plane_delays = WaveSim.delays_from_focus_and_steer(plane_params)
	@test length(plane_delays) == 4
	@test issorted(plane_delays)
end

@testset "element sensitivities" begin
	# 1D array - Vector input
	# 2 elements in azimuth, 1 in elevation
	params_1d = WaveSimParameters(
		aperture_size = 0.01f0,
		transducer_pitch = 0.005f0,
		aperture_size_elevation = 0.0f0,
		element_sensitivities = Float32[0.5, 2.0]
	)
	delays_1d = WaveSim.delays_from_focus_and_steer(params_1d)
	_, _, _, _, apodization_vec, _, _, _ = WaveSim.init(delays_1d, params_1d)
	@test length(apodization_vec) == 2
	@test apodization_vec == Float32[0.5, 2.0]

	# 2D array - Matrix input
	# 2x2 elements
	params_2d = WaveSimParameters(
		aperture_size = 0.01f0,
		aperture_size_elevation = 0.01f0,
		transducer_pitch = 0.005f0,
		transducer_pitch_elevation = 0.005f0,
		element_sensitivities = Float32[1.0 2.0; 3.0 4.0] # (az=2, el=2)
	)
	delays_2d = WaveSim.delays_from_focus_and_steer(params_2d)
	_, _, _, _, apodization_vec_2d, _, _, _ = WaveSim.init(delays_2d, params_2d)
	@test length(apodization_vec_2d) == 4
	# Internal layout is (el, az), flattened as el1-az1, el2-az1, el1-az2, el2-az2
	# For input [1.0 2.0; 3.0 4.0] (az1=[1.0, 3.0], az2=[2.0, 4.0]):
	# el1-az1 = 1.0, el2-az1 = 3.0? No, az1=[1,3], el1 is first row.
	# Input matrix is M[az, el]. az1 = [1, 2], az2 = [3, 4] (if 2 columns)
	# Wait, [1.0 2.0; 3.0 4.0] in Julia:
	# Row 1: 1.0, 2.0
	# Row 2: 3.0, 4.0
	# Columns: col1=[1.0, 3.0], col2=[2.0, 4.0]
	# If az=cols, el=rows: az1 row is [1,2], az2 row is [3,4]? No.
	# Size is (2,2). M[1,1]=1, M[2,1]=3, M[1,2]=2, M[2,2]=4.
	# Az=dim1: az1=[1,2], az2=[3,4]? No, az1 is index 1 of dim 1.
	# M[1,1] is az1, el1. M[2,1] is az2, el1. M[1,2] is az1, el2. M[2,2] is az2, el2.
	# Output vec in (el, az) order: (el1, az1), (el2, az1), (el1, az2), (el2, az2)
	# (1,1), (2,1), (1,2), (2,2) -> 1.0, 2.0, 3.0, 4.0
	@test apodization_vec_2d == Float32[1.0, 2.0, 3.0, 4.0]

	# error case: wrong dimensions
	params_err = WaveSimParameters(
		aperture_size = 0.01f0,
		transducer_pitch = 0.005f0,
		element_sensitivities = Float32[1.0, 2.0, 3.0]
	)
	delays_err = WaveSim.delays_from_focus_and_steer(params_err)
	@test_throws ErrorException WaveSim.init(delays_err, params_err)
end

@testset "small simulation" begin
	sim_params = WaveSimParameters(
		tx_frequency = 1.0f6,
		pulse_cycles = 1.0f0,
		c = 1540.0f0,
		end_simulation_time = 1e-6,
		temporal_res = 1e-6,
		fov = SVector{2,Float32}(0.004f0, 0.004f0),
		spatial_res = SVector{2,Int}(2, 2),
		transducer_pitch = 0.001f0,
		transducer_pitch_elevation = 0.001f0,
		aperture_size = 0.001f0,
		aperture_size_elevation = 0.001f0,
		aperture_radius = Inf,
		aperture_radius_elevation = Inf,
		focus_depth = Inf,
		focus_depth_elevation = Inf,
		steer_angle = 0.0f0,
		steer_angle_elevation = 0.0f0,
		directivity_func = (θ, tx_frequency, c, transducer_pitch) -> 1.0f0
	)

	trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)
	images = WaveSim.wavesim(trans_delays, sim_params)
	@test size(images) == (2, 2, 2)
end

@testset "cuda backend" begin
	if !WaveSim.cuda_backend_available()
		@test_skip false
	else
		sim_params = WaveSimParameters(
			focus_depth = 0.03,
			steer_angle = 10.0,
			aperture_size = 0.01,
			temporal_res = 0.5e-6,
			spatial_res = @SVector [32, 48]
		)

		trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)
		cpu_images = WaveSim.wavesim(trans_delays, sim_params; backend = :cpu)
		gpu_images = WaveSim.wavesim(trans_delays, sim_params; backend = :cuda)

		@test size(cpu_images) == size(gpu_images)
		@test isapprox(cpu_images, gpu_images; atol = 5f-4, rtol = 5f-3)

		cpu_elapsed = @elapsed WaveSim.wavesim(trans_delays, sim_params; backend = :cpu)
		gpu_elapsed = @elapsed WaveSim.wavesim(trans_delays, sim_params; backend = :cuda)
		@test cpu_elapsed > 0.0
		@test gpu_elapsed > 0.0
	end
end

@testset "warntype" begin
	sim_params = WaveSimParameters(
		focus_depth = 0.03,
		steer_angle = 10.0,
		aperture_size = 0.01,
		temporal_res = 0.5e-6,
		spatial_res = @SVector [64, 128]
	)

	trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)
	transducers, tvec, pulse_length, apodization_vec, pix_xs, pix_ys, firing = WaveSim.init(trans_delays, sim_params)
	image = zeros(Float32, (sim_params.spatial_res[1], sim_params.spatial_res[2]))
	t = tvec[10]

	warntype_output = sprint() do io
		InteractiveUtils.code_warntype(io, WaveSim.simulate_one_time_step!, Tuple{typeof(image), typeof(t), typeof(trans_delays), typeof(pulse_length), typeof(sim_params.tx_frequency), typeof(sim_params.c), typeof(sim_params.spatial_res), typeof(sim_params.pulse_shape_func), typeof(apodization_vec), typeof(sim_params.directivity_func), typeof(sim_params.transducer_pitch), typeof(sim_params.attenuation_coefficient), typeof(pix_xs), typeof(pix_ys), typeof(transducers), typeof(sim_params.beamplot_axes), typeof(firing)})
	end

	@test !occursin("::Any", warntype_output)
	@test !occursin("Union{", warntype_output)

	elevation_sim_params = WaveSimParameters(
		focus_depth = 0.03,
		steer_angle = 10.0,
		aperture_size = 0.01,
		aperture_size_elevation = 0.01,
		beamplot_axes = :elevation_depth,
		temporal_res = 0.5e-6,
		spatial_res = @SVector [64, 128]
	)

	elevation_trans_delays = WaveSim.delays_from_focus_and_steer(elevation_sim_params)
	elevation_transducers, elevation_tvec, elevation_pulse_length, elevation_apodization_vec, elevation_pix_xs, elevation_pix_ys, elevation_firing = WaveSim.init(elevation_trans_delays, elevation_sim_params)
	elevation_image = zeros(Float32, (elevation_sim_params.spatial_res[1], elevation_sim_params.spatial_res[2]))
	elevation_t = elevation_tvec[10]

	elevation_warntype_output = sprint() do io
		InteractiveUtils.code_warntype(io, WaveSim.simulate_one_time_step!, Tuple{typeof(elevation_image), typeof(elevation_t), typeof(elevation_trans_delays), typeof(elevation_pulse_length), typeof(elevation_sim_params.tx_frequency), typeof(elevation_sim_params.c), typeof(elevation_sim_params.spatial_res), typeof(elevation_sim_params.pulse_shape_func), typeof(elevation_apodization_vec), typeof(elevation_sim_params.directivity_func), typeof(elevation_sim_params.transducer_pitch), typeof(elevation_sim_params.attenuation_coefficient), typeof(elevation_pix_xs), typeof(elevation_pix_ys), typeof(elevation_transducers), typeof(elevation_sim_params.beamplot_axes), typeof(elevation_firing)})
	end

	@test !occursin("::Any", elevation_warntype_output)
	@test !occursin("Union{", elevation_warntype_output)
end
