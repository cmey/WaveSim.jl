using Pkg
Pkg.instantiate()

using CairoMakie
using BenchmarkTools
using SpecialFunctions
using WaveSim

function exact_directivity(θ::Float32, tx_frequency::Float32, c::Float32, transducer_pitch::Float32)::Float32
    crosstalk_factor = 1.2f0
    element_surface_diameter = transducer_pitch * crosstalk_factor
    a = element_surface_diameter / 2.0f0
    λ = c / tx_frequency
    k = 2.0f0 * Float32(pi) / λ
    ka_sin_θ = k * a * sin(θ)
    if abs(ka_sin_θ) < 1f-5
        return 1.0f0
    end
    return Float32(2.0f0 * besselj1(ka_sin_θ) / ka_sin_θ)
end

function sum_directivity(f, angles, tx_frequency, c, pitch)
    acc = 0.0f0
    @inbounds for θ in angles
        acc += f(θ, tx_frequency, c, pitch)
    end
    return acc
end

function sum_bessel(f, xs)
    acc = 0.0f0
    @inbounds for x in xs
        acc += f(x)
    end
    return acc
end

function report(label::AbstractString, trial)
    estimate = minimum(trial)
    println(label)
    println("  min time = ", estimate.time, " ns")
    println("  allocs   = ", estimate.allocs)
    println("  memory   = ", estimate.memory, " bytes")
    return estimate
end

function benchmark_directivity_case(tx_frequency::Float32, pitch::Float32, angles, c::Float32)
    approx_time = @belapsed sum_directivity($WaveSim.default_directivity, $angles, $tx_frequency, $c, $pitch)
    exact_time = @belapsed sum_directivity($exact_directivity, $angles, $tx_frequency, $c, $pitch)
    return approx_time, exact_time
end

tx_frequency = 5.0f6
c = 1540.0f0
pitch = 205f-6
angles = Float32.(range(0.0f0, stop=Float32(pi) / 2.0f0, length=2001))
xs = Float32.(range(0.0f0, stop=2.5091944f0, length=2001))

sum_directivity(WaveSim.default_directivity, angles, tx_frequency, c, pitch)
sum_directivity(exact_directivity, angles, tx_frequency, c, pitch)
sum_bessel(WaveSim.besselj1_approx, xs)
sum_bessel(besselj1, xs)

approx_trial = @benchmark sum_directivity($WaveSim.default_directivity, $angles, $tx_frequency, $c, $pitch)
exact_trial = @benchmark sum_directivity($exact_directivity, $angles, $tx_frequency, $c, $pitch)
approx_bessel_trial = @benchmark sum_bessel($WaveSim.besselj1_approx, $xs)
exact_bessel_trial = @benchmark sum_bessel($besselj1, $xs)

approx_directivity_estimate = report("approx directivity", approx_trial)
exact_directivity_estimate = report("exact directivity", exact_trial)
approx_bessel_estimate = report("approx besselj1", approx_bessel_trial)
exact_bessel_estimate = report("SpecialFunctions.besselj1", exact_bessel_trial)

println("Directivity benchmark on 2001 angles in [0, 90] degrees")
println("Geometry: 5.0 MHz, 205 um pitch, crosstalk factor 1.2")
println("  speedup = ", round(minimum(exact_trial).time / minimum(approx_trial).time, digits = 3), "x")
println("  speedup = ", round(minimum(exact_bessel_trial).time / minimum(approx_bessel_trial).time, digits = 3), "x")

outdir = joinpath(@__DIR__, "..", "images", "directivity_study")
mkpath(outdir)

fig = Figure(size = (1100, 550), fontsize = 20)

ax1 = Axis(fig[1, 1], title = "Directivity sweep cost", ylabel = "ns / call", xticklabelsvisible = false)
barplot!(ax1, [1, 2], [approx_directivity_estimate.time, exact_directivity_estimate.time],
    color = [:dodgerblue, :black],
    width = 0.6)
text!(ax1, [1, 2], [approx_directivity_estimate.time, exact_directivity_estimate.time],
    text = ["approx", "SpecialFunctions"],
    align = (:center, :bottom),
    offset = (0, 12))

ax2 = Axis(fig[1, 2], title = "Bessel J1 cost", ylabel = "ns / call", xticklabelsvisible = false)
barplot!(ax2, [1, 2], [approx_bessel_estimate.time, exact_bessel_estimate.time],
    color = [:dodgerblue, :black],
    width = 0.6)
text!(ax2, [1, 2], [approx_bessel_estimate.time, exact_bessel_estimate.time],
    text = ["approx", "SpecialFunctions"],
    align = (:center, :bottom),
    offset = (0, 12))

save(joinpath(outdir, "directivity_benchmark.png"), fig)
println("Saved plot to ", joinpath(outdir, "directivity_benchmark.png"))

frequencies_mhz = Float32[3.0f0, 5.0f0, 8.93f0, 12.0f0]
pitches_um = Float32[102.0f0, 150.0f0, 205.0f0, 300.0f0]
approx_ns_per_call = zeros(Float32, length(pitches_um), length(frequencies_mhz))
exact_ns_per_call = zeros(Float32, length(pitches_um), length(frequencies_mhz))
speedup = zeros(Float32, length(pitches_um), length(frequencies_mhz))

for (j, pitch_um) in enumerate(pitches_um)
    pitch_m = pitch_um * 1f-6
    for (i, frequency_mhz) in enumerate(frequencies_mhz)
        tx_frequency_hz = frequency_mhz * 1f6
        approx_time, exact_time = benchmark_directivity_case(tx_frequency_hz, pitch_m, angles, c)
        approx_ns_per_call[j, i] = Float32(approx_time * 1e9 / length(angles))
        exact_ns_per_call[j, i] = Float32(exact_time * 1e9 / length(angles))
        speedup[j, i] = exact_ns_per_call[j, i] / approx_ns_per_call[j, i]
    end
end

fig2 = Figure(size = (1350, 700), fontsize = 18)
title_row = Label(fig2[0, 1:3], "Directivity benchmark sweep", fontsize = 24, tellwidth = false)
_ = title_row

ax3 = Axis(fig2[1, 1], title = "Approx ns/call", xlabel = "Frequency [MHz]", ylabel = "Pitch [µm]", xticks = (1:length(frequencies_mhz), string.(frequencies_mhz)), yticks = (1:length(pitches_um), string.(pitches_um)))
heatmap!(ax3, 1:length(frequencies_mhz), 1:length(pitches_um), approx_ns_per_call, colormap = :viridis)

ax4 = Axis(fig2[1, 2], title = "Exact ns/call", xlabel = "Frequency [MHz]", ylabel = "Pitch [µm]", xticks = (1:length(frequencies_mhz), string.(frequencies_mhz)), yticks = (1:length(pitches_um), string.(pitches_um)))
heatmap!(ax4, 1:length(frequencies_mhz), 1:length(pitches_um), exact_ns_per_call, colormap = :magma)

ax5 = Axis(fig2[1, 3], title = "Speedup exact / approx", xlabel = "Frequency [MHz]", ylabel = "Pitch [µm]", xticks = (1:length(frequencies_mhz), string.(frequencies_mhz)), yticks = (1:length(pitches_um), string.(pitches_um)))
heatmap!(ax5, 1:length(frequencies_mhz), 1:length(pitches_um), speedup, colormap = :turbo)

for j in axes(speedup, 1), i in axes(speedup, 2)
    text!(ax5, i, j, text = string(round(speedup[j, i], digits = 2)), align = (:center, :center), color = :black, fontsize = 14)
end

save(joinpath(outdir, "directivity_speedup_sweep.png"), fig2)
println("Saved plot to ", joinpath(outdir, "directivity_speedup_sweep.png"))
