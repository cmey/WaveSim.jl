using Pkg
Pkg.instantiate()

using CairoMakie
using SpecialFunctions
using WaveSim

function exact_directivity_Umchid_2009(theta::Float32, tx_frequency::Float32, c::Float32, transducer_pitch::Float32)::Float32
    crosstalk_factor = 1.2f0
    element_surface_diameter = transducer_pitch * crosstalk_factor
    a = element_surface_diameter / 2.0f0
    lambda = c / tx_frequency
    k = 2.0f0 * Float32(pi) / lambda
    ka_sin_theta = k * a * sin(theta)
    if abs(ka_sin_theta) < 1f-5
        return 1.0f0
    end
    return Float32(2.0f0 * besselj1(ka_sin_theta) / ka_sin_theta)
end

function summarize_case(label::AbstractString, tx_frequency::Float32, c::Float32, transducer_pitch::Float32)
    angles_deg = Float32.(range(-90.0f0, 90.0f0, length=2001))
    angles_rad = angles_deg .* (Float32(pi) / 180.0f0)
    exact = exact_directivity_Umchid_2009.(angles_rad, tx_frequency, c, transducer_pitch)
    approx = WaveSim.default_directivity.(angles_rad, tx_frequency, c, transducer_pitch)
    abs_err = abs.(approx .- exact)
    rel_err = abs_err ./ max.(abs.(exact), eps(Float32))
    x_max = 2.0f0 * Float32(pi) * (tx_frequency / c) * ((transducer_pitch * 1.2f0) / 2.0f0)
    println(label)
    println("  x_max = ", x_max)
    println("  max abs err over ±90° = ", maximum(abs_err))
    println("  max rel err over ±90° = ", maximum(rel_err))
    println("  max rel err over ±60° = ", maximum(rel_err[abs.(angles_deg) .<= 60.0f0]))
    println("  max rel err over ±30° = ", maximum(rel_err[abs.(angles_deg) .<= 30.0f0]))
    return angles_deg, exact, approx, abs_err
end

outdir = joinpath(@__DIR__, "..", "images", "directivity_study")
mkpath(outdir)

configs = [
    ("5.0 MHz, 205 µm pitch", 5.0f6, 1540.0f0, 205f-6),
    ("8.93 MHz, 102 µm pitch", 8.93f6, 1540.0f0, 102f-6),
    ("3.0 MHz, 208 µm pitch", 3.0f6, 1540.0f0, 208f-6),
]

case_data = [summarize_case(label, tx_frequency, c, pitch) for (label, tx_frequency, c, pitch) in configs]

fig = Figure(size = (1500, 1200), fontsize = 20)
for (row, ((label, _, _, _), (angles_deg, exact, approx, abs_err))) in enumerate(zip(configs, case_data))
    ax_directivity = Axis(fig[row, 1], title = label, xlabel = "Angle from normal [deg]", ylabel = "Directivity")
    lines!(ax_directivity, angles_deg, exact, label = "Exact J1", linewidth = 3, color = :black)
    lines!(ax_directivity, angles_deg, approx, label = "Polynomial approx", linewidth = 2, color = :dodgerblue)
    axislegend(ax_directivity, position = :rb)

    ax_error = Axis(fig[row, 2], xlabel = "Angle from normal [deg]", ylabel = "Absolute error", yscale = log10, title = "|approx - exact|")
    lines!(ax_error, angles_deg, max.(abs_err, 1f-12), linewidth = 2, color = :crimson)
end

save(joinpath(outdir, "directivity_vs_angle_abs_error.png"), fig)

x_values = Float32.(range(0f0, 3.0f0, length = 2001))
j1_exact = Float32.(besselj1.(x_values))
j1_approx = WaveSim.besselj1_approx.(x_values)
j1_abs_err = abs.(j1_approx .- j1_exact)

fig2 = Figure(size = (1400, 700), fontsize = 20)
ax1 = Axis(fig2[1, 1], title = "Bessel J1 approximation", xlabel = "x", ylabel = "J1(x)")
lines!(ax1, x_values, j1_exact, label = "Exact", linewidth = 3, color = :black)
lines!(ax1, x_values, j1_approx, label = "Polynomial approx", linewidth = 2, color = :dodgerblue)
axislegend(ax1, position = :rb)

ax2 = Axis(fig2[1, 2], title = "Absolute error", xlabel = "x", ylabel = "|error|", yscale = log10)
lines!(ax2, x_values, max.(j1_abs_err, 1f-12), linewidth = 2, color = :crimson)

save(joinpath(outdir, "besselj1_approx_abs_error.png"), fig2)

println("Saved plots to ", outdir)