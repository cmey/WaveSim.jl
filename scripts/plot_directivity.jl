# display directivity_func
using CairoMakie
using SpecialFunctions

fc = 6.25f6
c = 1540f0
pitch = 205f-6  # Tip: try with a bigger pitch (e.g. 1205f-6) to see the oscillations
thetas = range(-π/2, stop=π/2, length=100)
thetas_deg = rad2deg.(thetas)

function directivity_func(θ, tx_frequency, c, transducer_pitch)
    # account for mechanical crosstalk between elements
    crosstalk_factor = 1.2f0  # assume 20% crosstalk factor
    element_surface_diameter = transducer_pitch * crosstalk_factor
    # from Umchid 2009 Directivity Pattern Measurement of Ultrasound Transducers
    # https://www.thaiscience.info/Journals/Article/IABE/10892457.pdf
    a = element_surface_diameter / 2
    λ = c / tx_frequency
    k = 2.0f0 * π / λ
    ka_sin_θ = k * a * sin(θ)
    if abs(ka_sin_θ) < 1f-5
        return 1.0f0
    end
    h = convert(Float32, 2 * besselj1(ka_sin_θ) / ka_sin_θ)
    return h
end

D1 = directivity_func.(thetas, fc, c, pitch)

function directivity_func_tasinkevych(θ, tx_frequency, c, transducer_pitch)
    # account for mechanical crosstalk between elements
    crosstalk_factor = 1.2f0  # assume 20% crosstalk factor
    element_surface_diameter = transducer_pitch * crosstalk_factor
    # from Tasinkevych 2010 Element Directivity Influence In The Synthetic Focusing Algorithm For Ultrasound Imaging
    # https://www.ippt.pan.pl/repository/open/o460.pdf
    λ = c / tx_frequency
    p_d_l_s_t = π * element_surface_diameter / λ * sin(θ)
    if abs(p_d_l_s_t) < 1f-5
        return 1.0f0
    end
    h = convert(Float32, sin(p_d_l_s_t) / p_d_l_s_t * cos(θ))
    return h
end

D2 = directivity_func_tasinkevych.(thetas, fc, c, pitch)

fig, ax, plt = plot(thetas_deg, D1, label="Umchid 2009", color=:blue)
plot!(ax, thetas_deg, D2, label="Tasinkevych 2010", color=:red)
ax.xlabel = "Angle from normal [deg]"
ax.xticks = -90:30:90
xlims!(ax, -90, 90)
axislegend(ax)
save("dir.png", fig)
