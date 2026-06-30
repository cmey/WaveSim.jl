# This function creates bipolar (i.e. signed) log of the input. Good for displaying RF data.
function bilog(datain::AbstractArray, dbrange::Number=40)
    maxabs = maximum(abs, datain)
    return bilog(datain, maxabs, dbrange)
end

function bilog(datain::AbstractArray, maxabs::Number, dbrange::Number)
    datain_real = real.(datain)
    norm_mag = (abs.(datain_real) .+ eps(Float32)) / maxabs
    out = sign.(datain_real) .* (clamp.(10 * log10.(norm_mag), -dbrange, 0) .+ dbrange)

    return out
end

# Normalized logarithm (goes from -dbrange to 0). Good for displaying intensity maps.
function normlog(data, dbrange=50)
    max = maximum(data)
    if max == 0
        return fill(-dbrange, size(data))
    end
    normalized = data ./ max
    out = clamp.(10 * log10.(normalized), -dbrange, 0)

    return out
end