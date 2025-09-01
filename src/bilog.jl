# This function creates bipolar (i.e. signed) log of the input. Good for displaying RF data.
function bilog(datain, dbrange=40)
    datain = real.(datain)
    maxabs = abs.(maximum(datain))
    minabs = abs.(minimum(datain))

    if minabs > maxabs
        maxabs = minabs
    end

    norm_mag = (abs.(datain) .+ eps()) / maxabs
    out = sign.(datain) .* (clamp.(20 * log10.(norm_mag), -dbrange, 0) .+ dbrange)

    return out
end

# Normalized logarithm (goes from -dbrange to 0). Good for displaying intensity maps.
function normlog(data, dbrange=50)
    max = maximum(data)
    if max == 0
        return fill(-dbrange, size(data))
    end
    normalized = data ./ max
    out = clamp.(20 * log10.(normalized), -dbrange, 0)

    return out
end