# This function creates bipolar (i.e. signed) log of the input.

function bilog(datain, dbrange=90)
    datain = real.(datain)
    maxabs = abs.(maximum(datain))
    minabs = abs.(minimum(datain))

    if minabs > maxabs
        maxabs = minabs
    end

    norm_mag = (abs.(datain) + eps()) / maxabs
    out = sign.(datain) .* (clamp.(20 * log10.(norm_mag), -dbrange, 0) + dbrange)

    return out
end
