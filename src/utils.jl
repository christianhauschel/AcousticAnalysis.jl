function dbsum(level1::Float64, level2::Float64)
    return 10 * log10(10^(level1 / 10) + 10^(level2 / 10))
end

function dbsum(levels::Vector)
    return 10 * log10(sum(10 .^ (levels / 10)))
end

function dB(p; prefix = 20, P_REF=20e-6)
    return prefix * log10.(p ./ P_REF)
end

"""
    pressure2power(L_p; r=1, Q=1)

Converts a pressure level to a power level.

# Arguments
- `L_p`: pressure level [dB]
- `r`: radius of the sphere [m]
- `Q`: radiation pattern (sphere: 1, semi-sphere: 2, quarter-sphere: 4, eighth-sphere: 8)

# Literature
[1] https://www.linkedin.com/pulse/acoustics-spl-vs-swl-sound-pressure-level-power-chris-jones/
"""
function pressure2power(L_p; r=1, Q=1)
    return L_p .+ abs.(10 * log10.(Q ./ (4π .* r.^2)))
end