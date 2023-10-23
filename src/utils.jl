function dbsum(level1::Float64, level2::Float64)
    return 10 * log10(10^(level1 / 10) + 10^(level2 / 10))
end

function dbsum(levels::Vector)
    return 10 * log10(sum(10 .^ (levels / 10)))
end

function dB(p; prefix = 20, P_REF=20e-6)
    return prefix * log10.(p ./ P_REF)
end