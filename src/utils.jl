function dbsum(level1::Float64, level2::Float64)
    return 10 * log10(10^(level1 / 10) + 10^(level2 / 10))
end

function dbsum(levels::Vector)
    return 10 * log10(sum(10 .^ (levels / 10)))
end

function dB(p; prefix=20, P_REF=20e-6)
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
    return L_p .+ abs.(10 * log10.(Q ./ (4π .* r .^ 2)))
end


"""
Remove the DC offset from a vector.
"""
function remove_dc_offset(p::Vector)
    return p .- mean(p)
end
function remove_dc_offset(pth::PressureTimeHistory)
    return PressureTimeHistory(remove_dc_offset(pth.p), timestep(pth), pth.t0)
end

"""
Normalize the pressure time history to -1 or 1
"""
function normalize(p::Vector{Float64}; offset_peak_dB=0.0)
    p_max = maximum(p)
    p_min = minimum(p)
    factor = 1.0
    if abs(p_max) > abs(p_min)
        factor = 1 / p_max
    else
        factor = abs(1 / p_min)
    end
    factor *= 10^(offset_peak_dB / 10)
    return p .* factor, factor
end
function normalize(pth::PressureTimeHistory; offset_peak_dB=0.0)
    p, factor = normalize(pth.p; offset_peak_dB=offset_peak_dB)
    return PressureTimeHistory(p, timestep(pth), pth.t0), factor
end
"""
Normalizes a vector of pressure time histories to -1 or 1.

Uses the minimal factor of all pressures, such that the signals are comparable.
"""
function normalize(ps::Vector{Vector{Float64}}; offset_peak_dB=0.0)
    # calculate factor for each vector
    factors = [AcousticAnalysis.normalize(p; offset_peak_dB=offset_peak_dB)[2] for p in ps]

    # find the min factor
    min_factor = minimum(factors)

    # normalize based on min factor 
    return [p .* min_factor for p in ps]
end
function normalize(pths::Vector{PressureTimeHistory{true, Vector{Float64}, Float64, Float64}}; offset_peak_dB=0.0)
    ps = [pth.p for pth in pths]
    ps = normalize(ps; offset_peak_dB=offset_peak_dB)
    return [PressureTimeHistory(p, timestep(pths[1]), pths[1].t0) for p in ps]
end
function normalize(pths::Vector{PressureTimeHistory{false, Vector{Float64}, Float64, Float64}}; offset_peak_dB=0.0)
    ps = [pth.p for pth in pths]
    ps = normalize(ps; offset_peak_dB=offset_peak_dB)
    return [PressureTimeHistory(p, timestep(pths[1]), pths[1].t0) for p in ps]
end
function normalize(pths::Vector{PressureTimeHistory}; offset_peak_dB=0.0)
    ps = [pth.p for pth in pths]
    ps = normalize(ps; offset_peak_dB=offset_peak_dB)
    return [PressureTimeHistory(p, timestep(pths[1]), pths[1].t0) for p in ps]
end