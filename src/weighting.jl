function dB2dBA(f, dB)
    h(a) = a^2.0 + f.^2.0
    R = (12200.0^2.0 * f.^4.0) / (h(20.6) * sqrt(h(107.7) * h(737.9)) * h(12200))
    weighting = 2.0 + 20.0 * log10.(R)
    return dB + weighting
end

"""
    aweighting(pth)

Compute the A-weighted time history of a pressure time history.
"""
function aweighting(pth::PressureTimeHistory)::PressureTimeHistory
    p = pressure(pth)
    fs = 1 / timestep(pth)

    wa = pyimport("waveform_analysis") # TODO: Julia solution

    pA = wa.weighting_filters.A_weight(p, fs)
    return PressureTimeHistory(pA, 1/fs)
end
