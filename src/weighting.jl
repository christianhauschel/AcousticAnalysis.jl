function dB2dBA(f, dB)
    h(a) = a^2.0 + f.^2.0
    R = (12200.0^2.0 * f.^4.0) / (h(20.6) * sqrt(h(107.7) * h(737.9)) * h(12200))
    weighting = 2.0 + 20.0 * log10.(R)
    return dB + weighting
end

"""
Compute the A-weighted time history of a pressure time history.
"""
function aweighting(pth::AbstractPressureTimeHistory)::AbstractPressureTimeHistory
    return _aweighting(pth.p, 1/pth.dt)
end

"""
Compute the A-weighted time history of a pressure vector.
"""
function _aweighting(p::Vector, fs)::Vector
    wa = pyimport("waveform_analysis") # TODO: Julia solution
    return wa.weighting_filters.A_weight(p, fs)
end