import DSP 

"""
Resampling of a signal to a new sampling frequency.

# Returns 
- resampled signal
- `t_resampled`: time vector of resampled signal
"""
function _resample(p::Vector, fs, fs_resample; t0=0.0)
    p_resampled = DSP.Filters.resample(p, fs_resample / fs)
    dt = 1 / fs_resample
    t_resampled = Vector(0:dt:dt*(length(p_resampled)-1)) .+ t0
    return p_resampled, t_resampled
end
function _resample(p::Vector, t::Vector, fs_resample)
    dt = t[2] - t[1]
    fs = 1 / dt
    p_resampled = DSP.Filters.resample(p, fs_resample / fs)

    dt = 1 / fs_resample
    t_resampled = Vector(0:dt:dt*(length(p_resampled)-1))
    return p_resampled, t_resampled
end

"""
Resampling of a pressure time history to a new sampling frequency.

# Returns
- resampled pressure time history
"""
function resample(pth::AbstractPressureTimeHistory, fs_resample)
    p_resampled, t_resampled = _resample(pth.p, 1/pth.dt, fs_resample)
    return PressureTimeHistory(p_resampled, t_resampled[2] - t_resampled[1], t_resampled[1])
end

"""
Selects pressure time history between t_start and t_end.
"""
function time_range(pth::AbstractPressureTimeHistory, t_start, t_end; suppress_warnings=true)
    t = Vector(time(pth))
    p, t = _time_range(pth.p, t, t_start, t_end; suppress_warnings=suppress_warnings)
    return PressureTimeHistory(p, t[2] - t[1], t[1])
end

"""
Selects signal between t_start and t_end.
"""
function _time_range(p::Vector, t::Vector, t_start, t_end; suppress_warnings=true)

    if t_start < t[1]
        t_start = t[1]
        if suppress_warnings == false
            println("Warning: t_start < t[1]. Setting t_start = t[1].")
        end
    end
    if t_end > t[end]
        t_end = t[end]
        if suppress_warnings == false
            println("Warning: t_end > t[end]. Setting t_end = t[end].")
        end
    end

    i0 = findfirst(t -> t >= t_start, t)
    i1 = findfirst(t -> t >= t_end, t)
    return p[i0:i1], t[i0:i1]
end
