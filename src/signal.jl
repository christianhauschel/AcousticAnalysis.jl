import DSP 

"""
Resampling of a signal/PressureTimeHistory to a new sampling frequency.

# Returns 
- resampled signal/PressureTimeHistory
- `t_resampled`: time vector of resampled signal
"""
function resample(p::Vector, fs, fs_resample; t0=0.0)
    p_resampled = DSP.Filters.resample(p, fs_resample / fs)
    dt = 1 / fs_resample
    t_resampled = Vector(0:dt:dt*(length(p_resampled)-1)) .+ t0
    return p_resampled, t_resampled
end
function resample(p::Vector, t::Vector, fs_resample)
    dt = t[2] - t[1]
    fs = 1 / dt
    p_resampled = DSP.Filters.resample(p, fs_resample / fs)

    dt = 1 / fs_resample
    t_resampled = Vector(0:dt:dt*(length(p_resampled)-1))
    return p_resampled, t_resampled
end
function resample(pth::PressureTimeHistory, fs_resample)
    p_resampled, t_resampled = resample(pth.p, 1/pth.dt, fs_resample)
    return PressureTimeHistory(p_resampled, t_resampled[2] - t_resampled[1], t_resampled[1])
end

"""
Selects signal/PressureTimeHistory between t_start and t_end.
"""
function time_range(pth::PressureTimeHistory, t_start, t_end; suppress_warnings=true)
    t = Vector(time(pth))
    p, t = time_range(pth.p, t, t_start, t_end; suppress_warnings=suppress_warnings)
    return PressureTimeHistory(p, t[2] - t[1], t[1])
end
function time_range(p::Vector, t::Vector, t_start, t_end; suppress_warnings=true)

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


# t = Vector(LinRange(0, 1, 1000))
# dt = t[2] - t[1]
# fs = 1 / dt
# p = sin.(2π*10*t)

# fs_resample = 25
# p_resampled, t_resampled = resample(p, t, fs_resample)

# pth = PressureTimeHistory(p, dt, 0.0)

# pth_selected = time_range(pth, 0.2, 1);

# # pth_resampled = resample(pth, fs_resample)

# plot_history([pth, pth_selected])