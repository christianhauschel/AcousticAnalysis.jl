using HDF5
using WAV

"""
Save a `PressureTimeHistory` to an HDF5 file. If `group` is specified, the
`PressureTimeHistory` will be saved to a group in the HDF5 file. Otherwise, it
will be saved to the root of the file.
"""
function save_h5(pth::AbstractPressureTimeHistory, fname::String; group=nothing, mode="w")
    h5open(fname, mode) do f
        if group !== nothing
            g = create_group(f, group)
        else
            g = f
        end
        write(g, "p", pth.p)
        write(g, "t0", pth.t0)
        write(g, "fs", 1 / pth.dt)
    end
end

"""
Save a vector of `PressureTimeHistory` to an HDF5 file. Each `PressureTimeHistory`
will be saved to a group in the HDF5 file. The group names are specified by
`groups`.
"""
function save_h5(pths::Vector, fname::String, groups::Vector{String})
    h5open(fname, "w") do f
        for (pth, group) in zip(pths, groups)
            save_h5(pth, fname, group=group, mode="cw")
        end
    end
end

"""
Load a `PressureTimeHistory` from an HDF5 file. If `group` is specified, the
`PressureTimeHistory` will be loaded from a group in the HDF5 file. Otherwise,
it will be loaded from the root of the file.
"""
function load_h5(fname; group=nothing)::AbstractPressureTimeHistory
    p = nothing
    t0 = nothing
    dt = nothing
    h5open(fname, "r") do f
        if group !== nothing
            g = f[group]
        else
            g = f
        end
        p = read(g, "p")
        t0 = read(g, "t0")
        dt = 1 / read(g, "fs")
    end
    return PressureTimeHistory(p, dt, t0)
end

"""
Load a vector of `PressureTimeHistory` from an HDF5 file. Each `PressureTimeHistory`
will be loaded from a group in the HDF5 file. The group names are specified by
`groups`.
"""
function load_h5(fname::String, groups::Vector{String})::Vector
    pths = Vector{PressureTimeHistory}(undef, length(groups))
    for (i, group) in enumerate(groups)
        pths[i] = load_h5(fname, group=group)
    end
    return pths
end

"""
Load a `PressureTimeHistory` from a WAV file.

# Arguments
- `fname::String`: Path to the WAV file.
- `calibration_factor::Number=1.0`: Calibration factor to apply to the pressure data. 
- `t0::Number=0.0`: Time offset to apply to the pressure data.
- `format::String="native"`: Format of the WAV file. See `wavread` for details.
- `kwargs...`: Additional keyword arguments to pass to `wavread`.

# Returns 
- `pth`: `PressureTimeHistory`
"""
function load_wav(fname; calibration_factor=1.0, t0=0.0, format="native", kwargs...)
    p, fs, nbits, opt = wavread(fname, format=format, kwargs...)
    p = p[:, 1]
    p *= calibration_factor
    dt = 1 / fs
    return PressureTimeHistory(p, dt, t0)#, Int(nbits), opt
end

"""
Save a pressure vector to a WAV file.
"""
function _save_wav(p::Vector, fs::Int, fname::String; kwargs...)
    wavwrite(p, fname, Fs=fs; kwargs...)
end

"""
Save a `PressureTimeHistory` to a WAV file.
"""
function save_wav(pth::AbstractPressureTimeHistory, fname::String; kwargs...)
    _save_wav(pth.p, Int(round(1.0/pth.dt)), fname; kwargs...)
end