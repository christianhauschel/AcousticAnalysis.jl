using HDF5
using WAV

function save_h5(pth::PressureTimeHistory, fname::String; group=nothing, mode="w")
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
function save_h5(pths::Vector, fname::String, groups::Vector{String})
    h5open(fname, "w") do f
        for (pth, group) in zip(pths, groups)
            save_h5(pth, fname, group=group, mode="cw")
        end
    end
end

function load_h5(fname; group=nothing)
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
function load_h5(fname::String, groups::Vector{String})
    pths = Vector{PressureTimeHistory}(undef, length(groups))
    for (i, group) in enumerate(groups)
        pths[i] = load_h5(fname, group=group)
    end
    return pths
end

function load_wav(fname; calibration_factor=2.56e-7, t0=0.0)::PressureTimeHistory
    p, fs = wavread(fname, format="native")
    p *= calibration_factor
    dt = 1 / fs
    return PressureTimeHistory(p, dt, t0)
end

save_wav(p::Vector, fs, fname::String) = wavwrite(p, fname, Fs=fs)
function save_wav(pth::PressureTimeHistory, fname::String)
    save_wav(pth.p, 1/pth.dt, fname)
end