function narrowband_spectrum(p::Vector, fs::Number; window=:hanning, type=:amplitude, aweighting=false)
    N = length(p)
    if window == :hanning
        window = DSP.hanning(N)
    elseif window === nothing
        window = ones(N)
    end
    p_win = p .* window
    X = fft(p_win)

    n_freq = div(N, 2) + 1
    f_nyquist = fs / 2.0
    f = LinRange(0, f_nyquist, n_freq)
    Δf = f[2] - f[1]

    X_ss = X[1:div(N, 2)+1]
    X_ss *= 2
    X_ss[1] = X_ss[1] / 2.0

    if type == :amplitude
        X_ss = abs.(X_ss) / sum(window)
        X_dB = 20.0 * log10.(X_ss / P_REF)
    else
        X_ss = abs2.(X_ss) / (fs * sum(window .^ 2))
        X_dB = 10.0 * log10.(X_ss * Δf / P_REF^2)
    end

    if aweighting
        X_dB = dB2dBA.(f, X_dB)
        X_dB = clean_data(f, X_dB)
    end

    return f, X_dB, X_ss
end
function narrowband_spectrum(p::PressureTimeHistory; window=:hanning, type=:amplitude, aweighting=false)
    return narrowband_spectrum(p.p, 1.0/p.dt; window=window, type=type, aweighting=aweighting)
end