
"""
    clean_data(x, y)

Clean data by removing NaN and Inf values and interpolating linearly.

# Arguments
- `x`: x data
- `y`: y data

# Returns
- `y`: cleaned y data
"""
function clean_data(x, y)
    mask = .!isinf.(y) .& .!isnan.(y)
    _y = y[mask]
    _x = x[mask]

    return linear(_x, _y, x)
end

function plot_spectrogram(
    pth::PressureTimeHistory;
    fname=nothing,
    dpi=300,
    figsize=(6, 5),
    t_window=0.05,
    onesided=true,
    title="Spectrogram",
    ylim=nothing
)
    p = pressure(pth)
    fs = 1 / timestep(pth)

    pplt = pyimport("proplot")
    spectro = spectrogram(p / P_REF, Int(round(fs * t_window)); fs=fs, onesided=onesided)


    fig, ax = pplt.subplots(figsize=figsize)
    cm = ax[1].pcolormesh(spectro.time, spectro.freq, 10 * log10.(spectro.power))
    cb = plt.colorbar(cm, label="Power [dB]")
    if ylim !== nothing
        ax[1].set(ylim=ylim)
    else
        ax[1].set(ylim=(0, fs / 2))
    end
    ax[1].set(
        xlabel="t [s]",
        ylabel="f [Hz]",
        title=title,
    )
    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end
    return fig
end

function plot_narrowband_spectrum(
    pth::PressureTimeHistory;
    fname=nothing,
    dpi=300,
    figsize=(6, 3),
    bpf=nothing,
    n_bpf=10,
    xmin=nothing,
    xmax=nothing,
    ymin=0,
    ymax=nothing,
    title="Spectrum",
    xlabel="f [Hz]",
    ylabel=nothing,
    type=:amplitude,
    window=:hanning,
    label=nothing,
    aweighting=false,
    show_oaspl=true,
    lw=1,
    alpha=1,
)
    p = pressure(pth)
    fs = 1 / timestep(pth)

    if xmin === nothing
        try
            xmin = minimum([100, minimum(bpf) * 0.5])
        catch
            xmin = 100
        end
    end
    if xmax === nothing
        xmax = fs / 2.0
    end

    if show_oaspl
        if label === nothing
            label = ""
        else
            label *= ", "
        end
        label *= f"{OASPL(pth):0.2f} dB"
    end

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
        X_ss_amp = abs.(X_ss) / sum(window)
        X_dB = 20.0 * log10.(X_ss_amp / P_REF)
    else
        X_psd = abs2.(X_ss) / (fs * sum(window .^ 2))
        X_dB = 10.0 * log10.(X_psd * Δf / P_REF^2)
    end

    if aweighting
        X_dB = dB2dBA.(f, X_dB)
        X_dB = clean_data(f, X_dB)
    end


    pplt = pyimport("proplot")

    fig, ax = pplt.subplots(figsize=figsize)

    if label === nothing
        ax[1].plot(f, X_dB, lw=lw, alpha=alpha)
    else
        ax[1].plot(f, X_dB, label=label, lw=lw, alpha=alpha)
    end

    # select X_dB > xmin, < xmax
    _X_dB = X_dB[(f.>xmin).&(f.<xmax)]
    if ymax === nothing

        ymax = maximum(_X_dB) + 2
    end
    if ymin === nothing
        ymin = minimum(_X_dB) - 2
    end

    if bpf !== nothing
        for i = 1:n_bpf
            ax[1].vlines(bpf * i, ymin, ymax, color="k", lw=0.5)
        end
    else
        bpf = 62.5
    end

    if ylabel === nothing
        if type == :amplitude
            ylabel = L"$L_p$ (Amp) [dB re 20 μPa]"
        else
            ylabel = L"$L_p$ [dB re 20 μPa / $Δf^{1/2}_{ref}$]"
        end
    end



    ax[1].set(
        title=title,
        xscale="log",
        xlabel=xlabel,
        ylabel=ylabel,
        ylim=(ymin, ymax),
        xlim=(xmin, xmax),
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)

    end

    return fig
end
function plot_narrowband_spectrum(
    list_pth::Vector;
    fname=nothing,
    dpi=300,
    figsize=(6, 3),
    bpf=[],
    n_bpf=10,
    xmin=nothing,
    xmax=nothing,
    ymin=0,
    ymax=nothing,
    title="Spectrum",
    xlabel="f [Hz]",
    ylabel=nothing,
    type=:amplitude,
    window=:hanning,
    label=nothing,
    aweighting=false,
    lw=1,
    alpha=1,
)
    pplt = pyimport("proplot")

    fig, ax = pplt.subplots(figsize=figsize)

    list_ymax = []
    list_ymin = []
    list_fs = []

    for k in 1:length(list_pth)
        fs = 1 / timestep(list_pth[k])
        push!(list_fs, fs)
    end

    if length(bpf) > 0
        plot_bpf = true
    else
        plot_bpf = false
        push!(bpf, 62.5)
    end


    if xmin === nothing
        xmin = minimum([100, minimum(bpf) * 0.5])
    end
    if xmax === nothing
        xmax = maximum(list_fs) / 2.0
    end


    for k in 1:length(list_pth)
        p = pressure(list_pth[k])
        fs = 1 / timestep(list_pth[k])



        N = length(p)


        if window == :hanning
            _window = DSP.hanning(N)
        elseif window === nothing
            _window = ones(N)
        end
        p_win = p .* _window
        X = fft(p_win)

        n_freq = div(N, 2) + 1
        f_nyquist = fs / 2.0
        f = LinRange(0, f_nyquist, n_freq)
        Δf = f[2] - f[1]

        X_ss = X[1:div(N, 2)+1]
        X_ss *= 2
        X_ss[1] = X_ss[1] / 2.0

        if type == :amplitude
            X_ss_amp = abs.(X_ss) / sum(_window)
            X_dB = 20.0 * log10.(X_ss_amp / P_REF)
        else
            X_psd = abs2.(X_ss) / (fs * sum(_window .^ 2))
            X_dB = 10.0 * log10.(X_psd * Δf / P_REF^2)
        end

        if aweighting
            X_dB = dB2dBA.(f, X_dB)
            X_dB = clean_data(f, X_dB)
        end

        if label === nothing
            ax[1].plot(f, X_dB, lw=lw, alpha=alpha)
        else
            ax[1].plot(f, X_dB, label=label[k], lw=lw, alpha=alpha)
        end

        _X_dB = X_dB[(f.>xmin).&(f.<xmax)]
        if ymax === nothing
            push!(list_ymax, maximum(_X_dB) + 2)
        end
        if ymin === nothing
            push!(list_ymin, minimum(_X_dB) - 2)
        end
    end


    if ymax === nothing
        ymax = maximum(list_ymax)
    end
    if ymin === nothing
        ymin = minimum(list_ymin)
    end

    if plot_bpf
        for k in 1:length(list_pth)
            for i = 1:n_bpf
                ax[1].vlines(bpf[k] * i, ymin, ymax, color=f"C{k-1}", lw=0.5)
            end
        end
    end


    if ylabel === nothing
        if type == :amplitude
            ylabel = L"$L_p$ (Amp) [dB re 20 μPa]"
        else
            ylabel = L"$L_p$ [dB re 20 μPa / $Δf^{1/2}_{ref}$]"
        end
    end

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    list_fs = [1 / timestep(pth) for pth in list_pth]

    ax[1].set(
        title=title,
        xscale="log",
        xlabel=xlabel,
        ylabel=ylabel,
        ylim=(ymin, ymax),
        xlim=(xmin, xmax),
    )

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)

    end

    return fig
end

function plot_history(
    list_pth::Vector;
    title="Time History",
    fname=nothing,
    dpi=300,
    figsize=(6, 2),
    xlabel="t [s]",
    ylabel="p [Pa]",
    label=nothing,
    lw=1,
    alpha=1,
    without_dc_offset=false,
)
    pplt = pyimport("proplot")

    if without_dc_offset
        list_pth = remove_dc_offset.(list_pth)
    end

    fig, ax = pplt.subplots(figsize=figsize)

    for i = 1:length(list_pth)
        p = pressure(list_pth[i])
        t = Vector(time(list_pth[i]))

        if label === nothing
            ax[1].plot(t, p, lw=lw, alpha=alpha)
        else
            ax[1].plot(t, p, label=label[i], lw=lw, alpha=alpha)
        end
    end

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    ax[1].set(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )
    fig

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end

    return fig
end
function plot_history(
    pth::PressureTimeHistory;
    t=nothing,
    title="Time History",
    fname=nothing,
    dpi=300,
    figsize=(6, 2),
    xlabel="t [s]",
    ylabel="p [Pa]",
    label=nothing,
    lw=1,
    alpha=1,
    without_dc_offset=false,
)
    if without_dc_offset
        pth = remove_dc_offset(pth)
    end

    p = pressure(pth)
    fs = 1 / timestep(pth)

    pplt = pyimport("proplot")

    fig, ax = pplt.subplots(figsize=figsize)

    N = length(p)

    Δt = 1.0 / fs
    if t === nothing
        t = Vector(Δt * (0:N-1))
    end

    if label === nothing
        ax[1].plot(t, p, lw=lw, alpha=alpha)
    else
        ax[1].plot(t, p, label=label, lw=lw, alpha=alpha)
    end


    ax[1].set(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end

    return fig
end

function plot_proportional_spectrum(
    pth::PressureTimeHistory;
    n=3,
    aweighting=false,
    figsize=(6, 3),
    fname=nothing,
    xmin=62.5,
    xmax=nothing,
    ymin=nothing,
    ymax=nothing,
    title=nothing,
    ylabel=nothing,
    xlabel=nothing,
    label=nothing,
    dpi=300,
    show_oaspl=true,
    lw=1,
    alpha=1,
)
    if show_oaspl
        if label === nothing
            label = ""
        else
            label *= ", "
        end
        label *= f"{OASPL(pth):0.2f} dB"
    end



    psd = PowerSpectralDensityAmplitude(pth)

    pbs = ProportionalBandSpectrum(ExactProportionalBands{n}, psd)
    cbands = center_bands(pbs)
    pbs_level = 10 * log10.(pbs / P_REF^2)

    if xmax === nothing
        xmax = cbands[end]
    end
    if xmin === nothing
        xmin = cbands[1]
    end

    if aweighting
        pbs_level = dB2dBA.(cbands, pbs_level)
        pbs_level = clean_data(cbands, pbs_level)
    end

    _pbs_level = pbs_level[(cbands.>xmin).&(cbands.<xmax)]
    if ymax === nothing
        ymax = maximum(_pbs_level) + 2
    end
    if ymin === nothing
        ymin = maximum([minimum(_pbs_level) - 2, 0])
    end


    xlim = (xmin, xmax)
    ylim = (ymin, ymax)

    if title === nothing
        if aweighting
            title = "A-Weighted 1/3-Octave Band Spectrum"
        else
            title = "1/3-Octave Band Spectrum"
        end
    end

    if ylabel === nothing
        if aweighting
            ylabel = L"$L_{pA}$ [dB re 20 μPa]"
        else
            ylabel = L"$L_p$ [dB re 20 μPa]"
        end
    end

    if xlabel === nothing
        xlabel = "f [Hz]"
    end

    pplt = pyimport("proplot")
    fig, ax = pplt.subplots(figsize=figsize)

    if label === nothing
        ax[1].plot(cbands, pbs_level, lw=lw, alpha=alpha)
    else
        ax[1].plot(cbands, pbs_level, label=label, lw=lw, alpha=alpha)
    end

    ax[1].set(
        xscale="log",
        xlim=xlim,
        ylim=ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end

    return fig
end
function plot_proportional_spectrum(
    list_pth::Vector;
    n=3,
    aweighting=false,
    figsize=(6, 3),
    fname=nothing,
    xmin=62.5,
    xmax=nothing,
    ymin=nothing,
    ymax=nothing,
    title=nothing,
    ylabel=nothing,
    xlabel=nothing,
    label=nothing,
    dpi=300,
    lw=1,
    alpha=1,
)
    pplt = pyimport("proplot")
    fig, ax = pplt.subplots(figsize=figsize)

    list_xmax = []
    list_xmin = []
    list_ymax = []
    list_ymin = []

    pth1 = list_pth[1]
    psd1 = PowerSpectralDensityAmplitude(pth1)
    pbs1 = ProportionalBandSpectrum(ExactProportionalBands{n}, psd1)
    cbands = center_bands(pbs1)
    if xmax === nothing
        push!(list_xmax, cbands[end])
    else
        push!(list_xmax, xmax)
    end
    if xmin === nothing
        push!(list_xmin, cbands[1])
    else
        push!(list_xmin, xmin)
    end
    xmax = maximum(list_xmax)
    xmin = minimum(list_xmin)
    xlim = (xmin, xmax)

    for i = 1:length(list_pth)
        pth = list_pth[i]
        psd = PowerSpectralDensityAmplitude(pth)

        pbs = ProportionalBandSpectrum(ExactProportionalBands{n}, psd)
        cbands = center_bands(pbs)
        pbs_level = 10 * log10.(pbs / P_REF^2)

        if aweighting
            pbs_level = dB2dBA.(cbands, pbs_level)
            pbs_level = clean_data(cbands, pbs_level)
        end

        if label === nothing
            ax[1].plot(cbands, pbs_level, lw=lw, alpha=alpha)
        else
            ax[1].plot(cbands, pbs_level, label=label[i], lw=lw, alpha=alpha)
        end

        _pbs_level = pbs_level[(cbands.>xmin).&(cbands.<xmax)]
        if ymax === nothing
            push!(list_ymax, maximum(_pbs_level) + 2)
        else
            push!(list_ymax, ymax)
        end
        if ymin === nothing
            push!(list_ymin, maximum([minimum(_pbs_level) - 2, 0]))
        else
            push!(list_ymin, ymin)
        end
    end


    ylim = (minimum(list_ymin), maximum(list_ymax))

    if title === nothing
        if aweighting
            title = "A-Weighted 1/3-Octave Band Spectrum"
        else
            title = "1/3-Octave Band Spectrum"
        end
    end

    if ylabel === nothing
        if aweighting
            ylabel = L"$L_{pA}$ [dB re 20 μPa]"
        else
            ylabel = L"$L_p$ [dB re 20 μPa]"
        end
    end

    if xlabel === nothing
        xlabel = "f [Hz]"
    end

    ax[1].set(
        xscale="log",
        xlim=xlim,
        ylim=ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end

    return fig
end