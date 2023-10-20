function plot_spectrogram(
    p, fs;
    fname=nothing,
    dpi=300,
    figsize=(6, 5),
    t_window=0.05,
    onesided=true,
    title="Spectrogram",
    ylim=nothing
)

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
end

function plot_spectrum(
    p::Vector, fs::Number;
    fname=nothing, 
    dpi=300, 
    figsize=(6, 3),
    bpf=nothing,
    n_bpf=10,
    ymin=0, 
    ymax=nothing,
    title="Spectrum", 
    xlabel="f [Hz]", 
    ylabel=nothing,
    type=:amplitude,
    window=:hanning,
    label=nothing,
)   
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
        X_psd = abs2.(X_ss) / (fs * sum(window.^2))
        X_dB = 10.0 * log10.(X_psd * Δf / P_REF^2)
    end

    pplt = pyimport("proplot")
   
    fig, ax = pplt.subplots(figsize=figsize)

    if label === nothing
        ax[1].plot(f, X_dB)
    else
        ax[1].plot(f, X_dB, label=label)
    end

    if ymax === nothing
        ymax = maximum(X_dB) + 0.1 * mean(X_dB)
    end
    if ymin === nothing
        ymin = minimum(X_dB) - 0.1 * mean(X_dB)
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
            ylabel = L"$L_p$" * " [dB re 20 μPa]"
        else
            ylabel = "PSD [dB re 20 μPa]"
        end
    end

    ax[1].set(
        title=title,
        xscale="log",
        xlabel=xlabel,
        ylabel=ylabel,
        ylim=(ymin, ymax),
        xlim=(minimum([100, bpf * 0.8]), fs / 2.0),
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)

    end
end
function plot_spectrum(
    list_p::Vector{Vector{Float64}}, 
    list_fs::Vector;
    fname=nothing, 
    dpi=300, 
    figsize=(6, 3),
    bpf=[],
    n_bpf=10,
    ymin=0, 
    ymax=nothing,
    title="Spectrum", 
    xlabel="f [Hz]", 
    ylabel=nothing,
    type=:amplitude,
    window=:hanning,
    label=nothing,
)   
    pplt = pyimport("proplot")
    
    fig, ax = pplt.subplots(figsize=figsize)

    list_ymax = []
    list_ymin = []

    for k in 1:length(list_p)
        p = list_p[k]
        fs = list_fs[k]

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
            X_psd = abs2.(X_ss) / (fs * sum(_window.^2))
            X_dB = 10.0 * log10.(X_psd * Δf / P_REF^2)
        end

        if label === nothing
            ax[1].plot(f, X_dB)
        else
            ax[1].plot(f, X_dB, label=label[k])
        end


        if ymax === nothing
            push!(list_ymax, maximum(X_dB) + 0.1 * mean(X_dB))
        end
        if ymin === nothing
            push!(list_ymin, minimum(X_dB) - 0.1 * mean(X_dB))
        end
    end

    if length(bpf) > 0
        plot_bpf = true
    else
        plot_bpf = false
        push(bpf, 62.5)
    end

    if ymax === nothing
        ymax = maximum(list_ymax)
    end
    if ymin === nothing
        ymin = minimum(list_ymin)
    end

    if plot_bpf
        for k in 1:length(list_fs)
            for i = 1:n_bpf
                ax[1].vlines(bpf[k] * i, ymin, ymax, color=f"C{k-1}", lw=0.5)
            end
        end
    end


    if ylabel === nothing
        if type == :amplitude
            ylabel = L"$L_p$" * " [dB re 20 μPa]"
        else
            ylabel = "PSD [dB re 20 μPa]"
        end
    end

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    ax[1].set(
        title=title,
        xscale="log",
        xlabel=xlabel,
        ylabel=ylabel,
        ylim=(ymin, ymax),
        xlim=(minimum([100, minimum(bpf) * 0.8]), maximum(list_fs) / 2.0),
    )

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)

    end
end

function plot_history(
    list_p::Vector{Vector{Float64}},
    list_fs::Vector;
    list_t=nothing,
    title="Time History",
    fname=nothing,
    dpi=300,
    figsize=(6, 2),
    xlabel="t [s]",
    ylabel="p [Pa]"
)
    pplt = pyimport("proplot")

    fig, ax = pplt.subplots(figsize=figsize)

    for i = 1:length(list_p)
        p = list_p[i]
        N = length(p)

        Δt = 1.0 / list_fs[i]
        if list_t === nothing
            t = Vector(Δt * (0:N-1))
        else
            t = list_t[i]
        end

        ax[1].plot(t, p)
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
end
function plot_history(
    p::Vector,
    fs;
    t=nothing,
    title="Time History",
    fname=nothing,
    dpi=300,
    figsize=(6, 2),
    xlabel="t [s]",
    ylabel="p [Pa]"
)
    pplt = pyimport("proplot")

    fig, ax = pplt.subplots(figsize=figsize)

    N = length(p)

    Δt = 1.0 / fs
    if t === nothing
        t = Vector(Δt * (0:N-1))
    end

    ax[1].plot(t, p)
    

    ax[1].set(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )
    fig

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end
end