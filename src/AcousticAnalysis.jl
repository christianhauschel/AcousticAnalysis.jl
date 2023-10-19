module AcousticAnalysis
using PyPlot, PyCall
using DSP
using Statistics, NaNStatistics
using LinearAlgebra
using PyFormattedStrings
using FFTW

export plot_spectrum, plot_spectrogram, P_REF

const P_REF = 20e-6

function plot_spectrogram(
    p, fs;
    fname=nothing, 
    dpi=300, 
    figsize=(6, 5),
    t_window=0.05, 
    onesided=true,
    title="Spectrogram",
    ylim=nothing,
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
    p, fs;
    fname=nothing, dpi=300, figsize=(6, 3), 
    bpf=nothing, 
    n_bpf=4, 
    ymin=0, ymax=nothing, 
    title="Spectrum", xlabel="f [Hz]", ylabel="Amplitude [dB]",
)   
    N = length(p)
    window = DSP.hanning(N)
    p_win = p .* window

    X = fft(p_win)
    X_ss = 2.0 * X
    X_ss[1] = X_ss[1] / 2.0
    X_ss_amp = abs.(X_ss)
    X_ss_amp_N = X_ss_amp / sum(window)

    n_freq = length(X_ss_amp)
    f_nyquist = fs / 2.0
    f = LinRange(0, f_nyquist, n_freq)
    Δf = f[2] - f[1]

    X_amp_N_dB = 20.0 * log10.(X_ss_amp_N / P_REF)


    pplt = pyimport("proplot")
    fig, ax = pplt.subplots(figsize=figsize)

    ax[1].plot(f, X_amp_N_dB)
    # spec, freqs, _ = ax.magnitude_spectrum(
    #     p / P_REF,
    #     Fs=fs,
    #     scale="dB",
    #     sides="onesided",
    #     label="loading",
    #     # window=DSP.hanning(n_τ),
    # )

    if ymax === nothing
        ymax = maximum(X_amp_N_dB)
    end
    if ymin === nothing
        ymin = minimum(X_amp_N_dB)
    end

    if bpf !== nothing
        for i = 1:n_bpf
            ax.vlines(bpf * i, y_min, y_max, color="k", lw=0.5)
        end
    else
        bpf = 62.5
    end

    ax[1].set(
        title=title,
        xscale="log",
        xlabel=xlabel,
        ylabel=ylabel,
        ylim=(ymin, ymax),
        xlim=(minimum([100, bpf * 0.8]), fs / 2.0),
    )

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)

    end
end

end
