module AcousticAnalysis
using PyPlot, PyCall
using DSP
using Statistics, NaNStatistics
using LinearAlgebra
using PyFormattedStrings

export plot_spectrum, plot_spectrogram, plot_msst, P_REF

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
    fname=nothing, dpi=300, figsize=(7, 3), bpf=nothing, n_bpf=4, y_lim=(-20, 60)
)

    # Spectrum
    pplt = pyimport("proplot")
    fig, ax = plt.subplots(figsize=figsize, tight_layout=true)
    spec, freqs, _ = ax.magnitude_spectrum(
        p / P_REF,
        Fs=fs,
        scale="dB",
        sides="onesided",
        label="loading",
        # window=DSP.hanning(n_τ),
    )

    if bpf !== nothing
        y_min = y_lim[1]
        y_max = y_lim[2]
        for i = 1:n_bpf
            ax.vlines(bpf * i, y_min, y_max, color="k", lw=0.5)
        end
    else
        bpf = 50
    end

    ax.set(
        xscale="log",
        title="Spectrum",
        ylim=y_lim,
        xlim=(minimum([100, bpf * 0.8]), fs / 2.0),
    )
    ax.legend()

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)

    end
end

end
