using AcousticAnalysis
using WAV
using Statistics
using PyFormattedStrings

# read wav file
s, fs = wavread("data/sound.wav")

s = s[1:end,1]
fs = Float64(fs)




plot_spectrogram(s, fs)
plot_spectrum(s, fs)


using PyPlot, PyCall
using DSP, FFTW

window = DSP.hanning(length(s))

s_win = s .* window

X = fft(s_win)
X_ss = 2.0 * X
X_ss[1] = X_ss[1] / 2.0

N = length(s)

X_ss_amp = abs.(X_ss)
X_ss_amp_N = X_ss_amp / N

n_freq = length(X_ss_amp)
f_nyquist = fs / 2.0
f = LinRange(0, f_nyquist, n_freq)
Δf = f[2] - f[1]

X_amp_N_dB = 20.0 * log10.(X_ss_amp_N / P_REF)

pplt = pyimport("proplot")
fig, ax = pplt.subplots(figsize=(7,3))
ax[1].plot(f, X_amp_N_dB)
ax[1].set(
    title=f"Spectrum",
    xlabel="f [Hz]",
    ylabel=L"$|X_{ss}|/N$",
    ylim=(0, maximum(X_amp_N_dB)*1.05),
    xlim=(50, f_nyquist),
    xscale="log",
    # ultitle=r"$\Delta f = $" * f"{Δf} Hz",
)
fig