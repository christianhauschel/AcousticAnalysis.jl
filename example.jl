using AcousticAnalysis
using WAV
using Statistics
using PyFormattedStrings

# read wav file
s, fs = wavread("data/sound.wav")

s = s[1:end,1]
fs = Float64(fs)


plot_spectrogram(s, fs; fname="plots/spectrogram.png")
plot_spectrum(s, fs; fname="plots/spectrum.png", ymax=100)
