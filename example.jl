using Revise
using AcousticAnalysis
using WAV
using Statistics
using PyFormattedStrings

using PyPlot, PyCall
pplt = pyimport("proplot")
pplt.close("all")

# read wav file
s1, fs1 = wavread("data/00.wav")
s2, fs2 = wavread("data/01.wav")

s1 = s1[400:end,1]
s2 = s2[400:end,1]
fs1 = Float64(fs1)
fs2 = Float64(fs2)
p1 = s1 .- mean(s1)
p2 = s2 .- mean(s2)

bpf1=100
bpf2=66.6667

# plot_spectrogram(p1, fs1; ylim=(0, 1000), fname="plots/spectrogram.png", t_window=0.5)
# plot_spectrum(p1, fs1; bpf=bpf1, fname="plots/spectrum.png")

plot_spectrum([p1,p2], [fs1,fs2]; bpf=[bpf1,bpf2], fname="plots/spectrum.png", label=["00", "01"])