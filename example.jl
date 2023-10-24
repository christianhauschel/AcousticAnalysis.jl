using Revise
using AcousticAnalysis
using WAV
using Statistics
using PyFormattedStrings
using PyPlot, PyCall
pplt = pyimport("proplot")
plt = pyimport("matplotlib.pyplot")
pplt.close("all")

# read wav file
s1, fs1 = wavread("data/01.wav")
s2, fs2 = wavread("data/00.wav")

s1 = s1[400:end,1]
s2 = s2[400:end,1]
fs1 = Float64(fs1)
fs2 = Float64(fs2)
p1 = s1 .- mean(s1)
p2 = s2 .- mean(s2)

pth1 = PressureTimeHistory(p1, 1/fs1)
pth2 = PressureTimeHistory(p2, 1/fs2)

bpf1=66.6667
bpf2=100

# plot_spectrogram(pth1; ylim=(0, 1000), fname="out/spectrogram.png", t_window=0.5)
# plot_history(pth2; fname="out/time_history.png", label="2")
# plot_narrowband_spectrum(pth1; bpf=bpf1, fname="out/spectrum_psd.png", type=:psd, aweighting=true)
# plot_narrowband_spectrum(p2, fs2; bpf=bpf2, fname="out/spectrum_psa.png")
# plot_narrowband_spectrum([pth1, pth2]; bpf=[bpf1,bpf2], fname="out/spectrum_narrowband.png", label=["00", "01"])
plot_history([pth1, pth2]; label=["1", "2"], fname="out/time_history.png")
# plot_proportional_spectrum(pth2; aweighting=false)
plot_proportional_spectrum([pth1, pth2]; aweighting=true, fname="out/spectrum_proportional_A.png", label=["1", "2"])

# OASPL(pth1)


pth1_A = aweighting(pth1)

# plot_narrowband_spectrum(pth1)
# plot_proportional_spectrum(pth1, fname="out/spectrum_proportional.png")
plot_narrowband_spectrum([pth1, pth1_A], label=[f"{OASPL(pth1):0.2f} dB", f"{OASPL(pth1_A):0.2f} dB(A)"], fname="out/spectrum_narrow.png")


L_p = OASPL(pth1)
L_W = pressure2power(L_p; r=1, Q=1)

