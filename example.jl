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
calibration_factor = 1 #256 * 1e-9
s1, fs1 = wavread("data/01.wav") #, format="native")
s2, fs2 = wavread("data/00.wav") #, format="native")

s1 = s1[300:end,1] * calibration_factor
s2 = s2[300:end,1] * calibration_factor
fs1 = Float64(fs1)
fs2 = Float64(fs2)


pth1 = PressureTimeHistory(s1, 1/fs1)
pth2 = PressureTimeHistory(s2, 1/fs2)
# pth1 = remove_dc_offset(pth1)
# pth2 = remove_dc_offset(pth2)

bpf1=66.6667
bpf2=100

# plot_spectrogram(pth2; ylim=(0, 1000), fname="out/spectrogram.png", t_window=0.5)
# plot_history(pth1; fname="out/time_history.png", label="2", alpha=0.5, lw=2, remove_dc_offset=true)
# fig1 = plot_narrowband_spectrum(pth1; bpf=bpf1, fname="out/spectrum_psd.png", type=:amplitude, aweighting=false)
# fig = plot_narrowband_spectrum(pth1; bpf=bpf1, fname="out/spectrum_psd.png", type=:amplitude, aweighting=false)
# fig2 = plot_narrowband_spectrum(p2, fs2; bpf=bpf2, fname="out/spectrum_psa.png")
# fig = plot_narrowband_spectrum([pth1, pth2]; bpf=[bpf1,bpf2], fname="out/spectrum_narrowband.png", label=["00", "01"], aweighting=true)
# plot_history([pth1, pth2]; label=["1", "2"], fname="out/time_history.png", without_dc_offset=true)
# fig = plot_propband_spectrum(pth2; aweighting=false)
fig = plot_propband_spectrum([pth1, pth2]; aweighting=true, fname="out/spectrum_proportional_A.png", label=["1", "2"], alpha=1, lw=1)

# OASPL(pth1)


# pth1_A = aweighting(pth1)

# plot_narrowband_spectrum(pth1)
# plot_propband_spectrum(highpass.([pth1, pth2], 1; attenuation_stopband=100))
# plot_narrowband_spectrum(pth1)
# plot_narrowband_spectrum([pth1, pth1_A], label=[f"{OASPL(pth1):0.2f} dB", f"{OASPL(pth1_A):0.2f} dB(A)"], fname="out/spectrum_narrow.png")

# algorithm to normalize pth1 to -1...1

pth1 = remove_dc_offset(pth1)
ps_norm = normalize([pth1, pth2]; offset_peak_dB=-1.0)

# wavwrite(pth1_norm.p, "out/01.wav", Fs=fs1)


# f, X_dB, X_ss_final = narrowband_spectrum(pth1; type=:psd, aweighting=true)

# L_p = OASPL(pth1)
# L_W = pressure2power(L_p; r=1, Q=1)

