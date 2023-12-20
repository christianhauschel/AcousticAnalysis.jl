
using AcousticAnalysis

fs = 10000.0
t = 0.0:1/fs:prevfloat(1.0)
f = 100
y = sin.(2pi * f * t) * 0.1

pth = PressureTimeHistory(y, 1 / fs, 0.0)

save_wav(pth, "out/test.wav"; nbits=64)

pth_load = load_wav("out/test.wav"; format="double")
maximum(pth_load.p)

plot_history(pth)

save_wav(pth, "out/test.wav")

using WAV
wavwrite(y, "out/test.wav", Fs=fs, nbits=64)


pth_new, nbits, opt = load_wav("data/01.wav"; calibration_factor=1, format="native")


y, fs, nbits, opt = wavread("data/hover_final.wav", format="native")
pth1 = PressureTimeHistory(y[:, 1] * 1e-6 / 20, 1 / fs, 0.0)

y, fs, nbits, opt = wavread("data/hover_raw.wav", format="native")
pth2 = PressureTimeHistory(y[:, 1] * 1e-9, 1 / fs, 0.0)

plot_history(pth2)

plot_history([pth1, pth2])

# pth_load = load_wav("out/test.wav"; calibration_factor=1, format="native")
# plot_history(pth_load)