using AcousticAnalysis
using WAV
using Statistics
using PyFormattedStrings
using NaNStatistics
using Roots
using PyPlot, PyCall
import DSP
pplt = pyimport("proplot")
plt = pyimport("matplotlib.pyplot")
pplt.close("all")


p, fs = wavread("data/matrice_hover_exp.wav")
p = p[:, 1]

pth = PressureTimeHistory(p, 1 / fs, 0.0)
pth = time_range(pth, 0.0, 1.0)

plot_msst(pth; fs_resample=3000, length_window=256, n_iterations=4, noise_floor_dB=20, y_max=400, max_ram=4.0, fname="msst.png")