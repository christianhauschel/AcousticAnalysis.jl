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




fig2 = plot_spectrogram(
    pth; y_max=400, noise_floor_dB=-100, t_window=0.1, x_min=nothing, x_max=1.0
)

pth = time_range(pth, 0.0, 1.0)
pth = resample(pth, 3000)

fig1 = plot_msst(
    pth; length_window=256, 
    n_iterations=5, noise_floor_dB=20, 
    y_max=400, max_ram=4.0
)