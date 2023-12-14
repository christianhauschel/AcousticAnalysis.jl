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


n = 2000
t = Vector(LinRange(0, 0.5, n))
dt = t[2] - t[1]
A1 = t -> 0.5 #* t + 0.5
A2 = t -> 1 # * t + 0.5
f1 = t -> -400 * t .^ 3 + 200
f2 = t -> 100 * t .^ 2 + 100
p = A1.(t) .* sin.(2 * pi .* f1.(t) .* t) .+ A2.(t) .* sin.(2 * pi .* f2.(t) .* t)

pth = PressureTimeHistory(p, dt, 0.0)

fname="doc/img/hist_msst.png"
fig3 = plot_history(time_range(pth, 0, 1); 
    lw=0.5, fname=fname, dpi=300,
    figsize=(6, 1)
)

fname="doc/img/spectro_msst.png"
fig2 = plot_spectrogram(
    pth; y_max=500, t_window=0.05, noise_floor_dB=60.0, 
    fname=fname, x_min=nothing, x_max=nothing,
    figsize=(6, 3.7)
)

fname="doc/img/msst.png"
fig1 = plot_msst(
    pth; length_window=300,
    n_iterations=8, noise_floor_dB=60,
    y_max=500, max_ram=4.0, fname=fname,
    figsize=(6, 3.7)
)