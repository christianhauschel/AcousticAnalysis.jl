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
p1 = s1[400:1000,1]
pth1 = PressureTimeHistory(p1, 1/fs1)
bpf1=66.6667

MSST = msst(pth1, 4; length_window=512)

# =====================
# Plot
# =====================

noise_floor = -50
T_msst_plot = 20.0 * log10.(abs.(MSST) / P_REF)
T_msst_plot_min = minimum(T_msst_plot[.!isinf.(T_msst_plot)])
T_msst_plot = ifelse.(isinf.(T_msst_plot), noise_floor, T_msst_plot)

t_msst = AcousticAnalysis.time(pth1)
f_msst = Vector(LinRange(0, fs1/2, size(MSST, 1)))

fig, ax = plt.subplots(figsize=(7, 3))

mesh = ax.pcolorfast(
    t_msst,
    f_msst,
    T_msst_plot,
    # shading="gouraud",
    vmin=maximum([T_msst_plot_min, noise_floor, maximum(T_msst_plot) - 45.0]),
    vmax=maximum(T_msst_plot),
    rasterized=true,
)  

cbar = fig.colorbar(mesh, ax=ax, label="dB")
ax.set(
    ylim=(f_msst[1], f_msst[end]),
    title="MSST",
    xlabel="t [s]",
    ylabel="f [Hz]",
)

fig