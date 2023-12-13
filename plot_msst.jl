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

# read wav file
s_old, fs_old = wavread("data/matrice_hover_exp.wav")
s_old, fs_old = wavread("data/matrice_hover_sim.wav")
s_old = s_old[:,1]
fs = 2000.0

s = DSP.Filters.resample(s_old, fs / fs_old)

max_ram = 3.0

pth0 = PressureTimeHistory(s[:, 1], 1 / fs)
t0 = Vector(AcousticAnalysis.time(pth0))

# id_min = 20000
# id_max = 25000

# p = s[id_min:id_max, 1]
# t  = t0[id_min:id_max]
p = s 
t  = t0

n_blades = 2
rpm = 2800
bpf = rpm / 60 * n_blades



function predict_ram(n)
    return (0.0001n^2 + 0.0164n + 7.6607) / 1000
end

function predict_n(ram)
    f = n -> predict_ram(n) - ram
    return find_zero(f, 1000)
end


if predict_ram(length(p)) > max_ram
    n_max = round(Int, predict_n(max_ram))
    p = p[1:n_max]

    println(f"Warning: RAM usage is predicted to be greater than {max_ram} GB. Restricting signal size to {n_max} samples.")
end


pth = PressureTimeHistory(p, 1 / fs, t[1])
t = Vector(AcousticAnalysis.time(pth))


n_iterations = 4
length_window = 256 # higher -> better resolution in frequency, lower in time

@time MSST = msst(pth, n_iterations; length_window=length_window);




y_min = nothing
y_max = bpf * 4
x_min = nothing
x_max = nothing

noise_floor = 20
T_msst_plot = 20.0 * log10.(abs.(MSST) / P_REF)
T_msst_plot_min = nanminimum(T_msst_plot)
T_msst_plot = ifelse.(isinf.(T_msst_plot), noise_floor, T_msst_plot)


f_msst = Vector(LinRange(0, fs / 2, size(MSST, 1)))

fig, ax = plt.subplots(figsize=(7, 3))

if y_min === nothing
    y_min = f_msst[1]
end
if y_max === nothing
    y_max = f_msst[end]
end

mesh = ax.pcolorfast(
    t,
    f_msst,
    T_msst_plot,
    vmin=maximum([T_msst_plot_min, noise_floor]),
    vmax=maximum(T_msst_plot),
    rasterized=true,
)

cbar = fig.colorbar(mesh, ax=ax, label="dB")
ax.set(
    ylim=(y_min, y_max),
    title="MSST",
    xlabel="t [s]",
    ylabel="f [Hz]",
)

fig