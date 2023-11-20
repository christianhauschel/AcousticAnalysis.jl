using AcousticAnalysis

fname = "pressure"
dir = "/home/hacs/projects/neaptide/flowunsteady/rotor_simple/data/hover_mid_FIXED/acoustics/runcase"

import FLOWNoise as FN
header, field = FN.read_wopwopoutput(fname; read_path=dir, tec=false)

n_grid = size(field, 1)
n_time = size(field, 3)
n_labels = size(field, 4)

pths_FU = Vector{PressureTimeHistory}(undef, n_grid)
for i in 1:n_grid 
    # calculate fs 
    t = field[i,1,:,1]
    fs = 1/(t[2]-t[1])

    p_total = field[i,1,:,4]

    pths_FU[i] = PressureTimeHistory(p_total, 1/fs)
end

# plot_spectrogram(pths_total[1]; ylim=(0, 1000), t_window=0.1)

plot_narrowband_spectrum(pths_FU; bpf=100 * ones(n_grid), lw=0.5)

# plot_history(pths_total; lw=0.5)