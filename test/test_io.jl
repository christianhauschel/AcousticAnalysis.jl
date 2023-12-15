using AcousticAnalysis

# Test 
fs = 2000
pth1 = PressureTimeHistory(5*(rand(20001).-0.7), 1/fs, 0.0)
pth2 = PressureTimeHistory(10*(rand(20001).-0.7), 1/fs, 0.0)

pths = [pth1, pth2]


plot_history(pths)
pths_norm = remove_dc_offset.(pths)
pths_norm = normalize(pths_norm; offset_peak_dB=-1.0)
plot_history(pths_norm)

save_wav.(pths_norm, ["out/01.wav", "out/02.wav"])
save_h5(pths_norm, "out/test.h5", ["01", "02"])