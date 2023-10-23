#  %%

import numpy as np
from scipy.signal import zpk2tf, freqs


z = [0, 0, 0, 0]
p = [
    -2*np.pi*20.598997057568145,
    -2*np.pi*20.598997057568145,
    -2*np.pi*12194.21714799801,
    -2*np.pi*12194.21714799801,
    -2*np.pi*107.65264864304628,
    -2*np.pi*737.8622307362899,
]
k = 7390100803.660323

b, a = zpk2tf(z, p, k)
# k /= abs(freqs(b, a, [2*np.pi*1000])[1][0])

f = np.logspace(0, 5, 1000)
omega = 2*np.pi*f
w, h = freqs(b, a, omega)

import proplot as pplt

fig, ax = pplt.subplots(figsize=(6,4))
ax.plot(w, 20*np.log10(abs(h)))
ax.format(
    # xlocator="log",
    xscale="log",
    ylim=(-50, 20),
    xlim=(10, 100000),
    # ylocator="log",
    xlabel="Frequency [Hz]",
    ylabel="Magnitude [dB]",
    title="Frequency Response",
)




# %%