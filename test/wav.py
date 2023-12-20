# %%

# read wav file using scipy 

import scipy.io.wavfile as wavfile
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt

# %%

fs = 10000
t = np.arange(0, 1, 1/fs)
f = 100
y = 0.5*np.sin(2 * np.pi * f * t)


# wave wav 
wavfile.write('out/test_python.wav', fs, y)

fs1, y1 = wavfile.read('out/test_python.wav')


# %%

import h5py
with h5py.File('out/test_h5.h5', 'w') as f:
    f.create_dataset('p', data=y)
    f.create_dataset('fs', data=fs)
    f.create_dataset('t0', data=0.0)
    g = f.create_group("correction")
    g.create_dataset('factor_microphone', data=1.0)
    g = f.create_group('units')
    g.create_dataset('p', data='Pa', dtype='S2')
    g.create_dataset('fs', data='Hz', dtype='S2')
    g.create_dataset('t0', data='s', dtype='S1')


# %%