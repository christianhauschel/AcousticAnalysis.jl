using DSP

"""
    highpass(pth::PressureTimeHistory; n_poles=10, attenuation_ripple=0.5, attenuation_stopband=60)

Highpass filter a pressure time history.

# Arguments
- `pth`: pressure time history
- `n_poles`: number of poles
- `attenuation_ripple`: attenuation in the passband [dB]
- `attenuation_stopband`: attenuation in the stopband [dB]
"""
function highpass(pth::PressureTimeHistory, f_highpass; n_poles=10, attenuation_ripple=0.5, attenuation_stopband=60.0)
    fs = 1 / timestep(pth)
    responsetype = Highpass(f_highpass; fs=fs)
    designmethod = Elliptic(n_poles, attenuation_ripple, attenuation_stopband)
    filter = digitalfilter(responsetype, designmethod)
    p_filter = DSP.filt(filter, pressure(pth))
    return PressureTimeHistory(p_filter, timestep(pth))
end