module AcousticAnalysis
using PyPlot, PyCall
using DSP
using Statistics, NaNStatistics
using LinearAlgebra
using PyFormattedStrings
using FFTW
using AcousticMetrics
using AcousticMetrics: starttime, timestep, time, pressure, frequency, halfcomplex, W_A, AbstractNarrowbandSpectrum, OASPL
using AcousticMetrics: PowerSpectralDensityAmplitude, MSPSpectrumAmplitude, ProportionalBandSpectrum, ExactProportionalBands, center_bands, lower_bands, upper_bands, AbstractPressureTimeHistory
using FLOWMath: linear

P_REF = 20e-6
export P_REF

export PressureTimeHistory, pressure, time, timestep, starttime, frequency, OASPL

include("signal.jl")
export resample, time_range, _resample, _time_range

include("spectrum.jl")
export narrowband_spectrum, propband_spectrum, _narrowband_spectrum, _propband_spectrum

include("msst.jl")
export msst, _msst

include("filter.jl")
export highpass, _highpass

include("utils.jl")
export dbsum, dB, pressure2power, remove_dc_offset, normalize, _normalize, _remove_dc_offset

include("weighting.jl")
export dB2dBA, aweighting, _aweighting

include("plots.jl")
export plot_narrowband_spectrum, plot_spectrogram, plot_history, plot_propband_spectrum, plot_msst

include("io.jl")
export save_h5, load_h5, load_wav, save_wav, _save_wav

end
