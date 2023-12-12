module AcousticAnalysis
using PyPlot, PyCall
using DSP
using Statistics, NaNStatistics
using LinearAlgebra
using PyFormattedStrings
using FFTW
using AcousticMetrics
using AcousticMetrics: starttime, timestep, time, pressure, frequency, halfcomplex, W_A, AbstractNarrowbandSpectrum, OASPL
using AcousticMetrics: PowerSpectralDensityAmplitude, MSPSpectrumAmplitude, ProportionalBandSpectrum, ExactProportionalBands, center_bands, lower_bands, upper_bands
using FLOWMath: linear

export PressureTimeHistory, pressure, time, timestep, starttime, frequency, OASPL

P_REF = 20e-6
export P_REF

include("spectrum.jl")
export narrowband_spectrum

include("filter.jl")
export highpass

include("utils.jl")
export dbsum, dB, pressure2power, remove_dc_offset, normalize

include("weighting.jl")
export dB2dBA, aweighting

include("plots.jl")
export plot_narrowband_spectrum, plot_spectrogram, plot_history, plot_proportional_spectrum

include("msst.jl")
export msst

end
