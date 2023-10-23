module AcousticAnalysis
using PyPlot, PyCall
using DSP
using Statistics, NaNStatistics
using LinearAlgebra
using PyFormattedStrings
using FFTW
using AcousticMetrics
using AcousticMetrics: starttime, timestep, time, pressure, frequency, halfcomplex, p_ref, W_A, AbstractNarrowbandSpectrum, OASPL
using AcousticMetrics: PowerSpectralDensityAmplitude, MSPSpectrumAmplitude, ProportionalBandSpectrum, ExactProportionalBands, center_bands, lower_bands, upper_bands
using FLOWMath: linear
export PressureTimeHistory

include("utils.jl")
export dbsum, dB

include("weighting.jl")
export dB2dBA, aweighting


export pressure, time, timestep, starttime, frequency, OASPL


include("plots.jl")
export plot_narrowband_spectrum, plot_spectrogram, plot_history, plot_proportional_spectrum


end
