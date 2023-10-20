module AcousticAnalysis
using PyPlot, PyCall
using DSP
using Statistics, NaNStatistics
using LinearAlgebra
using PyFormattedStrings
using FFTW

export P_REF

const P_REF = 20e-6

include("plots.jl")
export plot_spectrum, plot_spectrogram, plot_history


end
