using LinearAlgebra
using FFTW
using Statistics
using Base.Threads
using Roots

"""
    msst(pth::PressureTimeHistory, length_window::Int, n_iterations::Int)

Computes the multisynchrosqueezing transform (MSST) of the signal using the 
Expression (31)-based Algorithm [1].

# Arguments
- `pth::PressureTimeHistory`: Pressure time history
- `length_window::Int`: Length of the window function
- `n_iterations::Int`: Number of iterations

# Retuns
- `MSST::Array{ComplexF64,2}`: Multisynchrosqueezing transform

# References
[1] G. Yu, Z. Wang, and P. Zhao, “Multisynchrosqueezing Transform,” IEEE
    Transactions on Industrial Electronics, vol. 66, pp. 5441–5455,
    Jul. 2019, doi: 10.1109/TIE.2018.2868296.
    
[2] G. Yu, “A multisynchrosqueezing-based high-resolution time-frequency
    analysis tool for the analysis of non-stationary signals,” Journal of
    Sound and Vibration, vol. 492, p. 115813, Oct. 2020,
    doi: 10.1016/j.jsv.2020.115813.
"""
function _msst(pth::PressureTimeHistory, n_iterations::Int; length_window::Int=512)
    return _msst(pth.p, n_iterations; length_window=length_window)
end
function _msst(x::Vector, n_iterations::Int; length_window::Int=512)

    N = length(x)
    hlength = length_window + 1 - length_window % 2
    ht = Vector(LinRange(-0.5, 0.5, hlength))

    # Gaussian window
    h = exp.(-π / 0.32^2 .* ht .^ 2)
    hrow = length(h)
    Lh = div(hrow - 1, 2)

    t = 1:N
    tcol = length(t)

    n_half = round(Int, N / 2)

    tfr = zeros(ComplexF64, N, tcol)
    omega = zeros(Float64, n_half, tcol-1)
    omega2 = zeros(Int64, n_half, tcol)
    Ts = zeros(ComplexF64, n_half, tcol)

    # ---------------------------------
    # Block 0 (0.3 s)
    # ---------------------------------
    # t0 = time()
    for icol in 1:tcol
        ti = t[icol]
        tau = -min(
            n_half - 1,
            Lh,
            ti - 1
        ):min(
            n_half - 1,
            Lh,
            N - ti
        )
        indices = mod.(N .+ tau, N) .+ 1

        rSig = x[ti.+tau]

        tfr[indices, icol] = rSig .* conj.(h[Lh.+tau.+1])
    end


    tfr = fft(tfr, 1) # columns as vectors for FFT
    tfr = tfr[1:n_half, :]



    # ---------------------------------
    # Block 1 (1.2 s)
    # ---------------------------------
    # t0 = time()
    for i in 1:n_half
        omega[i, :] = diff(
            unwrap(
                angle.(tfr[i, :])
            )
        ) .* N / 2 / π
    end

    omega = hcat(omega, omega[:, end])
    omega = round.(Int64, omega)
    omega = Int64.(omega)

    neta, nb = size(tfr)

    # ---------------------------------
    # Block 2 (27 s)
    # ---------------------------------
    # t0 = time()

    

    if n_iterations > 1
        for kk = 1:n_iterations-1
            for b = 1:nb
                for eta = 1:neta
                    k = omega[eta, b]
                    if 1 <= k <= neta
                        omega2[eta, b] = omega[k, b]
                    end
                end
            end
            omega = copy(omega2)
        end
    else
        omega2 = copy(omega)
    end


    # ---------------------------------
    # Block 3 (12 s)
    # ---------------------------------
    # t0 = time()


    for b = 1:nb
        for eta = 1:neta
            if abs(tfr[eta, b]) > 0.0001
                k = omega2[eta, b]
                if 1 <= k <= neta
                    Ts[k, b] += tfr[eta, b]
                end
            end
        end
    end


    # STFT = tfr / (N / 2) # Short-Time Fourier transform
    MSST = Ts / (N / 2) # Multisynchrosqueezing transform

    return MSST
end


function msst(p::Vector, t::Vector; fs_resample=nothing, n_iterations=4, length_window=256, max_ram=3.0)
    fs = 1 / (t[2] - t[1])

    if fs_resample === nothing
        fs_resample = fs
    end
    p_resampled, t_resampled = resample(p, fs, fs_resample; t0=t[1])

    """
    Empirical function to predict RAM usage given n.
    """
    function predict_ram(n)
        return (0.0001n^2 + 0.0164n + 7.6607) / 1000
    end

    """
    Empirical function to predict n given the RAM.
    """
    function predict_n(ram)
        f = n -> predict_ram(n) - ram
        return find_zero(f, 1000)
    end


    if predict_ram(length(p_resampled)) > max_ram
        n_max = round(Int, predict_n(max_ram))
        p_resampled = p_resampled[1:n_max]

        println(f"Warning: RAM usage is predicted to be greater than {max_ram} GB. Restricting signal size to {n_max} samples.")
    end


    pth = PressureTimeHistory(p_resampled, 1 / fs_resample, t[1])

    MSST = _msst(pth, n_iterations; length_window=length_window)

    return MSST, t_resampled, p_resampled
end
function msst(pth::PressureTimeHistory; fs_resample=nothing, n_iterations=4, length_window=256, max_ram=3.0)
    return msst(pth.p, Vector(time(pth)); fs_resample=fs_resample, n_iterations=n_iterations, length_window=length_window, max_ram=max_ram)
end