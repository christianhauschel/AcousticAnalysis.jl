using LinearAlgebra
using FFTW
using Statistics
using Base.Threads
using Roots
using DSP

n = 10
t = Vector(LinRange(0, 1, n))
x = sin.(2π * 10 * t)
length_window=2
n_iterations=2

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
omega2 = zeros(Float64, n_half, tcol)
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
omega = round.(Int, omega)
omega = Int.(omega)

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

omega2 = Int.(omega2)

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


TFT = tfr / (N / 2)
MSST = Ts / (N / 2)
