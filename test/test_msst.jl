using AcousticAnalysis
using DelimitedFiles
using Statistics
using CSV, DataFrames


function readfile(fname; delimiter=',', dtype=ComplexF64)
    lines = readlines(fname)
    _data = split(lines[1], delimiter)
    m = length(_data)
    n = length(lines)
    data = zeros(dtype, n, m)
    for i in 1:n
        _data = split(lines[i], delimiter)
        for j = 1:m
            data[i,j] = parse(dtype, _data[j])
        end
    end
    return data
end

msst2_ref = readfile("../kunyu_msst2/msst2.csv")


# n = size(msst2_ref, 2)
# m = size(msst2_ref, 1)
n = 10
t = Vector(LinRange(0, 1, n))
p = sin.(2π * 10 * t)

msst2 = AcousticAnalysis._msst_core(p, 2; length_window=2)


# rmse = sqrt(
#     sum((abs.(msst2_ref) .- abs.(msst2)).^2) / (m*n) / sum(abs.(msst2_ref) .^ 2)
# )

msst2_ref[1,1]
msst2[1,1]