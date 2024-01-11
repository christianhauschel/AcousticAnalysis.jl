function mel_scale(f::Float64)
    return 2595 * log10(1 + f / 700)
end

function log_scale(f::Float64)
    return log10(f)
end

function linear_scale(f::Float64)
    return f
end

"""
Bark scale.

# References
https://en.wikipedia.org/wiki/Bark_scale
"""
function bark_scale(f::Float64)
    return 13 * atan(0.00076 * f) + 3.5 * atan((f / 7500)^2)
end


"""
Equivalent Rectangular Bandwidth (ERB) scale.

# References
https://en.wikipedia.org/wiki/Equivalent_rectangular_bandwidth#Approximations
"""
function erb_scale(f::Float64)
    return 6.23 * f^2 + 93.39 * f + 28.52
end


