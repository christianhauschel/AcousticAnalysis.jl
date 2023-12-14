
"""
    clean_data(x, y)

Clean data by removing NaN and Inf values and interpolating linearly.

# Arguments
- `x`: x data
- `y`: y data

# Returns
- `y`: cleaned y data
"""
function clean_data(x, y)
    mask = .!isinf.(y) .& .!isnan.(y)
    _y = y[mask]
    _x = x[mask]

    return linear(_x, _y, x)
end

function plot_spectrogram(
    pth::PressureTimeHistory;
    fname=nothing,
    dpi=300,
    figsize=(6, 5),
    t_window=0.05,
    onesided=true,
    noise_floor_dB=-50.0,
    title="Spectrogram",
    x_min = nothing,
    x_max = nothing,
    y_min = 0.0,
    y_max=nothing,
)
    p = pressure(pth)
    fs = 1 / timestep(pth)

    pplt = pyimport("proplot")
    spectro = spectrogram(p / P_REF, Int(round(fs * t_window)); fs=fs, onesided=onesided)
    power_dB = 10 * log10.(spectro.power)

    fig, ax = pplt.subplots(figsize=figsize)
    cm = ax[1].pcolormesh(
        spectro.time, 
        spectro.freq, 
        power_dB,
        vmin=nanmaximum([nanminimum(power_dB), noise_floor_dB]),
        vmax=nanmaximum(power_dB),
    )
    cb = plt.colorbar(cm, label="Power [dB]")
    
    if x_min === nothing
        x_min = spectro.time[1]
    end
    if x_max === nothing
        x_max = spectro.time[end]
    end

    if y_max === nothing
        y_max = fs/2
    end
    if y_min === nothing
        y_min = 0.0
    end
    ax[1].set(
        ylim = (y_min, y_max),
        xlim = (x_min, x_max),
    )
    ax[1].set(
        xlabel="t [s]",
        ylabel="f [Hz]",
        title=title,
    )
    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end
    return fig
end

function plot_narrowband_spectrum(
    pth::PressureTimeHistory;
    fname=nothing,
    dpi=300,
    figsize=(6, 3),
    bpf=nothing,
    n_bpf=10,
    xmin=nothing,
    xmax=nothing,
    ymin=0,
    ymax=nothing,
    title="Spectrum",
    xlabel="f [Hz]",
    ylabel=nothing,
    type=:amplitude,
    window=:hanning,
    label=nothing,
    aweighting=false,
    show_oaspl=true,
    lw=1,
    alpha=1,
)
    p = pressure(pth)
    fs = 1 / timestep(pth)

    if xmin === nothing
        try
            xmin = minimum([100, minimum(bpf) * 0.5])
        catch
            xmin = 100
        end
    end
    if xmax === nothing
        xmax = fs / 2.0
    end

    if show_oaspl
        if label === nothing
            label = ""
        else
            label *= ", "
        end
        label *= f"{OASPL(pth):0.2f} dB"
    end

    f, X_dB, _ = narrowband_spectrum(
        pth; type=type, aweighting=aweighting, window=window
    )


    pplt = pyimport("proplot")

    fig, ax = pplt.subplots(figsize=figsize)

    if label === nothing
        ax[1].plot(f, X_dB, lw=lw, alpha=alpha)
    else
        ax[1].plot(f, X_dB, label=label, lw=lw, alpha=alpha)
    end

    # select X_dB > xmin, < xmax
    _X_dB = X_dB[(f.>xmin).&(f.<xmax)]
    if ymax === nothing

        ymax = maximum(_X_dB) + 2
    end
    if ymin === nothing
        ymin = minimum(_X_dB) - 2
    end

    if bpf !== nothing
        for i = 1:n_bpf
            ax[1].vlines(bpf * i, ymin, ymax, color="k", lw=0.5)
        end
    else
        bpf = 62.5
    end

    if ylabel === nothing
        if type == :amplitude
            ylabel = L"$L_p$ (Amp) [dB re 20 μPa]"
        else
            ylabel = L"$L_p$ [dB re 20 μPa / $Δf^{1/2}_{ref}$]"
        end
    end



    ax[1].set(
        title=title,
        xscale="log",
        xlabel=xlabel,
        ylabel=ylabel,
        ylim=(ymin, ymax),
        xlim=(xmin, xmax),
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)

    end

    return fig
end
function plot_narrowband_spectrum(
    list_pth::Vector;
    fname=nothing,
    dpi=300,
    figsize=(6, 3),
    bpf=[],
    n_bpf=10,
    xmin=nothing,
    xmax=nothing,
    ymin=0,
    ymax=nothing,
    title="Spectrum",
    xlabel="f [Hz]",
    ylabel=nothing,
    type=:amplitude,
    window=:hanning,
    label=nothing,
    aweighting=false,
    lw=1,
    alpha=1,
)
    pplt = pyimport("proplot")

    fig, ax = pplt.subplots(figsize=figsize)

    list_ymax = []
    list_ymin = []
    list_fs = []

    for k in 1:length(list_pth)
        fs = 1 / timestep(list_pth[k])
        push!(list_fs, fs)
    end

    if length(bpf) > 0
        plot_bpf = true
    else
        plot_bpf = false
        push!(bpf, 62.5)
    end


    if xmin === nothing
        xmin = minimum([100, minimum(bpf) * 0.5])
    end
    if xmax === nothing
        xmax = maximum(list_fs) / 2.0
    end


    for k in 1:length(list_pth)
        p = pressure(list_pth[k])
        fs = 1 / timestep(list_pth[k])

        f, X_dB, _ = narrowband_spectrum(
            p, fs; type=type, aweighting=aweighting, window=window
        )


        if label === nothing
            ax[1].plot(f, X_dB, lw=lw, alpha=alpha)
        else
            ax[1].plot(f, X_dB, label=label[k], lw=lw, alpha=alpha)
        end

        _X_dB = X_dB[(f.>xmin).&(f.<xmax)]
        if ymax === nothing
            push!(list_ymax, maximum(_X_dB) + 2)
        end
        if ymin === nothing
            push!(list_ymin, minimum(_X_dB) - 2)
        end
    end


    if ymax === nothing
        ymax = maximum(list_ymax)
    end
    if ymin === nothing
        ymin = minimum(list_ymin)
    end

    if plot_bpf
        for k in 1:length(list_pth)
            for i = 1:n_bpf
                ax[1].vlines(bpf[k] * i, ymin, ymax, color=f"C{k-1}", lw=0.5)
            end
        end
    end


    if ylabel === nothing
        if type == :amplitude
            ylabel = L"$L_p$ (Amp) [dB re 20 μPa]"
        else
            ylabel = L"$L_p$ [dB re 20 μPa / $Δf^{1/2}_{ref}$]"
        end
    end

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    list_fs = [1 / timestep(pth) for pth in list_pth]

    ax[1].set(
        title=title,
        xscale="log",
        xlabel=xlabel,
        ylabel=ylabel,
        ylim=(ymin, ymax),
        xlim=(xmin, xmax),
    )

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)

    end

    return fig
end

function plot_history(
    list_pth::Vector;
    title="Time History",
    fname=nothing,
    dpi=300,
    figsize=(6, 2),
    xlabel="t [s]",
    ylabel="p [Pa]",
    label=nothing,
    lw=1,
    alpha=1,
    without_dc_offset=false,
)
    pplt = pyimport("proplot")

    if without_dc_offset
        list_pth = remove_dc_offset.(list_pth)
    end

    fig, ax = pplt.subplots(figsize=figsize)

    for i = 1:length(list_pth)
        p = pressure(list_pth[i])
        t = Vector(time(list_pth[i]))

        if label === nothing
            ax[1].plot(t, p, lw=lw, alpha=alpha)
        else
            ax[1].plot(t, p, label=label[i], lw=lw, alpha=alpha)
        end
    end

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    ax[1].set(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )
    fig

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end

    return fig
end
function plot_history(
    pth::PressureTimeHistory;
    t=nothing,
    title="Time History",
    fname=nothing,
    dpi=300,
    figsize=(6, 2),
    xlabel="t [s]",
    ylabel="p [Pa]",
    label=nothing,
    lw=1,
    alpha=1,
    without_dc_offset=false,
)
    if without_dc_offset
        pth = remove_dc_offset(pth)
    end

    p = pressure(pth)
    fs = 1 / timestep(pth)

    pplt = pyimport("proplot")

    fig, ax = pplt.subplots(figsize=figsize)

    N = length(p)

    Δt = 1.0 / fs
    if t === nothing
        t = Vector(Δt * (0:N-1))
    end

    if label === nothing
        ax[1].plot(t, p, lw=lw, alpha=alpha)
    else
        ax[1].plot(t, p, label=label, lw=lw, alpha=alpha)
    end


    ax[1].set(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end

    return fig
end

function plot_propband_spectrum(
    pth::PressureTimeHistory;
    n=3,
    aweighting=false,
    figsize=(6, 3),
    fname=nothing,
    xmin=62.5,
    xmax=nothing,
    ymin=nothing,
    ymax=nothing,
    title=nothing,
    ylabel=nothing,
    xlabel=nothing,
    label=nothing,
    dpi=300,
    show_oaspl=true,
    lw=1,
    alpha=1,
)
    if show_oaspl
        if label === nothing
            label = ""
        else
            label *= ", "
        end
        label *= f"{OASPL(pth):0.2f} dB"
    end

    cbands, pbs_level, _ = propband_spectrum(pth; aweighting=aweighting)


    if xmax === nothing
        xmax = cbands[end]
    end
    if xmin === nothing
        xmin = cbands[1]
    end
    _pbs_level = pbs_level[(cbands.>xmin).&(cbands.<xmax)]
    if ymax === nothing
        ymax = maximum(_pbs_level) + 2
    end
    if ymin === nothing
        ymin = maximum([minimum(_pbs_level) - 2, 0])
    end


    xlim = (xmin, xmax)
    ylim = (ymin, ymax)

    if title === nothing
        if aweighting
            title = "A-Weighted 1/3-Octave Band Spectrum"
        else
            title = "1/3-Octave Band Spectrum"
        end
    end

    if ylabel === nothing
        if aweighting
            ylabel = L"$L_{pA}$ [dB re 20 μPa]"
        else
            ylabel = L"$L_p$ [dB re 20 μPa]"
        end
    end

    if xlabel === nothing
        xlabel = "f [Hz]"
    end

    pplt = pyimport("proplot")
    fig, ax = pplt.subplots(figsize=figsize)

    if label === nothing
        ax[1].plot(cbands, pbs_level, lw=lw, alpha=alpha)
    else
        ax[1].plot(cbands, pbs_level, label=label, lw=lw, alpha=alpha)
    end

    ax[1].set(
        xscale="log",
        xlim=xlim,
        ylim=ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end

    return fig
end
function plot_propband_spectrum(
    list_pth::Vector;
    n=3,
    aweighting=false,
    figsize=(6, 3),
    fname=nothing,
    xmin=62.5,
    xmax=nothing,
    ymin=nothing,
    ymax=nothing,
    title=nothing,
    ylabel=nothing,
    xlabel=nothing,
    label=nothing,
    dpi=300,
    lw=1,
    alpha=1,
)
    pplt = pyimport("proplot")
    fig, ax = pplt.subplots(figsize=figsize)

    list_xmax = []
    list_xmin = []
    list_ymax = []
    list_ymin = []

    cbands, _, _ = propband_spectrum(list_pth[1]; aweighting=aweighting)

    if xmax === nothing
        push!(list_xmax, cbands[end])
    else
        push!(list_xmax, xmax)
    end
    if xmin === nothing
        push!(list_xmin, cbands[1])
    else
        push!(list_xmin, xmin)
    end
    xmax = maximum(list_xmax)
    xmin = minimum(list_xmin)
    xlim = (xmin, xmax)

    for i = 1:length(list_pth)
        pth = list_pth[i]
        cbands, pbs_level, _ = propband_spectrum(pth; aweighting=aweighting)

        if label === nothing
            ax[1].plot(cbands, pbs_level, lw=lw, alpha=alpha)
        else
            ax[1].plot(cbands, pbs_level, label=label[i], lw=lw, alpha=alpha)
        end

        _pbs_level = pbs_level[(cbands.>xmin).&(cbands.<xmax)]
        if ymax === nothing
            push!(list_ymax, maximum(_pbs_level) + 2)
        else
            push!(list_ymax, ymax)
        end
        if ymin === nothing
            push!(list_ymin, maximum([minimum(_pbs_level) - 2, 0]))
        else
            push!(list_ymin, ymin)
        end
    end


    ylim = (minimum(list_ymin), maximum(list_ymax))

    if title === nothing
        if aweighting
            title = "A-Weighted 1/3-Octave Band Spectrum"
        else
            title = "1/3-Octave Band Spectrum"
        end
    end

    if ylabel === nothing
        if aweighting
            ylabel = L"$L_{pA}$ [dB re 20 μPa]"
        else
            ylabel = L"$L_p$ [dB re 20 μPa]"
        end
    end

    if xlabel === nothing
        xlabel = "f [Hz]"
    end

    ax[1].set(
        xscale="log",
        xlim=xlim,
        ylim=ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title
    )

    if label !== nothing
        ax[1].legend(ncols=1)
    end

    if fname !== nothing
        fig.savefig(fname, dpi=dpi)
    end

    return fig
end


# ==============================================================================
# MSST 
# ==============================================================================

function plot_msst(
    MSST::Matrix,
    t::Vector;
    y_min=nothing,
    y_max=nothing,
    noise_floor_dB=0.0,
    figsize=(6, 3.7),
    rasterized=true,
    fname=nothing,
    tight_layout=true,
)
    fs = 1 / (t[2] - t[1])

    T_msst_plot = 20.0 * log10.(abs.(MSST) / P_REF)
    T_msst_plot_min = nanminimum(T_msst_plot)
    T_msst_plot = ifelse.(isinf.(T_msst_plot), noise_floor_dB, T_msst_plot)


    f_msst = Vector(LinRange(0, fs / 2, size(MSST, 1)))

    fig, ax = plt.subplots(figsize=figsize, tight_layout=tight_layout)

    if y_min === nothing
        y_min = f_msst[1]
    end
    if y_max === nothing
        y_max = f_msst[end]
    end

    mesh = ax.pcolorfast(
        t,
        f_msst,
        T_msst_plot,
        vmin=maximum([T_msst_plot_min, noise_floor_dB]),
        vmax=maximum(T_msst_plot),
        rasterized=rasterized,
    )

    cbar = fig.colorbar(mesh, ax=ax, label="dB")
    ax.set(
        ylim=(y_min, y_max),
        title="MSST",
        xlabel="t [s]",
        ylabel="f [Hz]",
    )

    if fname !== nothing
        fig.savefig(fname)
    end

    return fig
end
function plot_msst(
    pth::PressureTimeHistory;
    y_min=nothing,
    y_max=nothing,
    noise_floor_dB=0.0,
    figsize=(6, 3.7),
    rasterized=true,
    fname=nothing,
    tight_layout=true,
    fs_resample=nothing, n_iterations=4, length_window=256, max_ram=3.0
)
    MSST, t, _ = msst(pth; fs_resample=fs_resample, n_iterations=n_iterations, length_window=length_window, max_ram=max_ram)
    return plot_msst(
        MSST, t; 
        y_min=y_min, y_max=y_max, noise_floor_dB=noise_floor_dB, 
        figsize=figsize, rasterized=rasterized, fname=fname, tight_layout=tight_layout
    )
end
function plot_msst(
    p::Vector,
    t::Vector;
    y_min=nothing,
    y_max=nothing,
    noise_floor_dB=0.0,
    figsize=(6, 3.7),
    rasterized=true,
    fname=nothing,
    tight_layout=true,
    fs_resample=nothing, n_iterations=4, length_window=256, max_ram=3.0
)
    MSST, t, _ = msst(p, t; fs_resample=fs_resample, n_iterations=n_iterations, length_window=length_window, max_ram=max_ram)
    return plot_msst(
        MSST,t; 
        y_min=y_min, y_max=y_max, noise_floor_dB=noise_floor_dB, 
        figsize=figsize, rasterized=rasterized, fname=fname, tight_layout=tight_layout
    )
end