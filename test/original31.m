function [Ts, tfr, omega2] = MSST_Y_new(x, hlength, num);
    % Computes the MSST (Ts)  of the signal x.
    % Expression (31)-based Algorithm.
    % INPUT
    %    x      :  Signal needed to be column vector.
    %    hlength:  The hlength of window function.
    %    num    :  iteration number.
    % OUTPUT
    %    Ts     :  The SST
    %    tfr     :  The STFT

    [xrow, xcol] = size(x);

    if (xcol ~= 1),
        error('X must be column vector');
    end;

    if (nargin < 3),
        error('At least 3 parameter is required');
    end

    hlength = hlength + 1 - rem(hlength, 2);
    ht = linspace(-0.5, 0.5, hlength);
    ht = ht';

    % Gaussian window
    h = exp(-pi / 0.32 ^ 2 * ht .^ 2);

    [hrow, ~] = size(h);
    Lh = (hrow - 1) / 2;

    N = xrow;
    t = 1:xrow;

    [~, tcol] = size(t);

    n_half = round(N / 2);

    tfr = zeros(n_half, tcol);
    omega = zeros(n_half, tcol - 1);
    omega2 = zeros(n_half, tcol);
    Ts = zeros(n_half, tcol);

    % ---------------------------------
    % Block 0
    % ---------------------------------
    for icol = 1:tcol,
        ti = t(icol);
        tau = -min([n_half - 1, Lh, ti - 1]):min([n_half - 1, Lh, N - ti]);
        indices = rem(N + tau, N) + 1;
        rSig = x(ti + tau, 1);
        tfr(indices, icol) = rSig .* conj(h(Lh + 1 + tau));
    end;

    tfr = fft(tfr); % columns as vectors for FFT
    tfr = tfr(1:n_half, :);

    % ---------------------------------
    % Block 1
    % ---------------------------------
    for i = 1:n_half
        omega(i, :) = diff(unwrap(angle(tfr(i, :)))) * N / 2 / pi;
    end

    omega(:, end + 1) = omega(:, end);
    omega = round(omega);

    [neta, nb] = size(tfr);

    % ---------------------------------
    % Block 2
    % ---------------------------------
    if num > 1

        for kk = 1:num - 1

            for b = 1:nb

                for eta = 1:neta
                    k = omega(eta, b);

                    if k >= 1 && k <= neta
                        omega2(eta, b) = omega(k, b);
                    end

                end

            end

            omega = omega2;
        end

    else
        omega2 = omega;
    end

    % ---------------------------------
    % Block 3
    % ---------------------------------
    for b = 1:nb %time
        % Reassignment step
        for eta = 1:neta %frequency

            if abs(tfr(eta, b)) > 0.0001 %you can set much lower value than this.
                k = omega2(eta, b);

                if k >= 1 && k <= neta
                    Ts(k, b) = Ts(k, b) + tfr(eta, b);
                end

            end

        end

    end

    tfr = tfr / (N / 2);
    Ts = Ts / (N / 2);
end
