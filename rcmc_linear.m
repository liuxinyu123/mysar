function y = rcmc_linear(signal_rd, nearest_range_vector, azimuth_frequence, fr, lambda, vr)

c = 3e8; %light speed
[nan, nrn] = size(signal_rd); 
y = zeros(nan, nrn); % initialize output

window = waitbar(0, 'Linear interpolation');
for i = 1 : nrn
    for j = 1 : nan
        delta_range_rcm = (lambda / vr)^2 / 8 * nearest_range_vector(i) * azimuth_frequence(j)^2;
        delta_n_rcm = 2 * delta_range_rcm / c * fr;
        delta_n_rcm_integer = floor(delta_n_rcm);
        delta_n_rcm_decimal = delta_n_rcm - floor(delta_n_rcm);
        
        if i + delta_n_rcm > nrn
            y(j, i) = 0;
        else
            y(j, i) = signal_rd(j, i + delta_n_rcm_integer) + delta_n_rcm_decimal *...
                (signal_rd(j, ceil(delta_n_rcm)) -  signal_rd(j, floor(delta_n_rcm)));
        end
    end
    waitbar(i / nrn);
end
close(window);
    