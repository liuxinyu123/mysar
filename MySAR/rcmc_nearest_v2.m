% Liu YK 2017/1/17

function y = rcmc_nearest_v2(signal_rd, nearest_range_vector, azimuth_frequence, fr, lambda, vr, shift_n) 

window = waitbar(0, 'Nearest neighbor interpolation'); 

c = 3e8;
[nan, nrn] = size(signal_rd);
y = zeros(nan, nrn);

nan_boundary = shift_n + round(nan / 2);
delta_range_rcm = (lambda / vr)^2 / 8 * (ones(nan,1) * nearest_range_vector) .* (azimuth_frequence'.^2 * ones(1,nrn));
azimuth_frequence_ref = azimuth_frequence(nan_boundary);
nearest_range_ref = nearest_range_vector(round(size(nearest_range_vector,2) / 2));
delta_range_ref =  (lambda / vr)^2 / 8 * nearest_range_ref * azimuth_frequence_ref^2;
delta_range = c / 2 / fr;
n_rcm = (delta_range_rcm - delta_range_ref * ones(nan, nrn)) / delta_range;
% n_rcm_decimal = n_rcm - floor(n_rcm);

for i = shift_n : nan_boundary
    for j = nrn : 1
        if j + n_rcm(i,j) < 0 
            y(i,j) = 0;
        else
            n_rcm_decimal = n_rcm(i,j) - ceil(n_rcm(i,j));
            if abs(n_rcm_decimal) < 0.5
                y(i,j) = signal_rd(i,j + ceil(n_rcm(i,j)));
            else
                y(i,j) = signal_rd(i, j + floor(n_rcm(i,j)));
            end
        end
    end
    waitbar((i - shift_n) / nan);
end
      
for i = 1 : shift_n - 1
    for j = 1 : nrn
        if j + n_rcm(i,j) > nrn
            y(i,j) = 0;
        else
            n_rcm_decimal = n_rcm(i,j) - floor(n_rcm(i,j));
            if n_rcm_decimal < 0.5
                y(i,j) = signal_rd(i,j+floor(n_rcm(i,j)));
            else
                y(i,j) = signal_rd(i,j+ceil(n_rcm(i,j)));
            end
        end
    end
    waitbar(1/2 + i / nan);
end

for i = nan_boundary + 1 : nan
    for j = 1 : nrn
        if ceil(j + n_rcm(i,j)) >= nrn
            y(i,j) = 0;
        else
            n_rcm_decimal = n_rcm(i,j) - floor(n_rcm(i,j));
            if n_rcm_decimal < 0.5
                y(i,j) = signal_rd(i,j+floor(n_rcm(i,j)));
            else
                y(i,j) = signal_rd(i,j+ceil(n_rcm(i,j)));
            end
        end
    end
    waitbar(1 - (nan - i) / nan);
end

close(window);