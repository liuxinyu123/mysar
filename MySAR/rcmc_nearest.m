% Liu YK 2017/1/12

function y = rcmc_nearest(signal_rd, nearest_range_vector, azimuth_frequence, fr, lambda, vr) 

c = 3e8;
[nan, nrn] = size(signal_rd);
y = zeros(nan, nrn);

window = waitbar(0, 'Nearest neighbor interpolation'); 
for i = 1 : nan
    for j = 1 : nrn
        delta_range_rcm = (lambda / vr)^2 / 8 * nearest_range_vector(j) * azimuth_frequence(i)^2;
        n_rcm = 2 * delta_range_rcm / c * fr;
        n_rcm_decimal = n_rcm - floor(n_rcm);
        
        if j + n_rcm > nrn
%             y(i, j) = y(i, round(nrn / 2));
            y(i, j) = 0;
        else
            if n_rcm_decimal < 0.5
                y(i, j) = signal_rd(i, j + floor(n_rcm));
            else
                y(i, j) = signal_rd(i, j + ceil(n_rcm));
            end
        end
    end
    waitbar(i / nan);
end
close(window);