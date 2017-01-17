% Liu YK 2017/1/16

function y =  rcmc_squint_nearest(signal_rd, nearest_range_vector, azimuth_frequence, doppler_frequence_center, fr, lambda, vr) 

c = 3e8;
[nan, nrn] = size(signal_rd);
y = zeros(nan, nrn);

window = waitbar(0, 'Nearest neighbor interpolation'); 
for i = 1 : nan
    for j = 1 : nrn
        migrate_coff_ref = sqrt(1 - (lambda * doppler_frequence_center / 2 / vr)^2);
        migrate_coff = sqrt(1 - (lambda * azimuth_frequence(i) / 2 / vr)^2);
        delta_range_rcm = nearest_range_vector(j) * (1 / migrate_coff  - 1 / migrate_coff_ref);
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