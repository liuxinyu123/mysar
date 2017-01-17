function y = rcmc_sinc(signal_rd, nearest_range_vector, azimuth_frequence, fr, lambda, vr)

c = 3e8; % light speed
sinc_core_length = 32;
[nan, nrn] = size(signal_rd);
y = zeros(nan, nrn);

% delta_range = c/2/fr;

% n_rcm=(1./sqrt(1-(f_azimuth*lambda/2/V).^2)'-1)*R0/delta_range;
% base_n_rcm = (1./sqrt(1-(f_azimuth(Na/2)*lambda/2/V).^2)-1)*R0/delta_range;
% n_rcm = n_rcm - base_n_rcm;
for i = 1 : nrn
    nearest_range = nearest_range_vector(i);
    delta_rcm = lambda^2 * nearest_range / 8 / vr^2 * azimuth_frequence.^2;
    n_rcm = 2 * delta_rcm / c * fr;
end

win = waitbar(0,'Sinc²åÖµ');
for i = 1:nan
    for j = sinc_core_length:nrn
        rcm = n_rcm(i);
        if rcm > 0
            delta_n_rcm = rcm - floor(rcm);
        else            
            delta_n_rcm = rcm - ceil(rcm);
        end
        
        for k = -sinc_core_length/2:sinc_core_length/2 - 1
            if (k+j+ceil(rcm) > nrn) || (k+j+floor(rcm) <= 0)
                y(i,j) = y(i,j) + signal_rd(i,nrn) * sinc(k+rcm);
            else
                y(i,j) = y(i,j) + signal_rd(i,j+floor(rcm)+k) * sinc(delta_n_rcm+k);        
            end
        end
    end
    waitbar(i/nan);
end
close(win);



