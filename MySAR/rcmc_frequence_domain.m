% Liu YK 2017/1/13

function y = rcmc_frequence_domain(signal_rd, nearest_range_vector, azimuth_frequence, fr, lambda, vr) 

c = 3e8; % light speed
[nan, nrn] = size(signal_rd);
y = zeros(nan, nrn);

nearest_range_ref = nearest_range_vector(round(size(nearest_range_vector,2)  / 2));
alpha = [-nrn/2 : nrn/2 - 1] / nrn ;
% alpha = (circshift(alpha', shift_n))';
delta_range = c/2/fr; 
n_rcm=(1./sqrt(1-(azimuth_frequence*lambda/2/vr).^2)'-1) * nearest_range_ref / delta_range;
base_n_rcm = (1./sqrt(1-(azimuth_frequence(round(nan/2))*lambda/2/vr).^2)-1)*nearest_range_ref/delta_range;
n_rcm = n_rcm - base_n_rcm;

for i = 1:nan
    y(i,:) = fftshift(ifft(fftshift(fft(fftshift(signal_rd(i,:))) .* exp(1i*2*pi*n_rcm(i)*alpha))));
end 