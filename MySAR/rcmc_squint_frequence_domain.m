% Liu YK 2017/1/16

function y = rcmc_squint_frequence_domain(signal_rd, nearest_range_vector, azimuth_frequence, fr, lambda, vr) 

c = 3e8;
[nan, nrn] = size(signal_rd);
y = zeros(nan, nrn);

range_frequence = [-nrn/2 : nrn/2 - 1] / nrn * fr;
% range_frequence = [0:nrn-1] / nrn * fr;
% delta_range = c / 2 / fr;
magrate_coff = sqrt(1 - (lambda * azimuth_frequence / 2 / vr).^2);
magrate_coff_ref = sqrt(1 - (lambda * azimuth_frequence(round(size(azimuth_frequence,2)/2)) / 2 / vr).^2);
% nearest_range_ref = nearest_range_vector(round(size(nearest_range_vector,2)  / 2));
n_rcm = (ones(nan,1) * nearest_range_vector) .* ((1 ./ magrate_coff - 1 / magrate_coff_ref)' * ones(1,nrn));
delay_time_rcm = 2 * n_rcm / c;
 
window_function = kaiser(nrn,2.5);
corr_function = exp(1i * 2 * pi .* delay_time_rcm .* (ones(nan,1) * range_frequence)) .* (ones(nan,1) * window_function.') ;
y = iftx(ftx(signal_rd) .* corr_function ) ;
% for i = 1:nan
%     y(i,:) = fftshift(ifft(fftshift(fft(fftshift(signal_rd(i,:))) .* exp(1i*2*pi*n_rcm(i)*alpha))));
% end 