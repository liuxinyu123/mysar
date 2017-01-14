function y = rcmc_frequence_domain(signal_rd, nearest_range_vector, azimuth_frequence, fr, lambda, vr) 

c = 3e8; % light speed
[nan, nrn] = size(signal_rd);
% y = zeros(nan, nrn);
alpha = [-nrn/2 : nrn/2 - 1] / nrn;

delta_rcm = (lambda / vr)^2 / 8 * (ones(nan,1) * nearest_range_vector) .* (azimuth_frequence' * ones(1,nrn));
n_rcm = 2 * delta_rcm / c * fr;
corr_function = exp(1i * 2 * pi * (ones(nan,1) * alpha ) .* n_rcm);
y = iftx(ftx(signal_rd) .* corr_function);