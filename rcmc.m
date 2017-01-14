function y = rcmc(signal_rd, nearest_range_vector, azimuth_frequence, fr, lambda, vr)

c = 3e8; % light speed
[nan, nrn] = size(signal_rd);
y = zeros(nan, nrn);


for i = 1 : nrn
    nearest_range = nearest_range_vector(i);
    delta_rcm = lambda^2 * nearest_range / 8 / vr^2 * azimuth_frequence.^2;
    n_rcm = 2 * delta_rcm / c * fr;
    xx = 1 : nan;
    xxq = i + n_rcm;
    result = interp1(xx, signal_rd(:,i).', xxq);
    y(:,i) = result.';
end
