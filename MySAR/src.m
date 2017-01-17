% Liu YK 2017/1/14

function y = src(signal_rfaf, azimuth_frequence, range_frequence, vr, nearest_range_ref, fc )

c = 3e8;
lambda = c / fc;
[nan, nrn] = size(signal_rfaf);
migrate_coff = sqrt(1 - (lambda * azimuth_frequence / 2 / vr).^2);
k_src = 2 * vr^2 * fc^3 / c / nearest_range_ref ./ (azimuth_frequence' * ones(1,nrn)).^2  .* (migrate_coff' * ones(1,nrn)).^3;
src_function = exp(-1i * pi * (ones(nan,1) * range_frequence).^2  ./ k_src);
y = iftx(signal_rfaf .* src_function);
