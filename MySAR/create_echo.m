% Liu YK 2017/1/11
function echo = create_echo(scene_center_range, range_width, azimuth_width,...
    targets, prf, fr, lambda, kr, tp, vr, theta, beta )

c = 3e8;  % the light speed
targets_numbers = size(targets, 1); %the numbers of targets

scene_center_x = scene_center_range * cos(theta); % scene center location in range direction
scene_center_y = scene_center_range * sin(theta); % scene center location in azimuth direction

range_left = scene_center_x - range_width / 2; % the left boundary of range
range_right = scene_center_x + range_width / 2; % the right boundary of range
azimuth_bottom = scene_center_y - azimuth_width / 2; % the bottom boundary of azimuth 
azimuth_top = scene_center_y + azimuth_width / 2; % the top boundary of azimuth 

azimuth_begin = (azimuth_bottom - range_right * tan(theta + beta / 2)) / vr; % the beginning time in azimuth direction
azimuth_end = (azimuth_top - range_left * tan(theta - beta / 2)) / vr; % the ending time in azimuth direction
nan = round((azimuth_end - azimuth_begin) * prf); % the sampling numbers in azimuth direcion
% nan = 2 ^ nextpow2(nan); % for fft
slow_time = [-nan / 2 : nan / 2 - 1] / prf;  % slow time vector 

range_near = range_left / cos(theta - beta / 2);  % the nearest range in range direction
range_far = range_right / cos(theta + beta / 2); % the farthest range in range direction
range_mid = (range_far + range_near) / 2; % the middle range in range direction
nrn = round((2 * (range_far - range_near) / c + tp) * fr);  % the sampling numbers in range direction
% nrn = 2 ^ nextpow2(nrn); % for fft
fast_time = [-nrn / 2 : nrn / 2 - 1] / fr + 2 * range_mid / c; % fast time vector

echo = zeros(nan, nrn);  %initialize echo

for i = 1 : targets_numbers
    x = targets(i,1);
    y = targets(i,2);
    rcs = targets(i,3);
    
    beam_center_cross_time = (y - scene_center_y + (scene_center_x - x) * tan(theta)) / vr;
    beam_cross_begin_time = beam_center_cross_time - x * (tan(theta + beta / 2) - tan(theta)) / vr;
    beam_cross_end_time = beam_center_cross_time + x * (tan(theta) - tan(theta - beta / 2)) / vr;
    
    instant_range = sqrt(x^2 + (y - vr * slow_time).^2); 
    tau = ones(nan,1) * fast_time - 2 * instant_range' * ones(1,nrn) / c; % fast time - delay time
    range_phase = pi * kr * tau.^2;
    azimuth_phase = -4 * pi / lambda * (instant_range' * ones(1,nrn));
    range_limit = (abs(tau) < tp / 2);
    azimuth_limit = (((slow_time > beam_cross_begin_time) &...
        (slow_time < beam_cross_end_time))' * ones(1,nrn));
    
    echo = echo + rcs * exp(1i * range_phase) .* exp(1i * azimuth_phase) .* range_limit .* azimuth_limit;
end

    

