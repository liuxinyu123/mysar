% Liu YK 2017/1/12
function [echo, nearest_range_vector] = create_echo_v2(scene_center_range,...
    range_width, azimuth_width, targets, prf, fr, lambda, kr, tp, vr, theta, beta )

c = 3e8;  % the light speed
targets_numbers = size(targets, 1); %the numbers of targets

scene_center_x = scene_center_range * cos(theta); % scene center location in range direction
% scene_center_y = 0; % scene center location in azimuth direction
zero_time_x = scene_center_x; % beam center x location in zero time 
zero_time_y = scene_center_range * sin(theta); % beam center y location in zero time 

range_left = scene_center_x - range_width / 2; % the left boundary of range
range_right = scene_center_x + range_width / 2; % the right boundary of range
% azimuth_bottom = scene_center_y - azimuth_width / 2; % the bottom boundary of azimuth 
% azimuth_top = scene_center_y + azimuth_width / 2; % the top boundary of azimuth 

azimuth_begin = (-(azimuth_width / 2 + range_right * (tan(theta + beta / 2) -tan(theta)))...
    - (zero_time_y + range_width / 2 * tan(theta))) / vr; % the beginning time in azimuth direction
azimuth_end = ((azimuth_width / 2 + range_left * (tan(theta) - tan(theta - beta / 2)))...
    - (zero_time_y - range_width / 2 * tan(theta))) / vr; % the ending time in azimuth direction
nan = round((azimuth_end - azimuth_begin) * prf); % the sampling numbers in azimuth direcion
% nan = 2 ^ nextpow2(nan); % for fft
% slow_time = [-nan / 2 : nan / 2 - 1] / prf;  % slow time vector 
slow_time = linspace(azimuth_begin,azimuth_end,nan);

range_near = range_left / cos(theta - beta / 2);  % the nearest range in range direction
range_far = range_right / cos(theta + beta / 2); % the farthest range in range direction
range_mid = (range_far + range_near) / 2; % the middle range in range direction
nrn = round((2 * (range_far - range_near) / c + tp) * fr);  % the sampling numbers in range direction
% nrn = 2 ^ nextpow2(nrn); % for fft
fast_time = [-nrn / 2 : nrn / 2 - 1] / fr + 2 * range_mid / c; % fast time vector
nearest_range_vector = c * fast_time / 2 * cos(theta);

echo = zeros(nan, nrn);  %initialize echo

for i = 1 : targets_numbers
    x = targets(i,1);
    y = targets(i,2);
    rcs = targets(i,3);
    nearest_range_absolute_time = y / vr;
    
    beam_center_cross_time = (y - zero_time_y + (zero_time_x - x) * tan(theta)) / vr;
    beam_cross_begin_time = beam_center_cross_time - x * (tan(theta + beta / 2) - tan(theta)) / vr;
    beam_cross_end_time = beam_center_cross_time + x * (tan(theta) - tan(theta - beta / 2)) / vr;
    
    instant_range = sqrt(x^2 + (vr * (slow_time - nearest_range_absolute_time)).^2); 
    tau = ones(nan,1) * fast_time - 2 * instant_range' * ones(1,nrn) / c; % fast time - delay time
    range_phase = pi * kr * tau.^2;
    azimuth_phase = -4 * pi / lambda * (instant_range' * ones(1,nrn));
    range_limit = (abs(tau) < tp / 2);
    azimuth_limit = (((slow_time > beam_cross_begin_time) &...
        (slow_time < beam_cross_end_time))' * ones(1,nrn));
    
    echo = echo + rcs * exp(1i * range_phase) .* exp(1i * azimuth_phase) .* range_limit .* azimuth_limit;
end

    

