% Liu YK 2017/1/11
% broadside SAR and low squint

clc, clear, close all;

%=====================paremeter setting====================%
c = 3e8; %light speed
fc = 5.3e9; % C waveband
lambda = c / fc; %wave length
vr = 150; % radar speed
theta = 45 / 180 * pi; % squint angle
tp = 2.5e-6; % pulse duration
kr = 20e12; % chirp rate
br = abs(kr * tp); % signal bandwidth in range direction
scene_center_range = 20e3;  
d = 3.4; % radar length in azimuth direction
beta = lambda / d; % wave beam width in angle

fr = 60e6; %sampling rate in range direction
ka = -2*vr^2*cos(theta)^2/scene_center_range/lambda; % doppler rate
lsar = scene_center_range * cos(theta) * (tan(theta + beta / 2) -...
    tan(theta - beta / 2)); % synthetic aperture length
tsar = lsar / vr; % synthetic aperture time
ba = abs(ka * tsar);
prf = 120;

range_width = 200;
azimuth_width = 200;
% xc = scene_center_range * cos(theta);
% yc = scene_center_range * sin(theta);
xc = scene_center_range * cos(theta);
yc = 0;

k1 = 0.4;
k2 = 0.6;
targets = [ xc     yc       1
%             xc + k1 * range_width / 2    yc              1
%             xc                           yc - tan(theta) * k1 * range_width / 2    1
            xc + k1 * range_width / 2    yc + k2 * azimuth_width / 2      1
            xc + k1 * range_width / 2    yc - k2 * azimuth_width / 2      1
            xc - k1 * range_width / 2    yc + k2 * azimuth_width / 2      1
            xc - k1 * range_width / 2    yc - k2 * azimuth_width / 2      1
            ];
            
disp('plot scene targets...');
figure;
plot(targets(:,1),targets(:,2),'r*');
xlabel('range direction -(m)');
ylabel('azimuth direction -(m)');
xlim([xc - range_width/2 xc + range_width/2]);
ylim([yc - azimuth_width/2 yc + azimuth_width/2]);
title('scene targets');

%==================== generate echo===================%
disp('generate the radar echo signal...');
[echo, nearest_range_vector] = create_echo_v2(scene_center_range, range_width,...
    azimuth_width, targets, prf, fr, lambda, kr, tp, vr, theta, beta);
disp('plotting the graph of echo signal...');
figure;
colormap('gray');
subplot(221);
imagesc(real(echo));
xlabel('range direction(sampling points)');
ylabel('azimuth direction(sampling points)');
title('real section of raw echo signal');
subplot(222);
imagesc(imag(echo));
xlabel('range direction(sampling points)');
ylabel('azimuth direction(sampling points)');
title('imaginary section of raw echo signal');
subplot(223);
imagesc(255 - abs(echo));
xlabel('range direction(sampling points)');
ylabel('azimuth direction(sampling points)');
title('attitude of raw echo signal');
subplot(224);
imagesc(angle(echo));
xlabel('range direction(sampling points)');
ylabel('azimuth direction(sampling points)');
title('phase of raw echo signal');


[nan, nrn] = size(echo);
doppler_frequence_center = 2 * vr * sin(theta) / lambda;
n_blur = round(doppler_frequence_center / prf);
doppler_frequence_center_base = doppler_frequence_center - n_blur * prf;
shift_n = round(doppler_frequence_center_base / prf * nan);
range_frequence = [-nrn/2 : nrn/2 - 1] / nrn * fr;
azimuth_frequence = doppler_frequence_center + [-nan/2 : nan/2 - 1] / nan * prf;
azimuth_frequence = (circshift(azimuth_frequence', shift_n))';

%========================= range compression=========================%
disp('range pulse compression');
window_function = kaiser(nrn, 2.5).'; % kaiser window
range_compression_function = exp(1i * pi / kr * range_frequence.^2) .* window_function;
signal_range_compression = iftx(ftx(echo) .* (ones(nan,1) * range_compression_function));

disp('plotting the graph of signal after range compression');
figure;
colormap('gray');
subplot(121);
imagesc(real(signal_range_compression));
xlabel('Range direction(sampling points)');
ylabel('Azimuth direction(sampling points)');
title('real section of signal after range compression');
subplot(122);
imagesc(255 - abs(signal_range_compression));
xlabel('Range direction(sampling points)');
ylabel('Azimuth direction(sampling points)');
title('attitude of signal after range compression');

disp('azimuth fft');
signal_rtaf = fty(signal_range_compression);
% signal_rfaf = ftx(signal_rtaf);
% signal_src = src(signal_rfaf, azimuth_frequence, range_frequence, vr, scene_center_range, fc );
disp('azimuth fft finished');
  
disp('plotting the graph of signal after azimuth fft');
figure;
colormap('gray');
subplot(121);
imagesc(real(signal_rtaf));
title('real section of signal after azimuth fft');
xlabel('Range direction(sampling points)');
ylabel('Azimuth direction(Hz)');
% ylim([max(azimuth_frequence) min(azimuth_frequence)]);
subplot(122);
imagesc(255 - abs(signal_rtaf));
title('attitude of signal after azimuth fft');
xlabel('Range direction(sampling points)');
ylabel('Azimuth direction(Hz)');
% ylim([max(azimuth_frequence) min(azimuth_frequence)]); 

%=================  SRC ============================%
disp('src...');
signal_src = src(ftx(signal_rtaf), azimuth_frequence, range_frequence, vr,  scene_center_range * cos(theta), fc );
disp('src finished');
signal_rtaf = signal_src;

% ========================RCMC====================================
disp('rcmc...');
% signal_rcmc = rcmc_nearest(signal_rtaf, nearest_range_vector, azimuth_frequence, fr, lambda, vr);
%  signal_rcmc = rcmc_nearest_v2(signal_rtaf, nearest_range_vector, azimuth_frequence, fr, lambda, vr, shift_n);
% signal_rcmc = rcmc_frequence_domain(signal_rtaf, nearest_range_vector, azimuth_frequence, fr, lambda, vr);
signal_rcmc = rcmc_squint_frequence_domain(signal_rtaf, nearest_range_vector, azimuth_frequence, fr, lambda, vr);
disp('rcmc finished');

disp('plotting the graph of the signal after rcmc');
figure;
colormap('gray');
subplot(121);
imagesc(real(signal_rcmc));
title('real section of signal after rcmc');
xlabel('Range direction(sampling points)');
ylabel('Azimuth direction(Hz)');
subplot(122);
imagesc(255 - abs(signal_rcmc));
title('attitude of signal after rcmc');
xlabel('Range direction(sampling points)');
ylabel('Azimuth direction(Hz)');


%========================== azimuth compression============================
disp('azimuth compression');
window_function = kaiser(nan, 5);
window_function = (circshift(window_function,shift_n));
% azimuth_frequence = doppler_frequence_center + [-nan/2 : nan/2 - 1] / nan * prf;
azimuth_compression_function = exp(1i * pi / ka * azimuth_frequence.^2).' .* window_function;
signal_final = ifty(signal_rcmc .* (azimuth_compression_function * ones(1, nrn)));
disp('azimuth compression finished');

disp('plotting the final signal of the gargets');
figure;
colormap('gray');
subplot(121);
mesh(abs(signal_final));
title('the targets in scene');

subplot(122);
imagesc(255 - abs(signal_final));
title('the targets in scene');
xlabel('Range direction(sampling points)');
ylabel('Azimuth direction(sampling points)');