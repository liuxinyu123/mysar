% Liu YK 2017/1/11

function y = ifty(x) % azimuth ifft
y = fftshift(fft(fftshift(x)));