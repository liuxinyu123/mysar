% Liu YK 2017/1/11

function y = fty(x) % azimuth fft
y = fftshift(fft(fftshift(x)));