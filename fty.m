function y = fty(x) % azimuth fft
y = fftshift(fft(fftshift(x)));