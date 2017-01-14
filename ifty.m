function y = ifty(x) % azimuth ifft
y = fftshift(fft(fftshift(x)));