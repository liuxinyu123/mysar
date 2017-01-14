function y = ftx(x) % range fft
y = fftshift(fft(fftshift(x.'))).';