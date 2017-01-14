function y = iftx(x) % range ifft
y = fftshift(ifft(fftshift(x.'))).';