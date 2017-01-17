% Liu YK 2017/1/11

function y = ftx(x) % range fft
y = fftshift(fft(fftshift(x.'))).';