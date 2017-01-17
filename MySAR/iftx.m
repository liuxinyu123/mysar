% Liu YK 2017/1/11

function y = iftx(x) % range ifft
y = fftshift(ifft(fftshift(x.'))).';