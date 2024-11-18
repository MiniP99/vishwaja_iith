function dfdz = ddz(f,kz)
nz = size(f,3);
fhat = fft(f,[],3);
ghat = 1i*kz.*fhat;
ghat(:,:,nz/2+1) = 0;
dfdz = real(ifft(ghat,[],3));
