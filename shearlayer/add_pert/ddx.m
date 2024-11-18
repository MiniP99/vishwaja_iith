function dfdx = ddx(f,kx)
nx = size(f,1);
fhat = fft(f,[],1);
ghat = 1i*kx.*fhat;
ghat(nx/2+1,:,:) = 0;
dfdx = real(ifft(ghat,[],1));
