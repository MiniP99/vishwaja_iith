function dfdy = ddy(f,ky)
ny = size(f,2);
fhat = fft(f,[],2);
ghat = 1i*ky.*fhat;
ghat(:,ny/2+1,:) = 0;
dfdy = real(ifft(ghat,[],2));
