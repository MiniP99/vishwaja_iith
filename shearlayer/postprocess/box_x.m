function [f_fil_x_indices] = box_x(f,Nx)

%Lx=150;      Ly=200;      Lz=75; %Mc=0.2,0.4
%nx_les=180;  ny_les=240;  nz_les=84;
%nx_dns=1024; ny_dns=1448; nz_les=512;
%N=LES filter width/DNS filter width = (Lx/nx_les)/(Lx/nx_dns) Eg:1024/180=6 
%Nx=No. of points to be considered at one side of point-i = N/2 = 6/2 = 3 
%i-3,i-2,i-1,i,i+1,i+2,i+3

[nx, ny, nz] = size(f);

f_fil_x = f;

num_indices = floor(nx / (2 * Nx));
f_fil_x_indices = zeros(num_indices, ny, nz);

%Padded array for boundary conditions
f_pad = zeros(nx + 2 * Nx, ny, nz);
f_pad(Nx+1:Nx+nx, :, :)  = f;                         %Copy interior points
for i = 1:Nx
    f_pad(Nx+1-i, :, :)  = f_pad(Nx+1-i+nx, :, :);    %periodic conditions-left
    f_pad(Nx+nx+i, :, :) = f_pad(Nx+nx+i-nx, :, :);   %periodic conditions-right
end

index = 1;
%Apply box filter
for i = Nx+1:Nx+nx
    if mod(i-(Nx+1), 2*Nx) == 0  %%Calculating only at 2*Nx th point
        f_fil_x(i-Nx,:,:)    = mean(f_pad(i-Nx:i+Nx, :, :),1); %filter in x-direc
        f_fil_x_indices(index, :, :) = f_fil_x(i-Nx, :, :);
        index=index+1;
    end 
end

end
