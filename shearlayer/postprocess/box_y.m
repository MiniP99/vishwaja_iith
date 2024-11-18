function [f_fil_y_indices] = box_y(f,Ny)
[nx, ny, nz] = size(f);

f_fil_y = f;
num_indices = floor(ny / (2 * Ny));
f_fil_y_indices = zeros(nx, num_indices, nz);


f_pad= zeros(nx, ny + 2 * Ny, nz);
f_pad(:, Ny+1:Ny+ny, :) = f;  %Copy interior
for j = 1:Ny
    f_pad(:, Ny+1-j, :) = f_pad(:, Ny+1+j, :);   %left
    f_pad(:, Ny+ny+j, :) = f_pad(:, Ny+ny-j, :); %right
end

index = 1;
%Apply box filter
for j = Ny+1:Ny+ny
    if mod(j-(Ny+1), 2*Ny) == 0  %%Calculating only at 2*Ny th point
       f_fil_y(:, j-Ny, :) = mean(f_pad(:, j-Ny:j+Ny, :),2);
       f_fil_y_indices(:, index, :) = f_fil_y(:, j-Ny, :);
       index=index+1;
    end 
end

end
