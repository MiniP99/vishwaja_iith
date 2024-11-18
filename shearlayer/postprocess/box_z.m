function [f_fil_z_indices] = box_z(f,Nz)
[nx, ny, nz] = size(f);

f_fil_z = f;

num_indices = floor(nz / (2 * Nz));
f_fil_z_indices = zeros(nx, ny, num_indices);

f_pad= zeros(nx, ny, nz + 2 * Nz);
f_pad(:, :, Nz+1:Nz+nz) = f;                      %Copy interiors
for k = 1:Nz
    f_pad(:, :, Nz+1-k) = f_pad(:, :, Nz+1-k+nz); %left
    f_pad(:, :, Nz+nz+k) = f_pad(:, :, Nz+nz+k-nz);%right
end

index=1;
%Apply box filter
for k = Nz+1:Nz+nz
    if mod(k-(Nz+1), 2*Nz) == 0  %%Calculating only at 2*Nz th point
       f_fil_z(:, :, k-Nz) = mean(f_pad(:, :, k-Nz:k+Nz),3);
       f_fil_z_indices(:, :, index) = f_fil_z(:, :, k-Nz);
       index=index+1;
    end
end

end
