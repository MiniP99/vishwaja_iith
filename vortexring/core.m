function [eigenvalues]=core(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
[n_x,n_y]=size(dudx);

% Initialize a 3D array for the velocity gradient tensor
L = zeros(3, 3, n_x, n_y);

% Fill in the velocity gradient tensor for each point (example values)
for i = 1:n_x
    for j = 1:n_y
            % Define the velocity gradient tensor at (i, j, k)
            L(:,:,i,j) = [dudx(i,j), dudy(i,j), dudz(i,j); dvdx(i,j), dvdy(i,j), dvdz(i,j); ...
                            dwdx(i,j), dwdy(i,j), dwdz(i,j)]; 
    end
end

% Initialize arrays to store eigenvalues
eigenvalues = zeros(3, n_x, n_y);

% Compute eigenvalues for each tensor
for i = 1:n_x
    for j = 1:n_y
            % Compute eigenvalues of the velocity gradient tensor at (i, j, k)
            eig_values = eig(L(:,:,i,j));
            eigenvalues(:,i,j) = eig_values; % Store eigenvalues       
    end
end

%% Compute the imaginary parts of all eigenvalues
%imag_parts = imag(eigenvalues);
%
%% Initialize a matrix to store the indices where imaginary parts exist
%imag_indices = [];
%
%% Find the indices where the imaginary parts are non-zero
%for k = 1:3
%    [x_idx, y_idx] = find(squeeze(imag_parts(k, :, :)));
%    if ~isempty(x_idx)
%       % Concatenate the eigenvalue index, x index, and y index
%       imag_indices = [imag_indices; k * ones(size(x_idx)), x_idx, y_idx];
%    end
%end


% Compute the absolute imaginary parts of all eigenvalues
abs_imag_parts = abs(imag(eigenvalues));

% Initialize matrix to store the indices of the maximum imaginary values
max_imag_indices = zeros(3, 2);

% Find the indices of the maximum imaginary value for each eigenvalue
for k = 1:3
    [~, idx] = max(abs_imag_parts(k, 3:end, 3:end), [], 'all', 'linear');
    [x_idx, y_idx] = ind2sub([n_x, n_y], idx);
    max_imag_indices(k, :) = [x_idx, y_idx];
end

