function [dxidxij]=get_metrics(x,y,z,dxi,deta,dzeta)
[nx,ny,nz] = size(x);

disp('Finding metrics')
% Derivatives with respect to xi, eta, zeta
dxdxiij = zeros(nx, ny, nz, 9);
Inv_Jac = zeros(nx, ny, nz);

% Compute the derivatives
dxdxiij(:,:,:,1) = ddxi_compact10_gen(x,dxi);    % dxdxi
dxdxiij(:,:,:,2) = ddeta_compact10_gen(x,deta);  % dxdeta
dxdxiij(:,:,:,3) = ddzeta_compact10_gen(x,dzeta);% dxdzeta
dxdxiij(:,:,:,4) = ddxi_compact10_gen(y,dxi);    % dydxi
dxdxiij(:,:,:,5) = ddeta_compact10_gen(y,deta);  % dydeta
dxdxiij(:,:,:,6) = ddzeta_compact10_gen(y,dzeta);% dydzeta
dxdxiij(:,:,:,7) = ddxi_compact10_gen(z,dxi);    % dzdxi
dxdxiij(:,:,:,8) = ddeta_compact10_gen(z,deta);  % dzdeta
dxdxiij(:,:,:,9) = ddzeta_compact10_gen(z,dzeta);% dzdzeta

% Calculate the inverse Jacobian determinant
Inv_Jac = dxdxiij(:,:,:,1) .* (dxdxiij(:,:,:,5) .* dxdxiij(:,:,:,9) - dxdxiij(:,:,:,8) .* dxdxiij(:,:,:,6)) + ...
          dxdxiij(:,:,:,2) .* (dxdxiij(:,:,:,6) .* dxdxiij(:,:,:,7) - dxdxiij(:,:,:,9) .* dxdxiij(:,:,:,4)) + ...
          dxdxiij(:,:,:,3) .* (dxdxiij(:,:,:,4) .* dxdxiij(:,:,:,8) - dxdxiij(:,:,:,7) .* dxdxiij(:,:,:,5));

eps = 10^-32;

% Calculate the elements of the inverse Jacobian matrix
dxidxij = zeros(nx, ny, nz, 9);

dxidxij(:,:,:,1) =  (dxdxiij(:,:,:,5) .* dxdxiij(:,:,:,9) - dxdxiij(:,:,:,6) .* dxdxiij(:,:,:,8)) ./ (Inv_Jac + eps);
dxidxij(:,:,:,2) = -(dxdxiij(:,:,:,2) .* dxdxiij(:,:,:,9) - dxdxiij(:,:,:,3) .* dxdxiij(:,:,:,8)) ./ (Inv_Jac + eps);
dxidxij(:,:,:,3) =  (dxdxiij(:,:,:,2) .* dxdxiij(:,:,:,6) - dxdxiij(:,:,:,3) .* dxdxiij(:,:,:,5)) ./ (Inv_Jac + eps);
dxidxij(:,:,:,4) = -(dxdxiij(:,:,:,4) .* dxdxiij(:,:,:,9) - dxdxiij(:,:,:,6) .* dxdxiij(:,:,:,7)) ./ (Inv_Jac + eps);
dxidxij(:,:,:,5) =  (dxdxiij(:,:,:,1) .* dxdxiij(:,:,:,9) - dxdxiij(:,:,:,3) .* dxdxiij(:,:,:,7)) ./ (Inv_Jac + eps);
dxidxij(:,:,:,6) = -(dxdxiij(:,:,:,1) .* dxdxiij(:,:,:,6) - dxdxiij(:,:,:,3) .* dxdxiij(:,:,:,4)) ./ (Inv_Jac + eps);
dxidxij(:,:,:,7) =  (dxdxiij(:,:,:,4) .* dxdxiij(:,:,:,8) - dxdxiij(:,:,:,5) .* dxdxiij(:,:,:,7)) ./ (Inv_Jac + eps);
dxidxij(:,:,:,8) = -(dxdxiij(:,:,:,1) .* dxdxiij(:,:,:,8) - dxdxiij(:,:,:,2) .* dxdxiij(:,:,:,7)) ./ (Inv_Jac + eps);
dxidxij(:,:,:,9) =  (dxdxiij(:,:,:,1) .* dxdxiij(:,:,:,5) - dxdxiij(:,:,:,2) .* dxdxiij(:,:,:,4)) ./ (Inv_Jac + eps);

disp('Done finding metrics')
