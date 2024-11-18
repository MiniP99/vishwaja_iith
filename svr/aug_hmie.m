clear all,clc,clf
inputdir  = '/home/vsj/Codes/svr/';
outputdir  = '/home/vsj/Codes/svr/data_aug/';
mkdir(outputdir)
filename = 'TotalSample_IITH_240520_IITH_2nd_Eng.xlsx';
filepath = strcat(inputdir,filename);
data_old = readmatrix(filepath,'Sheet','2.Guess NextGen@Homo','Range','B7:D243');
data_new = readmatrix(filepath,'Sheet','2.Guess NextGen@Homo','Range','F7:H237');
%Training data
X_tra = [data_old(:,1), data_old(:,2)];
Y_tra = data_old(:,3);
Y_train_indices = find(Y_tra < 300);
X_train = [data_old(Y_train_indices,1), data_old(Y_train_indices,2)];
Y_train = data_old(Y_train_indices,3);

%Testing data
X_tes = [data_new(:,1), data_new(:,2)];
Y_tes = data_new(:,3);
Y_tes_indices = find(Y_tes < 300);
X_test  = [data_new(Y_tes_indices,1), data_new(Y_tes_indices,2)];
Y_test  = data_new(Y_tes_indices,3);
%%%% Train SVR model
diary(strcat(outputdir,'output_poly_train.txt'));
svrModel = fitrsvm(X_train, Y_train, 'KernelFunction', 'polynomial', 'KernelScale',1.4123, ...
           'Standardize',true,'BoxConstraint',127.11,'Epsilon',4.156, 'PolynomialOrder',3); %Train
conv = svrModel.ConvergenceInfo.Converged;
iter = svrModel.NumIterations;
disp(svrModel)
if conv==1
   disp([strcat('SVR Model converged in ',sprintf(' %d ', iter), ' iterations')]);
else
   disp([strcat('SVR Model not converged')]);
end
disp(svrModel.KernelParameters)
disp(strcat('Box Constraint:',num2str(svrModel.BoxConstraints(1,1))))
disp(strcat('Epsilon:',num2str(svrModel.Epsilon)))
%%% Predict on test set
Y_pred = predict(svrModel, X_train);
%%% Evaluate model performance
err = Y_pred - Y_train;
max_error = max(err);
mean_square_error = mean(err.^2)
rms_error         = sqrt(mean_square_error)
disp('=============================================================================')
diary off;


%%%% Step-2: Augmenting the data %%%%
x = X_train; y = Y_train;
rpm_min = min(x(:, 1));     rpm_max = max(x(:, 1));
torque_min = min(x(:, 2));  torque_max = max(x(:, 2));
n_augx = 64;  n_augy = 64;
rpm_new = linspace(rpm_min, rpm_max, n_augx);
torque_new = linspace(torque_min, torque_max, n_augy);
%for i=1:n_augx
%    for j=1:n_augy
%    rpm_grid(i,j) = rpm_new(i);
%    torque_grid(i,j) = torque_new(j);
%    end
%end
[torque_grid, rpm_grid] = meshgrid(torque_new, rpm_new);
x_new = [rpm_grid(:), torque_grid(:)];
y_new = predict(svrModel, x_new);
bsfc_grid = reshape(y_new, [n_augx,n_augy]);

figure,clf
n=256;
t1=linspace(min(min(bsfc_grid)),max(max(bsfc_grid)),n+1);
t1=t1(1:n);
contourf(rpm_grid, torque_grid, bsfc_grid,min(min( bsfc_grid)):5: max(max( bsfc_grid)),'ShowText','Off','LineColor','None');
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
set(ax, 'TickLabelInterpreter', 'latex');
xlab=xlabel('Speed$\,$(RPM)');set(xlab,'Fontsize',20,'interpreter','latex');
ylab=ylabel('Torque$\,$(Kgfm)');set(ylab,'Fontsize',20,'interpreter','latex');
tit=title('Predicted BSFC$\,$(g/KWhr)');set(tit,'Fontsize',20,'interpreter','latex');
colormap(jet(n));
xlim([min(min(rpm_grid)) max(max(rpm_grid))]); ylim([min(min(torque_grid)) max(max(torque_grid))]);
c=colorbar;c.FontSize=15;
set(c,'TickLabelInterpreter','latex')
screen2jpeg(strcat(outputdir,'cont_poly_aug.png'))

%%%% Step-3: Setting Threshold based on peak value%%%%
bsfc_min  = min(y_new);  bsfc_max  = max(y_new);
thres_bsfc = 20;
mask = find(y_new < bsfc_min + (thres_bsfc/100)*(bsfc_max-bsfc_min));
y_extracted = y_new(mask);
rpm_extracted = rpm_grid(mask);
torque_extracted = torque_grid(mask);
% Scatter plot of the extracted points
figure,clf
scatter(rpm_extracted, torque_extracted, 50, y_extracted, 'filled');box on
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
set(ax, 'TickLabelInterpreter', 'latex');
xlab=xlabel('Speed$\,$(RPM)');set(xlab,'Fontsize',20,'interpreter','latex');
ylab=ylabel('Torque$\,$(Kgfm)');set(ylab,'Fontsize',20,'interpreter','latex');
tit=title('Predicted BSFC$\,$(g/KWhr)');set(tit,'Fontsize',20,'interpreter','latex');
colormap(jet(n));caxis([bsfc_min bsfc_max])
xlim([min(min(rpm_grid)) max(max(rpm_grid))]); ylim([min(min(torque_grid)) max(max(torque_grid))]);
c=colorbar;c.FontSize=15;
set(c,'TickLabelInterpreter','latex')
screen2jpeg(strcat(outputdir,'sca_poly_aug_20p0.png'))


%% Step-4: Setting Threshold based on gradient%%%%
drpm = rpm_grid(2, 1) - rpm_grid(1, 1);  % Assuming uniform spacing in speed
dtorque = torque_grid(1, 2) - torque_grid(1, 1);  % Assuming uniform spacing in torque
[dBSFC_drpm] = ddxi_compact10_gen(bsfc_grid,drpm); %General boundary condition
[dBSFC_dtorque] = ddeta_compact10_gen(bsfc_grid,dtorque);  %General boundary condition
%[dBSFC_dtorque, dBSFC_drpm] = gradient(bsfc_grid, dtorque, drpm);  %in-buit function central difference
grad_mag = sqrt(dBSFC_drpm.^2 + dBSFC_dtorque.^2);
figure,clf;
%contourf(rpm_grid, torque_grid, grad_mag);
contourf(rpm_grid, torque_grid, grad_mag, min(min(grad_mag)):0.1: max(max(grad_mag)),'Linecolor','None');
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
set(ax, 'TickLabelInterpreter', 'latex');
xlab=xlabel('Speed$\,$(RPM)');set(xlab,'Fontsize',20,'interpreter','latex');
ylab=ylabel('Torque$\,$(Kgfm)');set(ylab,'Fontsize',20,'interpreter','latex');
tit=title('$|\nabla$BSFC$|$');set(tit,'Fontsize',20,'interpreter','latex');
colormap(jet(n));caxis([min(min(grad_mag)) max(max(grad_mag))])
xlim([min(min(rpm_grid)) max(max(rpm_grid))]); ylim([min(min(torque_grid)) max(max(torque_grid))]);
c=colorbar;c.FontSize=15;
set(c,'TickLabelInterpreter','latex')
screen2jpeg(strcat(outputdir,'grad_aug.png'))

% Set threshold on grad(bsfc)
grad_bsfc = grad_mag(:);
gmin = min(grad_bsfc); gmax = max(grad_bsfc);
thres_grad = 25;
grad_mask = find(grad_bsfc > gmin + (thres_grad/100)*(gmax-gmin));
grad_ext = grad_bsfc(grad_mask);
rpm_ext = rpm_grid(grad_mask);
torque_ext = torque_grid(grad_mask);
y_ext = bsfc_grid(grad_mask);
% Scatter plot of the extracted points
figure,clf
scatter(rpm_ext, torque_ext, 100, grad_ext, 'filled');box on;
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
set(ax, 'TickLabelInterpreter', 'latex');
xlab=xlabel('Speed$\,$(RPM)');set(xlab,'Fontsize',20,'interpreter','latex');
ylab=ylabel('Torque$\,$(Kgfm)');set(ylab,'Fontsize',20,'interpreter','latex');
tit=title('$|\nabla$BSFC$|$');set(tit,'Fontsize',20,'interpreter','latex');
colormap(jet(n));caxis([gmin gmax])
xlim([min(min(rpm_grid)) max(max(rpm_grid))]); ylim([min(min(torque_grid)) max(max(torque_grid))]);
c=colorbar;c.FontSize=15;
set(c,'TickLabelInterpreter','latex')
screen2jpeg(strcat(outputdir,'grad_bsfc_25p0.png'))

% Write data to a file
grad_extracted = grad_mag(mask); 
data_ext = [rpm_extracted, torque_extracted, y_extracted, grad_extracted];
header = {'Speed', 'TORQUE', 'BSFC', 'grad(BSFC)'};
fid=fopen(strcat(outputdir,'bsfc_20p0.txt'),'w');
fprintf(fid,'   %s     %s     %s     %s \r\n',header{:});
for k=1:size(data_ext,1)
 fprintf(fid, '    %4.2f',data_ext(k,:));
 fprintf(fid, '\r\n');
end

% Write data to a file
data_grad_ext = [rpm_ext, torque_ext, y_ext, grad_ext];
header = {'Speed', 'TORQUE', 'BSFC', 'grad(BSFC)'};
fid=fopen(strcat(outputdir,'grad_bsfc_25p0.txt'),'w');
fprintf(fid,'   %s     %s     %s     %s \r\n',header{:});
for k=1:size(data_grad_ext,1)
 fprintf(fid, '    %4.2f',data_grad_ext(k,:));
 fprintf(fid, '\r\n');
end

%Total data to a file
data_tot_ext = [rpm_grid(:), torque_grid(:), bsfc_grid(:), grad_mag(:)];
header = {'Speed', 'TORQUE', 'BSFC', 'grad(BSFC)'};
fid=fopen(strcat(outputdir,'aug_data_total.txt'),'w');
fprintf(fid,'   %s     %s     %s     %s \r\n',header{:});
for k=1:size(data_tot_ext,1)
 fprintf(fid, '    %4.2f',data_tot_ext(k,:));
 fprintf(fid, '\r\n');
end





%% Contour plot of the derivative with respect to speed
%figure,clf;
%contourf(rpm_grid, torque_grid, dBSFC_drpm);
%ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
%set(ax, 'TickLabelInterpreter', 'latex');
%xlab=xlabel('Speed$\,$(RPM)');set(xlab,'Fontsize',20,'interpreter','latex');
%ylab=ylabel('Torque$\,$(Kgfm)');set(ylab,'Fontsize',20,'interpreter','latex');
%tit=title('Derivative of BSFC w.r.t Speed');set(tit,'Fontsize',20,'interpreter','latex');
%colormap(jet(n));
%c=colorbar;c.FontSize=15;
%set(c,'TickLabelInterpreter','latex')
%screen2jpeg(strcat(outputdir,'rpm_poly_aug.png'))
%% Contour plot of the derivative with respect to torque
%figure,clf;
%contourf(rpm_grid, torque_grid, dBSFC_dtorque);
%ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
%set(ax, 'TickLabelInterpreter', 'latex');
%xlab=xlabel('Speed$\,$(RPM)');set(xlab,'Fontsize',20,'interpreter','latex');
%ylab=ylabel('Torque$\,$(Kgfm)');set(ylab,'Fontsize',20,'interpreter','latex');
%tit=title('Derivative of BSFC w.r.t Torque');set(tit,'Fontsize',20,'interpreter','latex');
%colormap(jet(n));
%c=colorbar;c.FontSize=15;
%set(c,'TickLabelInterpreter','latex')
%screen2jpeg(strcat(outputdir,'tor_poly_aug.png'))
% Contour plot of the magnitude of grad
