clear all,clc,clf
inputdir  = '/home/vsj/Codes/svr/';
outputdir  = '/home/vsj/Codes/svr/fig/';
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

% To optimize Hyperparameters
diary('output_poly_train_2000.txt');
% Define the SVR model with a fixed Gaussian kernel and customize hyperparameter optimization
svrModel_opt = fitrsvm(X_train, Y_train, ...
    'KernelFunction', 'Polynomial', ... % Fixed kernel function
    'OptimizeHyperparameters', 'auto', ... % Hyperparameters to optimize
    'HyperparameterOptimizationOptions', struct(...
        'AcquisitionFunctionName', 'expected-improvement-plus', ... % Optimization strategy
        'MaxObjectiveEvaluations', 2000, ... % Adjust this number as needed
        'KFold', 5)); % Number of folds for cross-validation (optional)

conv = svrModel_opt.ConvergenceInfo.Converged;
iter = svrModel_opt.NumIterations;
disp([strcat('SVR Model converged in ',sprintf(' %d ', iter), ' iterations')]);
disp(svrModel_opt.KernelParameters)
disp(svrModel_opt)
%%% Predict on test set
Y_pred = predict(svrModel_opt, X_train);
%%% Evaluate model performance
err = Y_pred - Y_train;
mean_square_error = mean(err.^2)
rms_error         = sqrt(mean_square_error)
disp('=============================================================================')
diary off;

figure,clf;
%n=size(X_train,1);
%plot(1:n, Y_train', 'sk','MarkerFacecolor','k');hold on
%plot(1:n, Y_pred','+r','LineWidth',2);
plot(Y_train, Y_pred,'or','LineWidth',2,'HandleVisibility','off');hold on
plot(Y_train,Y_train,'--k','Linewidth',1,'DisplayName','Slope=1')
ax=gca;ax.YAxis.FontSize = 18;ax.XAxis.FontSize = 18;
%leg=legend('Testing Data Actual', 'SVR Prediction');
leg.ItemTokenSize = [3 3];
legend show;legend('Location','north','box','off','FontSize',17,'Numcolumns',3)
title(['Polynomial : RMS Error = ', sprintf('%2.2f',rms_error)],'Fontsize',15);
%xlab=xlabel('Test data points');set(xlab,'Fontsize',15);
xlab=xlabel('Actual BSFC (g/KWhr)');set(xlab,'Fontsize',20);
ylab=ylabel('Predicted BSFC (g/KWhr)');set(ylab,'Fontsize',20);
screen2jpeg(strcat(outputdir,'svr_train_poly.png'))

%figure,clf
%speed    = X_test(:,1);
%torque = X_test(:,2);
%bsfc   = Y_test;
%[X,Y] = meshgrid(speed, torque);
%% Interpolate BSFC values onto the grid
%Z = griddata(speed, torque, bsfc, X, Y,'linear');
%% Create the contour plot
%n=256;
%t1=linspace(min(min(Z)),max(max(Z)),n+1);
%t1=t1(1:n);
%%contourf(X, Y, Z,'LineColor','None','Levellist',t1);hold on
%contourf(X, Y, Z,'LineColor','None');hold on
%colormap(jet);colorbar
%xlim([min(speed) max(speed)]);ylim([min(torque) max(torque)]);
%%caxis([min(bsfc) max(bsfc)]);
%caxis([ 218.8584 678.1259]);
%scatter3(speed,torque,bsfc,20,'MarkerEdgeColor','None')
%colormap(jet);colorbar
%xlab=xlabel('Speed (RPM)');set(xlab,'Fontsize',15);
%ylab=ylabel('Torque (Kgfm)');set(ylab,'Fontsize',15);
%title('Actual New Engine BSFC (g/KWhr)','Fontsize',15);
%screen2jpeg(strcat(outputdir,'Test_new_data.png'))

