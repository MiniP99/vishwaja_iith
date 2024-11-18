clear all,clc,clf
inputdir  = '/home/vsj/Codes/svr/';
outputdir  = '/home/vsj/Codes/svr/fig/hmie/';
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
diary('output_poly_train.txt');
%svrModel = fitrsvm(X_train, Y_train, 'KernelFunction', 'gaussian', 'KernelScale',1.5567, ...
%           'Standardize',true,'BoxConstraint',994.81,'Epsilon',1.1794 ); %Train
svrModel = fitrsvm(X_train, Y_train, 'KernelFunction', 'polynomial', 'KernelScale',1.4123, ...
           'Standardize',true,'BoxConstraint',127.11,'Epsilon',4.156, 'PolynomialOrder',3); %Train
%svrModel = fitrsvm(X_test, Y_test, 'KernelFunction', 'gaussian', 'KernelScale',0.87807, ...
%           'Standardize',true,'BoxConstraint',999.9,'Epsilon',0.025113 ); %Train
%svrModel = fitrsvm(X_test, Y_test, 'KernelFunction', 'polynomial', 'KernelScale',0.12441, ...
%           'Standardize',true,'BoxConstraint',0.006108,'Epsilon', 19.8, 'PolynomialOrder',3); %Train
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

figure,clf;
plot(Y_train, Y_pred,'or','LineWidth',2,'HandleVisibility','off');hold on
plot(Y_train,Y_train,'--k','Linewidth',1,'DisplayName','Slope=1')
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
set(ax, 'TickLabelInterpreter', 'latex');
xlab=xlabel('Actual BSFC$\,$(g/KWhr)');set(xlab,'Fontsize',18,'interpreter','latex');
ylab=ylabel('Predicted BSFC$\,$(g/KWhr)');set(ylab,'Fontsize',18,'interpreter','latex');
tit=title(['Polynomial : RMS Error = ', sprintf('%2.2f',rms_error)]);set(tit,'Fontsize',18,'interpreter','latex');
legend show;leg=legend('Location','northwest','box','off','FontSize',18,'Numcolumns',1,'Interpreter','latex');
screen2jpeg(strcat(outputdir,'poly_train.png'))

isSupportVector = svrModel.IsSupportVector;
supportVectors  = Y_train(isSupportVector, :);
notsupportVectors = Y_train(~isSupportVector, :);
n1=size(supportVectors,1);n2=size(notsupportVectors,1);
figure,clf;
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
%%% Predict on test set
Y_train_pred = resubPredict(svrModel);
L = resubLoss(svrModel)
txt1=strcat('Support Vectors:',sprintf('%d',n1));
txt2=strcat('Not Support Vectors:',sprintf('%d',n2));
plot(Y_train(isSupportVector, :), Y_pred(isSupportVector, :),'LineWidth',2,'DisplayName',txt1,'Marker','o','Color',myg,'Linestyle','None');hold on
plot(Y_train(~isSupportVector, :), Y_pred(~isSupportVector, :),'LineWidth',2,'DisplayName',txt2,'Marker','o','Color','r','Linestyle','None');hold on
plot(Y_train,Y_train,'--k','Linewidth',1,'DisplayName','Slope=1')
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
set(ax, 'TickLabelInterpreter', 'latex');
legend show;leg=legend('Location','northwest','box','off','FontSize',18,'Numcolumns',1,'Interpreter','latex');
title(['Polynomial : RMS Error = ', sprintf('%2.2f',rms_error)],'Fontsize',15);
xlab=xlabel('Actual BSFC$\,$(g/KWhr)');set(xlab,'Fontsize',18,'interpreter','latex');
ylab=ylabel('Predicted BSFC$\,$(g/KWhr)');set(ylab,'Fontsize',18,'interpreter','latex');
tit=title(['Polynomial : RMS Error = ', sprintf('%2.2f',rms_error)]);set(tit,'Fontsize',18,'interpreter','latex');
screen2jpeg(strcat(outputdir,'poly_train_support.png'))
