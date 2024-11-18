%%========================= Plot of SGS cofficients =============%% 
clear,clf
inputdir='/home/vsj/Codes/check_mgm_dyn/';
outputdir='/home/vsj/Codes/plot/myplots/';
outputdir1='/home/vsj/Codes/plot/myplots/grid_clocal/';
mkdir(outputdir1);
cases=["runs128","runs180","runs256"];
grid1=["128x180x64","180x240x84","256x360x128"];
%cases=["run20_25","run20","run20_75"];
%grid1=["L_{s}=2.5","L_{s}=5","L_{s}=7.5"];
caseid="mgm_dyn";
all_marks = {'>','o','^','s','d','x','v','+','<','p','h'};
Mc = [0.2,0.4,0.8,1.2,1.6,2.0]; c=sqrt(1.4);du = 2*Mc*c;
mc = ["02","04","08","12","16","20"];
Ly_mat = [200,200,100,100,80,80];
%ny = [180,240,360];         %%%-------------change  
ny = [180,240,360];         %%%-------------change  
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={myg,'r','b'};



%t_i=[2000,2700,2000];t_f=[4000,5400,6800];    %%MGM
t_i=[0,2800,4500];t_f=[0,7300,11700];    %%MGM

%% ---------------- Plot of averaged profiles along y at different Mc --------------------- %%
for i=6%:size(mc,2)
for k=2:size(cases,2)
step = t_i(k):100:t_f(k);
Ly=Ly_mat(i);dy = Ly/ny(k); y = linspace(-Ly/2,Ly/2,ny(k))'; 
c_local=zeros(ny(k),size(step,2));c_local_qsgs=zeros(ny(k),size(step,2));

for j= 1: size(step,2)
%fname=strcat(inputdir,cases(k),'/','run',mc(i),'/',caseid,'/cmodel_',sprintf('%05d',step(j)),'.dat');
fname=strcat(inputdir,cases(k),'/cmodel_',sprintf('%05d',step(j)),'.dat');
f=importdata(fname);
c_local(:,j)=f(:,2);c_local_qsgs(:,j)=f(:,3);
%c_local(:,j)=f(:,2)*(1.5^(4/3));c_local_qsgs(:,j)=f(:,3)*(1.5^(4/3));  %Sigma,mgm
%c_local(:,j)=f(:,2)*(1.5*1.5);c_local_qsgs(:,j)=f(:,3)*(1.5*1.5);   %AMD
end
size(c_local)
fname2 = strcat(inputdir,cases(k),'/',caseid,'_data','/','delta_tavg_',caseid,'.dat');
%f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
f2 = importdata(fname2); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);y_c=f2(idx,7);
c_local_tavg=sum(c_local,2)/size(step,2);c_local_qsgs_tavg=sum(c_local_qsgs,2)/size(step,2);
%plot(f(8:end-8,1)./del_tavg,c_local_qsgs_tavg(8:end-8),'Marker',all_marks{i},'Markerindices',8:15:length(y)-16,'color',col{i},'LineWidth',1.5,'DisplayName',strcat(sprintf('%1.1f',Mc(i))));hold on
plot(f(:,1)./del_tavg,c_local_qsgs_tavg,'Marker',all_marks{k},'Markerindices',1:15:length(y),'color',col{k},'LineWidth',1.8,'DisplayName',strcat(grid1(k)));hold on

txt1=strcat(sprintf('%1.4f',mean(c_local_qsgs_tavg)))
yline(mean(c_local_qsgs_tavg),'--r',txt1,'LineWidth',1,'HandleVisibility','off');hold on;
%c_local=zeros(ny(k),size(step,2));c_local_qsgs=zeros(ny(k),size(step,2));

end

end
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',20);
ylab=ylabel('\boldmath$C_{\epsilon}$','Interpreter','latex');set(ylab,'Fontsize',20);
xlim([-10 10])
%ylim([0 0.45])
legend show;legend('Location','Northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',6)
ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'c_eps_tavg_profiles_20',caseid,'.jpeg'));


