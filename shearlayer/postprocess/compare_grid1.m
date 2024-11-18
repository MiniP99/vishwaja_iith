clear,clf
inputdir='/scratch/peddamma.iith/';
outputdir='/home/peddamma.iith/Codes/plot/myplots/';
outputdir1='/home/peddamma.iith/Codes/plot/myplots/grid_len/';
mkdir(outputdir1);
sub_name="128";
cases=["runs128","runs128_len"];
grid1=["80x80x40","80x120x40"];
caseid="mgm_sta_c1";
all_marks = {'>','o','^','s','d','x','v','+','<','p','h'};
Mc = [0.2,0.4,0.8,1.2,1.6,2.0]; c=sqrt(1.4);du = 2*Mc*c;
Ly_mat = [80,120];ny = [180,256];         %%%-------------change  
delta_0=1;
del_the = 0.018;


col={'r','b','g','k','m','c'};
for i=6%:size(Mc,2)
figure(i+12),clf 
for j=1:size(cases,2)
Ly=Ly_mat(j);dy = Ly/ny(j); y = linspace(-Ly/2,Ly/2,ny(j))'; 
fname1 = strcat(inputdir,cases(j),'/',caseid,'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(i)),'_',caseid,'.dat');
f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
fname2 = strcat(inputdir,cases(j),'/',caseid,'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid,'.dat');
%f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
f2 = importdata(fname2); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);

subplot(2,2,1)
plot(y/del_tavg,f1(:,5),'color',col{j},'Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'DisplayName',grid1(j));hold on
xlabel('y/{\delta_\theta(\tau)}');xlim([-11 11]);
h=ylabel('\boldmath$ \sqrt{|R_{11}|}/{\Delta{u}}$');set(h,'Interpreter','latex');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')

subplot(2,2,2)
plot(y/del_tavg,f1(:,6),'color',col{j},'Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'DisplayName',grid1(j));hold on
xlabel('y/{\delta_\theta(\tau)}');xlim([-11 11]);
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
h=ylabel('\boldmath$ \sqrt{|R_{12}|}/{\Delta{u}}$');set(h,'Interpreter','latex');

subplot(2,2,3)
plot(y(14:end-14)/del_tavg,f1(14:end-14,8),'color',col{j},'Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'DisplayName',grid1(j));hold on
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
xlabel('y/{\delta_\theta(\tau)}');xlim([-11 11]);
h=ylabel('\boldmath$ \sqrt{|R_{22}|}/{\Delta{u}}$');set(h,'Interpreter','latex');

subplot(2,2,4)
plot(y/del_tavg,f1(:,11),'color',col{j},'LineStyle','-','Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'DisplayName',grid1(j));hold on;xlim([-7 7]);
plot(y/del_tavg,-f1(:,12),'color',col{j},'LineStyle','-.','Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'HandleVisibility','off');hold on;
xlabel('y/{\delta_\theta(\tau)}');xlim([-6 6]);
end
ax=gca;  ax.YRuler.Exponent = 0;
set(findall(gcf,'-property','FontSize'),'FontSize',6)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',2)
txt1 = ' Production(P)\rightarrow';
t = text(-5.5,0.0004,txt1,'Fontsize',4,'FontWeight','bold');
txt2 = 'Dissipation(\epsilon)\rightarrow';
t = text(-5,-0.0004,txt2,'Fontsize',4,'FontWeight','bold');
%xlim([-10 10]);
%xlabel('y/{\delta_\theta(\tau)}','FontWeight','bold');
%h=ylabel('\boldmath$ \sqrt{|R_{12}|}/{\Delta{u}}$');set(h,'Interpreter','latex');
%ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
legend show
leg=legend('Location','Southeast','Fontsize',4,'box','off');
leg.ItemTokenSize = [3 3];
screen2jpeg(strcat(outputdir1,'subplot_',sub_name,'_.png'));
end





%%-----------------momentum thickness-------------------
figure(10),clf
mc=["20"];
Mc=[2.0];
du=2*Mc*sqrt(1.4);
t_i=[180,160];t_f=[315,300];  %AMD
caseid="mgm_sta_c1";
for i=1:size(cases,2)
fname = strcat(inputdir,cases(i),'/postprocess/','delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
f=readmatrix(fname);
idx1=find(f(:,1)==t_i(i));idx2=find(f(:,1)==t_f(i));
tau=(f(1:idx2,1)*du)/delta_0;
del=f(1:idx2,2)/delta_0;
plot(tau,del,'Color',col{i},'LineWidth',1,'Marker',all_marks{i},'DisplayName',grid1(i));
hold on

tau=(f(idx1:idx2,1)*du)/delta_0;
del=f(idx1:idx2,2)/delta_0;
eq=polyfit(tau,del,1);
delta1 = polyval(eq,tau);
fittxt=strcat('Linearfit:y=',sprintf('%1.3f',eq(1,1)),'x+(',sprintf('%1.3f',eq(1,2)),')');
plot(tau,delta1,'--k','LineWidth',2,'DisplayName',fittxt);
hold on
end
legend show
title('momentum thickness vs time'), 
xlabel('(t*du)/\delta_{\theta}_0');ylabel('\delta_{\theta}/ \delta_{\theta}_0');
legend('Location','northwest','box','off');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
grid on;hAx=gca; set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir1,'mom_thick__',sub_name,'.png'));

