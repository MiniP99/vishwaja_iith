clear,clf
inputdir='/home/vsj/Codes/plot/';
outputdir='/home/vsj/Codes/plot/myplots/';
outputdir1='/home/vsj/Codes/plot/myplots/runs180_sta_sta_mgm/';
dir = '/home/vsj/Codes/runs180/';
mkdir(outputdir1)
%%%------------------------ subplot of R11,R12,R22,prod-diss of different models on single plot-------------------------
%caseid=["mgm_dyn","mgm_sta","mgm_sta_c1"];caseid_name=["mgm-dynamic","mgm-static","mgm-static(c=1)"];
caseid=["mgm_dyn","mgm_sta_c1"];caseid_name=["mgm-dynamic","mgm-static(c=1)"];
c={'r','b','g'};
%Mc = [0.2,0.4,0.8,1.2,1.6,2.0]; c=sqrt(1.4);du = 2*Mc*c;
Mc = [2.0]; c=sqrt(1.4);du = 2*Mc*c;
%Ly_mat = [200,200,100,100,80,80];ny = 240;         %%%-------------change  
Ly_mat = [80];ny = 240;         %%%-------------change  
col={'r','b','g','k','m','c'};
for i=1:size(Mc,2)
figure(i+3),clf 
for j=1:size(caseid,2)
Ly=Ly_mat(i);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 
fname1 = strcat(dir,caseid(j),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(i)),'_',caseid(j),'.dat');
f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
fname2 = strcat(dir,caseid(j),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(j),'.dat');
%f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
f2 = importdata(fname2); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);
sub_data= [f1(:,5),f1(:,6),f1(:,8)];
sub_txt = ["\sqrt{|R_{11}|}/{\Delta{u}}","\sqrt{|R_{12}|}/{\Delta{u}}","\sqrt{|R_{22}|}/{\Delta{u}}"];

for k=1:3
subplot(2,2,k)
plot(y/del_tavg,sub_data(:,k),'color',col{j},'LineStyle','-','LineWidth',1,'DisplayName',caseid_name(j));hold on
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
xlabel('y/{\delta_\theta(\tau)}');
%xlim([-10 10]);
h=ylabel(strcat('\boldmath$',sub_txt(k),'$'));set(h,'Interpreter','latex');
end

subplot(2,2,4)
plot(y/del_tavg,f1(:,11),'color',col{j},'LineStyle','-','LineWidth',1,'DisplayName',caseid_name(j));hold on;
xlim([-7 7]);
plot(y/del_tavg,-f1(:,12),'color',col{j},'LineStyle','-.','LineWidth',1,'HandleVisibility','off');hold on;
xlabel('y/{\delta_\theta(\tau)}');
xlim([-6 6]);
end
ax=gca;  ax.YRuler.Exponent = 0;
set(findall(gcf,'-property','FontSize'),'FontSize',6)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
txt1 = ' Production(P)\rightarrow';
t = text(-6,0.0004,txt1,'Fontsize',4,'FontWeight','bold');
txt2 = 'Dissipation(\epsilon)\rightarrow';
t = text(-6,-0.0003,txt2,'Fontsize',4,'FontWeight','bold');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
legend show
leg=legend('Location','Southeast','Fontsize',4,'box','off');
leg.ItemTokenSize = [3 3];
screen2jpeg(strcat(outputdir1,'compare_models_',sprintf('%1.1f',Mc(i)),'_.png'));
end

%%----------------------------------------- Single model at different mc subplot--------------------
for i=1:size(caseid,2)
figure(i+10),clf
for j=1:size(Mc,2)
Ly=Ly_mat(j);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 
fname1 = strcat(dir,caseid(i),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(j)),'_',caseid(i),'.dat');
f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
fname2 = strcat(dir,caseid(i),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(i),'.dat');
%f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
f2 = importdata(fname2); 
idx=find(abs(f2(:,1)-Mc(j))<0.001);del_tavg=f2(idx,2);

sub_data= [f1(:,5),f1(:,6),f1(:,8)];
sub_txt = ["\sqrt{|R_{11}|}/{\Delta{u}}","\sqrt{|R_{12}|}/{\Delta{u}}","\sqrt{|R_{22}|}/{\Delta{u}}"];

for k=1:3
subplot(2,2,k)
plot(y/del_tavg,sub_data(:,k),'color',col{j},'LineStyle','-','LineWidth',1,'DisplayName',strcat('M_c=',sprintf('%1.1f',Mc(j))));hold on;
xlabel('y/{\delta_\theta(\tau)}');
%xlim([-10 10]);
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
h=ylabel(strcat('\boldmath$',sub_txt(k),'$'));set(h,'Interpreter','latex');
end

subplot(2,2,4)
plot(y/del_tavg,f1(:,11),'color',col{j},'LineStyle','-','LineWidth',1,'DisplayName',strcat('M_c=',sprintf('%1.1f',Mc(j))));hold on;
plot(y/del_tavg,-f1(:,12),'color',col{j},'LineStyle','-.','LineWidth',1,'HandleVisibility','off');hold on;
xlabel('y/{\delta_\theta(\tau)}');
xlim([-6 6]);
end
ax=gca;  ax.YRuler.Exponent = 0;
set(findall(gcf,'-property','FontSize'),'FontSize',6)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
txt1 = ' Production(P)\rightarrow';
t = text(-6,0.0004,txt1,'Fontsize',4,'FontWeight','bold');
txt2 = 'Dissipation(\epsilon)\rightarrow';
t = text(-6,-0.0004,txt2,'Fontsize',4,'FontWeight','bold');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
legend show
leg=legend('Location','Northeast','Fontsize',4,'box','off');
leg.ItemTokenSize = [3 3];
screen2jpeg(strcat(outputdir1,'compare_mc_',caseid(i),'_.png'));
end
