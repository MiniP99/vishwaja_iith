clear,clf
inputdir='/home/vsj/Codes/plot/';
outputdir='/home/vsj/Codes/plot/myplots/';
outputdir1='/home/vsj/Codes/plot/myplots/runs180/tke_plots/';
dir = '/home/vsj/Codes/runs180/';
dir1 = '/home/vsj/dns_data_csl/';
mkdir(outputdir1)

kris_Rij_int = readmatrix(strcat(outputdir,'kris_rij_int.csv'));
kris_p12     = readmatrix(strcat(outputdir,'kris_intp12.csv'));
kris_diss    = readmatrix(strcat(outputdir,'kris_intdiss.csv'));
kris_dissprod= readmatrix(strcat(outputdir,'kris_intdissprod.csv'));
kris_pi1122  = readmatrix(strcat(outputdir,'kris_pi1122.csv'));
kris_pi11p12 = readmatrix(strcat(outputdir,'kris_pi11p12.csv'));
kris_pi12p12 = readmatrix(strcat(outputdir,'kris_pi12p12.csv'));
kris_pi22p12 = readmatrix(strcat(outputdir,'kris_pi22p12.csv'));

myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={myg,'r','b','none','none'};
col1={myg,'r','b','k','k'};
all_marks = {'d','^','v','s','p'};
all_marks_1 = {'o','v','x','*','p','>','h','<'};
%%Plotting Integrated Rij/R11
figure(20),clf
caseid=["R11","R22","R12"];caseid_name=["R_{11}","R_{22}","R_{12}"];
for i=1:3
if(i==1)
j=1;
else 
j=i*(i+1)+1;
end
disname=strcat(caseid_name(i),'-DNS');
scatter(kris_Rij_int(j:j+5,1),kris_Rij_int(j:j+5,2),'k','Marker',all_marks_1{i+5},'MarkerFaceColor','k','Displayname',disname);hold on
end                % [0.7 0.7 0.7]

caseid = ["amd_dyn","sigma_dyn","mgm_dyn","lad","base"];
caseid_name = ["AMD-Dyn","Sigma-Dyn","MGM-Dyn","LAD","Base"];
c={'r','b','g'};
for i=1:size(caseid,2)
fname = strcat(dir,caseid(i),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(i),'.dat');
f = readmatrix(fname,'ConsecutiveDelimitersRule','join'); 
%1-Mc,2-del_avg,3-R11,4-R12,5-R13,6-R22,7-R23,8-R33,9-ke,10-P12,11-diss,12-pi11,13-pi12,14-pi22
plot(f(:,1),f(:,3)./f(1,3),'r','Marker',all_marks{i},'LineStyle','none','DisplayName',caseid_name{i});hold on
plot(f(:,1),f(:,4)./f(1,3),'g','Marker',all_marks{i},'LineStyle','none','HandleVisibility','off');hold on
plot(f(:,1),f(:,6)./f(1,3),'m','Marker',all_marks{i},'LineStyle','none','HandleVisibility','off');hold on
end
%ylim([-0.4 1.3])
ylim([-0.4 1.7])
xticks([0.2,0.4,0.8,1.2,1.6,2.0]);
ax=gca;ax.YAxis.FontSize = 13;ax.XAxis.FontSize = 13;
ax=gca;  ax.YRuler.Exponent = 0; grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
xlab=xlabel('\boldmath$M_c$','Interpreter','latex');set(xlab,'Fontsize',17);
ylab=ylabel('\boldmath$R_{ij}/R_{11}(0.2)$', 'Interpreter','latex');set(ylab,'Fontsize',17);
legend show;legend('Location','northeast','box','off','FontSize',10,'Numcolumns',2)
set(findall(gcf,'-property','Markersize'),'Markersize',7)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','FontName'),'FontName','Times');
screen2jpeg(strcat(outputdir1,'Rey_int.jpeg'));

%-----------------------------------------------------------------------------------------------------
%%%P_12,\epsilon,\epsilon/P12,\pi_{11}/P_{12},\pi_{22}/P_{12},\pi_{12}/P_{12},\pi_{22}/\pi{11}
figure(1),clf

plot(kris_diss(:,1),kris_diss(:,2),'Marker','o','MarkerFaceColor','k','MarkerEdgecolor','k','Displayname','DNS');hold on
%plot(kris_pi22p12(:,1),kris_pi22p12(:,2).*kris_p12(:,2),'Marker','o','MarkerFaceColor','k','MarkerEdgecolor','k','Displayname','DNS');hold on
for i=1:size(caseid,2)
fname = strcat(dir,caseid(i),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(i),'.dat');
f = readmatrix(fname,'ConsecutiveDelimitersRule','join'); 
%1-Mc,2-del_avg,3-R11,4-R12,5-R13,6-R22,7-R23,8-R33,9-ke,10-P12,11-diss,12-pi11,13-pi12,14-pi22
plot(f(:,1),f(:,11),'k','Marker',all_marks{i},'MarkerFaceColor',col{i},'MarkerEdgeColor',col1{i},'DisplayName',caseid_name(i));hold on
end
ylim([0.0005 0.0055])
legend show;legend('Location','Northeast','box','off','Fontsize',12,'Fontname','Timesnewroman');
xticks([0.2,0.4,0.8,1.2,1.6,2.0]);ax=gca;  ax.YRuler.Exponent = 0; 
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
ax=gca;ax.YAxis.FontSize = 13;ax.XAxis.FontSize = 13;
xlab=xlabel('\boldmath$M_c$','Interpreter','latex');set(xlab,'Fontsize',17);
ylab=ylabel('\boldmath$\epsilon$', 'Interpreter','latex');set(ylab,'Fontsize',17);
%P_{12},\epsilon,\epsilon/P_{12},\pi_{11}/P_{12},\pi_{22}/P_{12},\pi_{12}/P_{12},\pi_{22}/\pi_{11}
set(findall(gcf,'-property','Markersize'),'Markersize',8)
set(findall(gcf,'-property','FontName'),'FontName','Times');
set(findall(gcf,'-property','Linestyle'),'Linestyle','none')
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
screen2jpeg(strcat(outputdir1,'diss.jpeg'));
%P12,diss,dissprod,pi11p12,pi22p12,pi12p12,pi22pi11


%-----------------------------------------------------------------------------------------------------
figure(2),clf
kris_R11 = readmatrix(strcat(outputdir,'kris_R11_peak.csv'));
kris_R22 = readmatrix(strcat(outputdir,'kris_R22_peak.csv'));
kris_R12 = readmatrix(strcat(outputdir,'kris_R12_peak.csv'));


Mc = [0.2,0.4,0.8,1.2,1.6,2.0]; c=sqrt(1.4);
du = 2*Mc*c;du=du';
R_max(:,1) = Mc ; 


plot(kris_R12(:,1),kris_R12(:,2),'Marker','o','MarkerFaceColor','k','MarkerEdgecolor','k','Displayname','DNS');hold on

for i=1:size(caseid,2)-1
for j=1:size(Mc,2)
fname = strcat(dir,caseid(i),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(j)),'_',caseid(i),'.dat');
f = readmatrix(fname,'ConsecutiveDelimitersRule','join'); 
R_max(j,2)= max(abs(f(:,6)));
end
plot(R_max(:,1),R_max(:,2),'k','Marker',all_marks{i},'MarkerFaceColor',col{i},'MarkeredgeColor',col1{i},'DisplayName',caseid_name(i));
hold on
end

for i=5
for j=1:4
fname = strcat(dir,caseid(i),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(j)),'_',caseid(i),'.dat');
f = readmatrix(fname,'ConsecutiveDelimitersRule','join'); 
R_max(j,2)= max(abs(f(:,6)));   %%5-R11,8-R22,6-R12
end
plot(R_max(1:4,1),R_max(1:4,2),'k','Marker',all_marks{i},'MarkerFaceColor',col{i},'MarkeredgeColor',col1{i},'DisplayName',caseid_name(i));hold on
end

ylim([0.05 0.11]);
legend show;legend('Location','Northeast','box','off','Fontsize',12,'Fontname','Timesnewroman');
xticks([0.2,0.4,0.8,1.2,1.6,2.0]);
ax=gca;  ax.YRuler.Exponent = 0; grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
ax=gca;ax.YAxis.FontSize = 13;ax.XAxis.FontSize = 13;
xlab=xlabel('\boldmath$M_c$','Interpreter','latex');set(xlab,'Fontsize',17);
ylab=ylabel('\boldmath$\sqrt{|R_{12}|}/\Delta u$', 'Interpreter','latex');set(ylab,'Fontsize',17);
set(findall(gcf,'-property','Linestyle'),'Linestyle','none')
set(findall(gcf,'-property','Markersize'),'Markersize',8)
set(findall(gcf,'-property','FontName'),'FontName','Times');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
screen2jpeg(strcat(outputdir1,'R12_peak.jpeg'));
%-------------------------------------------------------------------------------------------
Mc = [0.2,0.4,0.8,1.2,1.6,2.0]; c=sqrt(1.4);du = 2*Mc*c;
Ly_mat = [200,200,100,100,80,80];ny = 240;         %%%-------------change  

%%  1-ub,    2-util,  3-vbar,  4-vtil,  5-R11,  6-R12, 7-R13, 8-R22, 9-R23, 10-R33
%%  11-prod,12-eps,   13-trans, 14-baro, 15-pre_dil,  16-mass_ru,  17-mass_rv

%-----------------------------------------------------------------------------------------------
%%% subplot of R11,R12,R22,prod-diss of different models on single plot-------------------------
%col={'r',myg,'b',myv,myy,'k'};
%for i=1:size(Mc,2)-2
%figure(i+2),clf 
%for j=1:size(caseid,2)
%Ly=Ly_mat(i);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 
%fname1 = strcat(dir,caseid(j),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(i)),'_',caseid(j),'.dat');
%f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
%fname2 = strcat(dir,caseid(j),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(j),'.dat');
%f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
%idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);
%sub_data= [f1(:,5),f1(:,6),f1(:,8)];
%sub_txt = ["\sqrt{|R_{11}|}/{\Delta{u}}","\sqrt{|R_{12}|}/{\Delta{u}}","\sqrt{|R_{22}|}/{\Delta{u}}"];
%
%for k=1:3
%subplot(2,2,k)
%plot(y/del_tavg,sub_data(:,k),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-','LineWidth',1,'DisplayName',caseid_name(j));hold on
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%xlabel('y/{\delta_\theta(\tau)}');xlim([-10 10]);
%h=ylabel(strcat('\boldmath$',sub_txt(k),'$'));set(h,'Interpreter','latex');
%set(findall(gcf,'-property','Markersize'),'Markersize',3)
%end
%
%subplot(2,2,4)
%plot(y/del_tavg,f1(:,11),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-','LineWidth',1,'DisplayName',caseid_name(j));hold on;xlim([-7 7]);
%plot(y/del_tavg,-f1(:,12),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-.','LineWidth',1,'HandleVisibility','off');hold on;
%xlabel('y/{\delta_\theta(\tau)}');xlim([-6 6]);
%end
%ax=gca;  ax.YRuler.Exponent = 0;
%set(findall(gcf,'-property','Markersize'),'Markersize',3)
%set(findall(gcf,'-property','FontSize'),'FontSize',6)
%set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%txt1 = ' Production(P)\rightarrow';
%t = text(-6,0.0004,txt1,'Fontsize',4,'FontWeight','bold');
%txt2 = 'Dissipation(\epsilon)\rightarrow';
%t = text(-6,-0.0003,txt2,'Fontsize',4,'FontWeight','bold');
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%legend show
%leg=legend('Location','Southeast','Fontsize',4,'box','off');
%leg.ItemTokenSize = [3 3];
%screen2jpeg(strcat(outputdir1,'compare_models_',sprintf('%1.1f',Mc(i)),'_.jpeg'));
%end

%%% If single plot needed----------------
col={myg,'r','b','none','none'};
col1={myg,'r','b','k','k'};
all_marks = {'d','^','v','s','p'};
caseid_name=["AMD","Sigma","MGM"];
for i=4:size(Mc,2)
figure(i+2),clf
Ly=Ly_mat(i);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)';
figure(i+2),clf
for j=1:size(caseid,2)-2
fname1 = strcat(dir,caseid(j),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(i)),'_',caseid(j),'.dat');
f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
fname2 = strcat(dir,caseid(j),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(j),'.dat');
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);
plot(y/del_tavg,f1(:,8),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-','LineWidth',2,'DisplayName',caseid_name(j));hold on
%%plot(y/del_tavg,-f1(:,12),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-.','LineWidth',1.5,'HandleVisibility','off');hold on
end  %%5-R11,8-R22,6-R12,11-prod,12-diss


ny_dns=1448;
dy_dns = Ly/ny_dns; y_dns = linspace(-Ly/2,Ly/2,ny_dns)'; 
fname1_dns = strcat(dir1,'dns_pp_data/','vel_rey_nor_',sprintf('%1.1f',Mc(i)),'_','dns.dat');
f1_dns = readmatrix(fname1_dns,'ConsecutiveDelimitersRule','join');  
fname2_dns = strcat(dir1,'dns_pp_data/','delta_tavg_','dns.dat');
f2_dns = readmatrix(fname2_dns,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2_dns(:,1)-Mc(i)) < 0.001); %Mc,delta_theta
del_tavg_dns=f2_dns(idx,2);
plot(y_dns/del_tavg_dns,f1_dns(:,5),'marker','o','color','k','LineWidth',1.5,'Linestyle','none','Displayname','DNS','markerindices',1:20:length(y_dns));hold on


ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
xlim([-10 10]);
%ylim([0 0.15]);
xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',20);
h=ylabel('\boldmath$\sqrt{|R_{22}|}/{\Delta{u}}$');set(h,'Interpreter','latex');set(h,'Fontsize',20);
legend show;legend('Location','northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
ax=gca;  ax.YRuler.Exponent = 0; 
set(findall(gcf,'-property','Markersize'),'Markersize',6)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','FontName'),'FontName','Times');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%%txt1 = ' Production(P)\rightarrow';
%%t = text(-4.5,0.001,txt1,'Fontsize',9,'FontWeight','bold');
%%txt2 = 'Dissipation(\epsilon)\rightarrow';
%%t = text(-4.5,-0.001,txt2,'Fontsize',9,'FontWeight','bold');
screen2jpeg(strcat(outputdir1,'R22_profile_',sprintf('%1.1f',Mc(i)),'_.jpeg'));
end
%%%%%-----------------------------------------------------------------------------------------------------






%% Single model at different mc subplot--------------------
%all_marks = {'o','v','d','^','s','>','x','v','>','<','p','h'};
%for i=1%:size(caseid,2)
%figure(i+10),clf
%for j=1:size(Mc,2)-2
%Ly=Ly_mat(j);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 
%fname1 = strcat(dir,caseid(i),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(j)),'_',caseid(i),'.dat');
%f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
%fname2 = strcat(dir,caseid(i),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(i),'.dat');
%f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
%idx=find(abs(f2(:,1)-Mc(j))<0.001);del_tavg=f2(idx,2);
%
%sub_data= [f1(:,5),f1(:,6),f1(:,8)];
%sub_txt = ["\sqrt{|R_{11}|}/{\Delta{u}}","\sqrt{|R_{12}|}/{\Delta{u}}","\sqrt{|R_{22}|}/{\Delta{u}}"];
%
%for k=1:3
%subplot(2,2,k)
%plot(y/del_tavg,sub_data(:,k),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-','LineWidth',1,'DisplayName',strcat('M_c=',sprintf('%1.1f',Mc(j))));hold on;
%xlabel('y/{\delta_\theta(\tau)}');xlim([-10 10]);
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%h=ylabel(strcat('\boldmath$',sub_txt(k),'$'));set(h,'Interpreter','latex');
%set(findall(gcf,'-property','Markersize'),'Markersize',3)
%end
%
%subplot(2,2,4)
%plot(y/del_tavg,f1(:,11),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-','LineWidth',1,'DisplayName',strcat('M_c=',sprintf('%1.1f',Mc(j))));hold on;
%plot(y/del_tavg,-f1(:,12),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-.','LineWidth',1,'HandleVisibility','off');hold on;
%xlabel('y/{\delta_\theta(\tau)}');xlim([-6 6]);
%end
%ax=gca;  ax.YRuler.Exponent = 0;
%set(findall(gcf,'-property','Markersize'),'Markersize',3)
%set(findall(gcf,'-property','FontSize'),'FontSize',6)
%set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%txt1 = ' Production(P)\rightarrow';
%t = text(-6,0.0004,txt1,'Fontsize',4,'FontWeight','bold');
%txt2 = 'Dissipation(\epsilon)\rightarrow';
%t = text(-6,-0.0004,txt2,'Fontsize',4,'FontWeight','bold');
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%legend show
%leg=legend('Location','Northeast','Fontsize',4,'box','off');
%leg.ItemTokenSize = [3 3];
%screen2jpeg(strcat(outputdir1,'compare_mc_',caseid(i),'_.jpeg'));
%end
%

%---------if single plot needed------------

%all_marks = {'o','v','d','^','s','>','x','v','>','<','p','h'};
%for i=2:size(caseid,2)
%figure(i+10),clf
%for j=1:size(Mc,2)
%Ly=Ly_mat(j);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 
%fname1 = strcat(dir,caseid(i),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(j)),'_',caseid(i),'.dat');
%f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
%fname2 = strcat(dir,caseid(i),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(i),'.dat');
%f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
%idx=find(abs(f2(:,1)-Mc(j))<0.001);del_tavg=f2(idx,2);
%plot(y/del_tavg,f1(:,8),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-','LineWidth',1.5,'DisplayName',strcat('M_c=',sprintf('%1.1f',Mc(j))));hold on;
%%plot(y/del_tavg,-f1(:,12),'Marker',all_marks{j},'Markerindices',1:15:length(y),'color',col{j},'LineStyle','-.','LineWidth',1.5,'HandleVisibility','off');hold on;
%end          %%%5-R11,8-R22,6-R12,11-prod,12-diss
%xlim([-10 10]);
%ax=gca;  ax.YRuler.Exponent = 0; 
%set(findall(gcf,'-property','Markersize'),'Markersize',4)
%set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',15);
%h=ylabel('\boldmath$\sqrt{|R_{22}|}/{\Delta{u}}$');set(h,'Interpreter','latex');set(ylab,'Fontsize',15);
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%%txt1 = ' Production(P)\rightarrow';
%%t = text(-4.5,0.0015,txt1,'Fontsize',9,'FontWeight','bold');
%%txt2 = 'Dissipation(\epsilon)\rightarrow';
%%t = text(-4.5,-0.0006,txt2,'Fontsize',9,'FontWeight','bold');
%legend show;legend('Location','Northeast','box','off');
%screen2jpeg(strcat(outputdir1,'R22_profiles_',caseid(i),'_.jpeg'));
%end



%%% If single plot needed----------------
lin={':','-.','-'};
col={myb,myb,myb,'r','r','r'};
all_marks = {'none','none','none'};
caseid=["lad","amd_dyn_c0","amd_dyn"];
caseid_name=["LAD","AMD-Dyn-c_{\mu}=0","AMD-Dyn"];
figure(20),clf
for i=1:5:size(Mc,2)
Ly=Ly_mat(i);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)';

for j=1:size(caseid,2)
fname1 = strcat(dir,caseid(j),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(i)),'_',caseid(j),'.dat');
%fname1 = strcat(dir,caseid(j),'_data/','precorr_tavg_',sprintf('%1.1f',Mc(i)),'_',caseid(j),'.dat');
f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
fname2 = strcat(dir,caseid(j),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(j),'.dat');
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);

plot(y/del_tavg,f1(:,12),'Marker',all_marks{j},'Markerindices',1:15:length(y),'Markerfacecolor',col{i},'MarkerEdgeColor',col{i},'LineStyle',lin{j},'LineWidth',2,'color',col{i},'DisplayName',strcat(caseid_name(j),'-',sprintf('%1.1f',Mc(i))));hold on
end  %%5-R11,8-R22,6-R12,11-prod,12-diss


ny_dns=1448;
dy_dns = Ly/ny_dns; y_dns = linspace(-Ly/2,Ly/2,ny_dns)'; 
%fname1_dns = strcat(dir1,'dns_pp_data/','vel_rey_nor_',sprintf('%1.1f',Mc(i)),'_','dns.dat');
fname1_dns = strcat(dir1,'dns_pp_data/','diss_nor_',sprintf('%1.1f',Mc(i)),'_','dns.dat');
f1_dns = readmatrix(fname1_dns,'ConsecutiveDelimitersRule','join');  
fname2_dns = strcat(dir1,'dns_pp_data/','delta_tavg_','dns.dat');
f2_dns = readmatrix(fname2_dns,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2_dns(:,1)-Mc(i)) < 0.001); %Mc,delta_theta
del_tavg_dns=f2_dns(idx,2);
mar={'o','s','s','s','s','s'};
plot(y_dns/del_tavg_dns,f1_dns(:,2),'marker',mar{i},'Markerfacecolor','none','MarkerEdgeColor','k','LineWidth',1.5,'Linestyle','none','Displayname',strcat('DNS-',sprintf('%1.1f',Mc(i))),'markerindices',1:40:length(y_dns));hold on


ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
xlim([-10 10]);%ylim([0 0.15]);
xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',20);
%h=ylabel('\boldmath$\sqrt{|R_{11}|}/{\Delta{u}}$');set(h,'Interpreter','latex');set(h,'Fontsize',20);
h=ylabel('\boldmath$\epsilon\delta_\theta/\Delta{u}^3$');set(h,'Interpreter','latex');set(h,'Fontsize',20);
legend show;leg=legend('Location','northeast','box','off','Fontsize',6,'Fontname','Timesnewroman','numcolumns',1);
leg.ItemTokenSize = [10 10];
ax=gca;  ax.YRuler.Exponent = -4; 
set(findall(gcf,'-property','Markersize'),'Markersize',6)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','FontName'),'FontName','Times');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir1,'diss_profile_amd.jpeg'));
end
%%%%%-----------------------------------------------------------------------------------------------------

