%%%===================Static-Dynamic Comparision==============================%%%
%change outputdir1,caseid,caseid_name,caseid_image,col,col1,all_marks

clear,clf
inputdir='/home/vsj/Codes/plot/';
outputdir='/home/vsj/Codes/plot/myplots/';
outputdir1='/home/vsj/Codes/plot/myplots/runs180/sta_dyn_mgm_plots/';
dir = '/home/vsj/Codes/runs180/';
mkdir(outputdir1)

kris_Rij_int = readmatrix(strcat(outputdir,'kris_rij_int.csv'));
kris_p12     = readmatrix(strcat(outputdir,'kris_intp12.csv'));
kris_diss    = readmatrix(strcat(outputdir,'kris_intdiss.csv'));
kris_pi11p12 = readmatrix(strcat(outputdir,'kris_pi11p12.csv'));
kris_pi12p12 = readmatrix(strcat(outputdir,'kris_pi12p12.csv'));
kris_pi22p12 = readmatrix(strcat(outputdir,'kris_pi22p12.csv'));
kris_R11 = readmatrix(strcat(outputdir,'kris_R11_peak.csv'));
kris_R22 = readmatrix(strcat(outputdir,'kris_R22_peak.csv'));
kris_R12 = readmatrix(strcat(outputdir,'kris_R12_peak.csv'));


%-----------------------------------------------------------------------------------------------------
%%  1-ub,    2-util,  3-vbar,  4-vtil,  5-R11,  6-R12, 7-R13, 8-R22, 9-R23, 10-R33
%%  11-prod,12-eps,   13-trans, 14-baro, 15-pre_dil,  16-mass_ru,  17-mass_rv
caseid=["mgm_dyn","mgm_sta"];caseid_name=["MGM-Dynamic","MGM-Static"];caseid_image=["mgm"];
Mc = [0.2,0.4,0.8,1.2,1.6,2.0]; c=sqrt(1.4);du = 2*Mc*c;
Ly_mat = [200,200,100,100,80,80];ny = 240;         %%%-------------change  
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={'b','none'};  %change marker
col1={'b','k'};  
all_marks = {'v','v'};

%budget_label=["P_{12}","\epsilon","\pi_{11}/\pi_{11,0.2}","\pi_{12}/\pi_{12,0.2}","\pi_{22}/\pi_{22,0.2}","\sqrt{|R_{11}|}/\Delta u","\sqrt{|R_{12}|}/\Delta u","\sqrt{|R_{22}|}/\Delta u"];
budget_label=["P_{12}","\epsilon","|\pi_{11}|","|\pi_{12}|","|\pi_{22}|","\sqrt{|R_{11}|}/\Delta u","\sqrt{|R_{12}|}/\Delta u","\sqrt{|R_{22}|}/\Delta u"];

%budget_matsuno = [kris_p12(:,2),kris_diss(:,2),(kris_pi11p12(:,2).*kris_p12(:,2))./(kris_pi11p12(1,2).*kris_p12(1,2)),(kris_pi11p12(:,2).*kris_p12(:,2))./(kris_pi12p12(1,2).*kris_p12(1,2)),(kris_pi22p12(:,2).*kris_p12(:,2))./(kris_pi22p12(1,2).*kris_p12(1,2)),kris_R11(:,2),kris_R12(:,2),kris_R22(:,2)];
budget_matsuno = [kris_p12(:,2),kris_diss(:,2),kris_pi11p12(:,2).*kris_p12(:,2),kris_pi11p12(:,2).*kris_p12(:,2),kris_pi22p12(:,2).*kris_p12(:,2),kris_R11(:,2),kris_R12(:,2),kris_R22(:,2)];

im_txt_tke=["P12","diss","pi11","pi12","pi22"];
for k=1:5
figure(k),clf
plot(Mc,budget_matsuno(:,k),'Marker','o','MarkerFaceColor','k','MarkerEdgecolor','k','Displayname','DNS');hold on
for i=1:size(caseid,2)
fname = strcat(dir,caseid(i),'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid(i),'.dat');
f = readmatrix(fname,'ConsecutiveDelimitersRule','join'); 
%1-Mc,2-del_avg,3-R11,4-R12,5-R13,6-R22,7-R23,8-R33,9-ke,10-P12,11-diss,12-pi11,13-pi12,14-pi22
f_new = [f(:,10),f(:,11),-f(:,12),f(:,13),f(:,14)];
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
plot(f(:,1),f_new(:,k),'Marker',all_marks{i},'MarkerFaceColor',col{i},'MarkeredgeColor',col1{i},'Linestyle','none','Displayname',caseid_name(i));hold on
h=ylabel(strcat('\boldmath$',budget_label(k),'$'));set(h,'Interpreter','latex','Fontsize',20);
end
legend show;legend('Location','Northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
xticks([0.2,0.4,0.8,1.2,1.6,2.0]);ax=gca;  ax.YRuler.Exponent = 0; 
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
ax=gca;  ax.YRuler.Exponent = 0;
set(findall(gcf,'-property','Markersize'),'Markersize',10)
set(findall(gcf,'-property','Linestyle'),'Linestyle','none')
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
xlab=xlabel('\boldmath$M_c$','Interpreter','Latex');set(xlab,'Fontsize',20);
screen2jpeg(strcat(outputdir1,im_txt_tke(k),'_',caseid_image,'.jpeg'));
end

R_max(:,1) = Mc ;
im_txt=["R11","R12","R22"];
for k=6:8
figure(k),clf
plot(R_max(:,1),budget_matsuno(:,k),'Marker','o','MarkerFaceColor','k','MarkerEdgecolor','k','Displayname','DNS');hold on
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
for i=1:size(caseid,2)
for j=1:size(Mc,2)
fname = strcat(dir,caseid(i),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(j)),'_',caseid(i),'.dat');
f_R= readmatrix(fname,'ConsecutiveDelimitersRule','join'); % %%5-R11,8-R22,6-R12
R_max(j,2)= max(abs(f_R(:,5)));R_max(j,3)= max(abs(f_R(:,6)));R_max(j,4)= max(abs(f_R(:,8)));
end

plot(R_max(:,1),R_max(:,k-4),'Marker',all_marks{i},'MarkerFaceColor',col{i},'MarkerEdgeColor',col1{i},'Displayname',caseid_name(i));hold on
end
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
h=ylabel(strcat('\boldmath$',budget_label(k),'$'));set(h,'Interpreter','latex','Fontsize',20);
legend show;legend('Location','Northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
xticks([0.2,0.4,0.8,1.2,1.6,2.0]);ax=gca;  ax.YRuler.Exponent = 0; 
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
ax=gca;  ax.YRuler.Exponent = 0;
set(findall(gcf,'-property','Markersize'),'Markersize',10)
set(findall(gcf,'-property','Linestyle'),'Linestyle','none')
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
xlab=xlabel('\boldmath$M_c$','Interpreter','Latex');set(xlab,'Fontsize',20);
screen2jpeg(strcat(outputdir1,im_txt(k-5),'_peak_',caseid_image,'.jpeg'));

end
%%%%-----------------------------------------------------------------------------------------------------
%%plot of mach
kris_y0   = readmatrix(strcat(outputdir,'kris_y0.csv'));
kris_udel = readmatrix(strcat(outputdir,'kris_udel.csv'));
kris_mach = readmatrix(strcat(outputdir,'kris_mach.csv'));
kris_mom  = readmatrix(strcat(outputdir,'matsuno_fit.csv'));

matsuno_corr = [kris_y0(1:6,2),kris_udel(1:6,2),kris_mom(:,2)];matsuno_label=["\delta_y/\delta_{99}","U_\delta / \Delta u","\dot{\delta_{\theta}}/\dot{\delta}_{\theta,inc}"];

im_txt_scale=["y0","udel","mom_rate"];
for k=1:3
figure(k),clf
plot(Mc',matsuno_corr(:,k),'Marker','o','MarkerFaceColor','k','MarkerEdgecolor','k','Displayname','DNS');hold on
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
for i=1:size(caseid,2)
fname2 = strcat(dir,caseid(i),'_data/','delta_tavg_',caseid(i),'.dat'); %1-Mc,2-del,3-del_ome,4-del_99,5-del_y,6-del_dot,7-y_c,8-U_del
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join');
f2_new = [f2(:,5)./f2(:,4),f2(:,8),f2(:,6)/0.018];
plot(f2(:,1),f2_new(:,k),'Marker',all_marks{i},'MarkerFaceColor',col{i},'MarkerEdgeColor',col1{i},'Displayname',caseid_name(i));hold on
end
ylab=ylabel(strcat('\boldmath$',matsuno_label(k),'$'), 'Interpreter','latex');set(ylab,'Fontsize',20);
legend show;legend('Location','Northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
xticks([0.2,0.4,0.8,1.2,1.6,2.0]);ax=gca;  ax.YRuler.Exponent = 0; 
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
ax=gca;  ax.YRuler.Exponent = 0;
set(findall(gcf,'-property','Markersize'),'Markersize',10)
set(findall(gcf,'-property','Linestyle'),'Linestyle','none')
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
xlab=xlabel('\boldmath$M_c$','Interpreter','Latex');set(xlab,'Fontsize',20);
screen2jpeg(strcat(outputdir1,im_txt_scale(k),'_',caseid_image,'.jpeg'));
end

%%%%------------------------ subplot of R11,R12,R22,prod-diss of different models on single plot-------------------------
%line_style=['-','--','-.'];
%all_marks = {'o','>','d','s','^','v','x','*','+','<','p','h'};
%
%for i=1:size(Mc,2)
%figure(i+3),clf 
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
%plot(y/del_tavg,sub_data(:,k),'color',col{j},'Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'DisplayName',caseid_name(j));hold on
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%xlabel('y/{\delta_\theta(\tau)}');
%xlim([-10 10]);
%h=ylabel(strcat('\boldmath$',sub_txt(k),'$'));set(h,'Interpreter','latex');
%text(0.025,0.95,charlbl{k},'Units','normalized');
%end
%
%subplot(2,2,4)
%plot(y/del_tavg,f1(:,11),'color',col{j},'LineStyle','-','Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'DisplayName',caseid_name(j));hold on;xlim([-7 7]);
%plot(y/del_tavg,-f1(:,12),'color',col{j},'LineStyle','--','Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'HandleVisibility','off');hold on;
%xlabel('y/{\delta_\theta(\tau)}');
%xlim([-6 6]);
%text(0.025,0.95,charlbl{4},'Units','normalized');
%end
%ax=gca;  ax.YRuler.Exponent = 0;
%set(findall(gcf,'-property','FontSize'),'FontSize',6)
%set(findall(gcf,'-property','Fontweight'),'Fontweight','bold')
%set(findall(gcf,'-property','Markersize'),'Markersize',2)
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
%
%%----------------------------------------- Single model at different mc subplot--------------------

%for i=1:size(caseid,2)
%figure(i+10),clf
%for j=1:size(Mc,2)
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
%plot(y/del_tavg,sub_data(:,k),'color',col{j},'LineWidth',1,'Marker',all_marks{j},'Markerindices',1:15:length(y),'DisplayName',strcat('M_c=',sprintf('%1.1f',Mc(j))));hold on;
%xlabel('y/{\delta_\theta(\tau)}');
%xlim([-10 10]);
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%h=ylabel(strcat('\boldmath$',sub_txt(k),'$'));set(h,'Interpreter','latex');
%text(0.025,0.95,charlbl{k},'Units','normalized');
%end
%
%subplot(2,2,4)
%plot(y/del_tavg,f1(:,11),'color',col{j},'LineStyle','-','Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'DisplayName',strcat('M_c=',sprintf('%1.1f',Mc(j))));hold on;
%plot(y/del_tavg,-f1(:,12),'color',col{j},'LineStyle','--','Marker',all_marks{j},'Markerindices',1:15:length(y),'LineWidth',1,'HandleVisibility','off');hold on;
%xlabel('y/{\delta_\theta(\tau)}');
%xlim([-6 6]);
%text(0.025,0.95,charlbl{4},'Units','normalized');
%end
%ax=gca;  ax.YRuler.Exponent = 0;
%set(findall(gcf,'-property','FontSize'),'FontSize',6)
%set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%set(findall(gcf,'-property','Markersize'),'Markersize',2)
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
