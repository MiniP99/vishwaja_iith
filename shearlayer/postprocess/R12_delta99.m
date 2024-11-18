clear,clf
inputdir='/home/vsj/Codes/plot/';
outputdir='/home/vsj/Codes/plot/myplots/';
outputdir1='/home/vsj/dns_data_csl/plots/';
dir = '/home/vsj/Codes/runs180/';
dir1 = '/home/vsj/dns_data_csl/';
mkdir(outputdir1)
caseid=["base","lad","amd_dyn","sigma_dyn","mgm_dyn"];caseid_name=["Base","LAD","AMD-Dyn","Sigma-Dyn","MGM-Dyn"];
caseid_dns="dns";

%Mc = [0.2,0.4,0.8,1.2,1.6,2.0]; c=sqrt(1.4);
Mc = [0.2,1.2,1.6,2.0]; c=sqrt(1.4);
du = 2*Mc*c;du=du';
%Ly_mat = [200,200,100,100,80,80];ny = 240;ny_dns=1448;         %%%-------------change 
Ly_mat = [200,100,80,80];ny = 240;ny_dns=1448;         %%%-------------change 
 
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];

%col={'r',myg,'b',myv,myy,'k'}; %mach
%all_marks = {'o','v','d','^','s','>','x','v','>','<','p','h'};
col={'r',myv,myy,'k'}; %mach
all_marks = {'o','^','s','>','x','v','>','<','p','h'};

figure(1),clf
%%  1-ub,    2-util,  3-vbar,  4-vtil,  5-R11,  6-sqrt(R12)/du, 7-R13, 8-R22, 9-R23, 10-R33
%%  11-prod,12-eps,   13-trans, 14-baro, 15-pre_dil,  16-mass_ru,  17-mass_rv
for j=2:size(caseid,2)
figure(j),clf
for i=1:size(Mc,2)
Ly=Ly_mat(i);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 
fname1 = strcat(dir,caseid(j),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc(i)),'_',caseid(j),'.dat');
f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
fname2 = strcat(dir,caseid(j),'_data/','delta_tavg_',caseid(j),'.dat');
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001); %Mc,delta_theta,delta_omega,delta_99
del_tavg=f2(idx,2);
plot(y/del_tavg,f1(:,8),'color',col{i},'LineWidth',2,'DisplayName',strcat(sprintf('%1.1f',Mc(i))));hold on


Ly=Ly_mat(i);dy = Ly/ny_dns; y_dns = linspace(-Ly/2,Ly/2,ny_dns)'; 
fname1_dns = strcat(dir1,caseid_dns,'_pp_data/','vel_rey_nor_',sprintf('%1.1f',Mc(i)),'_',caseid_dns,'.dat');
f1_dns = readmatrix(fname1_dns,'ConsecutiveDelimitersRule','join');  
fname2_dns = strcat(dir1,caseid_dns,'_pp_data/','delta_tavg_',caseid_dns,'.dat');
f2_dns = readmatrix(fname2_dns,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2_dns(:,1)-Mc(i)) < 0.001); %Mc,delta_theta
del_tavg_dns=f2_dns(idx,2);
plot(y_dns/del_tavg_dns,f1_dns(:,5),'marker',all_marks{i},'color',col{i},'LineWidth',2,'Linestyle','none','HandleVisibility','off','markerindices',1:40:length(y_dns));hold on

end
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
xlim([-10 10])
xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',20);
ylab=ylabel('\boldmath$\sqrt{|R_{22}|}/{\Delta{u}}$','Interpreter','latex');set(ylab,'Fontsize',20);
legend show;legend('Location','Northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',6)
ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir1,'R22_',caseid(j),'.jpeg'));
end
%%%%-----------------------------------------------------------------------------------------------------


%%%%-----------------------------------------------------------------------------------------------------
%%plot of mach
%col={'r',myb,'g','m','y'};
%all_marks = {'o','s','d','^','v'};
%
%kris_y0   = readmatrix(strcat(outputdir,'kris_y0.csv'));
%kris_udel = readmatrix(strcat(outputdir,'kris_udel.csv'));
%kris_mach = readmatrix(strcat(outputdir,'kris_mach.csv'));
%kris_mom  = readmatrix(strcat(outputdir,'matsuno_fit.csv'));
%
%
%figure(2),clf
%txt_mach=["M_t","M_{tv}","M_\tau","M_g"];
%for k=4
%for i=1:size(caseid,2)
%fname = strcat(dir,caseid(i),'_data/','mach_tavg_',caseid(i),'.dat'); %1-Mc, 2-M_t, 3-M_tv, 4-M_tau, 5-M_g
%f = readmatrix(fname,'ConsecutiveDelimitersRule','join'); 
%scatter(f(:,1),f(:,k+1),'k','Marker',all_marks{i},'MarkerFaceColor',col{i},'Displayname',caseid_name(i));hold on
%end
%j=(6*k)-5;
%scatter(kris_mach(j:j+5,1),kris_mach(j:j+5,2),'k*','DisplayName','DNS');hold on
%xticks([0.2,0.4,0.8,1.2,1.6,2.0]);ax=gca;  ax.YRuler.Exponent = 0; 
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%xlabel('\boldmath$M_c$','Interpreter','latex');
%ylabel(strcat('\boldmath','$',txt_mach(k),'$'), 'Interpreter','latex');
%legend show;legend('Location','Northwest','box','off');
%set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%screen2jpeg(strcat(outputdir1,'M_g.png'));
%end
%
%

%%%----------------------
%figure(3),clf
%matsuno_corr = [kris_y0(1:6,2),kris_udel(1:6,2),kris_mom(:,2)];matsuno_label=["\delta_y/\delta_{99}","U_\delta / \Delta u","\dot{\delta}/\dot{\delta}_{inc}"];
%for k=6  %5:7
%for i=1:size(caseid,2)
%fname2 = strcat(dir,caseid(i),'_data/','delta_tavg_',caseid(i),'.dat'); %1-Mc,2-del,3-del_ome,4-del_99,5-del_y,6-del_dot,7-y_c,8-U_del
%f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join');
%f2_new = [f2(:,5)./f2(:,4),f2(:,8),f2(:,6)/0.018];
%scatter(f2(:,1),f2_new(:,k-4),'k','Marker',all_marks{i},'MarkerFaceColor',col{i},'Displayname',caseid_name(i));hold on
%end
%scatter(Mc',matsuno_corr(:,k-4),'k*','DisplayName','DNS');hold on
%xticks([0.2,0.4,0.8,1.2,1.6,2.0]);ax=gca;ax.YRuler.Exponent = 0; grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%xlabel('\boldmath$M_c$','Interpreter','latex');
%ylabel(strcat('\boldmath$',matsuno_label(k-4),'$'), 'Interpreter','latex');
%legend show;legend('Location','Northeast','box','off');
%set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%screen2jpeg(strcat(outputdir1,'delta_U.png'));
%end
