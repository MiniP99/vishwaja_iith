%%========================= Plot of SGS cofficients =============%% 

clear,clf
inputdir='/home/vsj/Codes/runs180/';
inputdir1='/home/vsj/Codes/runs180/smag_dyn_data/';
outputdir='/home/vsj/Codes/plot/myplots/runs180/Dynamic_plots_nor/';
mkdir(outputdir)
caseid='smag_dyn';
figure(1),clf
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={'r',myg,'b',myv,myy,'k'};
all_marks = {'o','v','d','^','s','>','x','v','>','<','p','h'};
mc=["02","04","08","12","16","20"];
Mc=[0.2,0.4,0.8,1.2,1.6,2.0];
Ly_mat=[200,200,100,100,80,80];
du=2*Mc*sqrt(1.4);
ny=240;


%t_i=[2700,2600,2000,1900,2600,2800];t_f=[5400,4900,5200,5400,6600,7300];    %%MGM
%t_i=[3300,2600,2800,2400,2400,2500];t_f=[6000,5800,6000,6000,7000,7000];    %%AMD
%t_i=[3300,3400,2800,2500,2500,3400];t_f=[5900,5700,6000,6000,6900,8400];    %%Sigma
t_i=[3300,0,0,0,0,3400];t_f=[5900,0,0,0,0,9600];    %%Smagorinsky

%% ---------------- Plot of averaged profiles along y at different Mc --------------------- %%
for i=1%:size(mc,2)
step = t_i(i):100:t_f(i);
Ly=Ly_mat(i);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 

for j= 1: size(step,2)
fname=strcat(inputdir,'run',mc(i),'/',caseid,'/cmodel_',sprintf('%05d',step(j)),'.dat');
f=importdata(fname);
%c_local(:,j)=f(:,2);c_local_qsgs(:,j)=f(:,3);
c_local(:,j)=f(:,2)*(1.5^(4/3));c_local_qsgs(:,j)=f(:,3)*(1.5^(4/3));  %Sigma,mgm
%c_local(:,j)=f(:,2)*(1.5*1.5);c_local_qsgs(:,j)=f(:,3)*(1.5*1.5);   %AMD
end
size(c_local)
fname2 = strcat(inputdir1,'delta_tavg_',caseid,'.dat');
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);y_c=f2(idx,7);
c_local_tavg=sum(c_local,2)/size(step,2);c_local_qsgs_tavg=sum(c_local_qsgs,2)/size(step,2);
%plot(f(8:end-8,1)./del_tavg,c_local_qsgs_tavg(8:end-8),'Marker',all_marks{i},'Markerindices',8:15:length(y)-16,'color',col{i},'LineWidth',1.5,'DisplayName',strcat(sprintf('%1.1f',Mc(i))));hold on
plot(f(:,1)./del_tavg,c_local_qsgs_tavg,'Marker',all_marks{i},'Markerindices',1:15:length(y),'color',col{i},'LineWidth',1.5,'DisplayName',strcat(sprintf('%1.1f',Mc(i))));
hold on

[value_yc,idx_yc] = min(abs(f(:,1)-y_c));  %%centreline y_c at max S
c_eps_centre(i)=c_local_tavg(idx_yc);c_epsT_centre(i)=c_local_qsgs_tavg(idx_yc);
dy = f(2,1)-f(1,1);
%c_eps_int(i) = abs(sum(c_local_tavg(26:end-26))*dy)/del_tavg; c_epsT_int(i) = abs(sum(c_local_qsgs_tavg(26:end-26))*dy)/del_tavg;
c_eps_int(i) = abs(sum(c_local_tavg)*dy)/del_tavg; c_epsT_int(i) = abs(sum(c_local_qsgs_tavg)*dy)/del_tavg;

y_nor(i)=(f(end,1)-f(1,1))/del_tavg;
%y_nor(i)=(f(end-26,1)-f(26,1))/del_tavg;

c_local=zeros(ny,size(step,2));c_local_qsgs=zeros(ny,size(step,2));


end
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',20);
ylab=ylabel('\boldmath$C_{\epsilon{T}}$','Interpreter','latex');set(ylab,'Fontsize',20);
xlim([-20 20])
%ylim([0 0.45])
legend show;legend('Location','Northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',6)
ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'c_epsT_tavg_profiles_',caseid,'.jpeg'));

c_eps_avg = c_eps_int./y_nor;c_epsT_avg = c_epsT_int./y_nor;

%%---write mach nos to a file----------
%fname = strcat(outputdir,'dyn_',caseid,'.dat');
%fileID = fopen(fname,'a');
%for i=1:size(Mc,2)
%fprintf(fileID,'%1.1f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f \r\n',Mc(i),c_eps_centre(i),c_epsT_centre(i),c_eps_int(i),c_epsT_int(i),c_eps_avg(i),c_epsT_avg(i));
%end
%fclose(fileID);

%%% ---------------- Plot of average values -------------------------------- %%
figure(2),clf  %change at fname,screen2jpeg
col={'r','b','g',myv,myg,'k'};
caseid = ["sigma_dyn","mgm_dyn","amd_dyn"];
const_name=["C_{\epsilon,c}","C_{\epsilon{T},c}","C_{\epsilon,int}","C_{\epsilon{T},int}","C_{\epsilon,avg}","C_{\epsilon{T},avg}"];
fname = strcat(outputdir,'dyn_',caseid(2),'.dat');
f=readmatrix(fname);
plot(f(:,1),f(:,6),'k','Marker',all_marks{1},'MarkerFaceColor',col{1},'MarkeredgeColor',col{1},'Linestyle','-','Linewidth',1.5,'DisplayName', const_name(5));hold on
plot(f(:,1),f(:,7),'k','Marker',all_marks{2},'MarkerFaceColor',col{2},'MarkeredgeColor',col{2},'Linestyle','-.','Linewidth',1.5,'DisplayName',const_name(6));hold on
%ylim([0.1 0.23])
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
ax=gca;  ax.YRuler.Exponent = 0;
legend show;legend('Location','Northeast','box','off','Fontsize',14);
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
xlab=xlabel('\boldmath$M_c$','interpreter','Latex');set(xlab,'Fontsize',20);
%ylab=ylabel(strcat('\boldmath$',const_name(5),'$'),'interpreter','latex');set(ylab,'Fontsize',20);
xticks([0.2,0.4,0.8,1.2,1.6,2.0]);
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',10)
screen2jpeg(strcat(outputdir,'c_avg_',caseid(2),'.jpeg'));



%%% ---------------- Plot of Integrated and centreline values -------------------------------- %%
%figure(2),clf
%myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
%caseid = ["amd_dyn","sigma_dyn","mgm_dyn"];
%caseid_name = ["AMD","Sigma","MGM"];
%col={'g','m','y'};
%all_marks = {'d','^','v'};
%line_sty = {'-','-.','--'};
%%%const_name=["c_{\epsilon,c}^{-2}","c_{\epsilon{T},c}^{-2}","c_{\epsilon,int}^{-2}","c_{\epsilon{T},int}^{-2}","c_{\epsilon,avg}^{-2}","c_{\epsilon{T},avg}^{-2}"];
%const_name=["c_{\epsilon,c}","c_{\epsilon{T},c}","c_{\epsilon,int}","c_{\epsilon{T},int}","c_{\epsilon,avg}","c_{\epsilon{T},avg}"];
%
%for i=1:size(caseid,2)
%fname = strcat(outputdir,'dyn_',caseid(i),'.dat');
%f=readmatrix(fname,'ConsecutiveDelimitersRule','join');
%plot(f(:,1),f(:,3),'k','Marker',all_marks{i},'MarkerFaceColor',col{i},'Linestyle',line_sty{i},'DisplayName',caseid_name(i));hold on
%end
%xlab=xlabel('\boldmath$M_c$','interpreter','latex');set(xlab,'Fontsize',15);
%ylab=ylabel(strcat('\boldmath$',const_name(2),'$'),'interpreter','latex');set(ylab,'Fontsize',15);
%ax=gca;  ax.YRuler.Exponent = 0;ylim([0 2])
%%legend show;legend('Position',[0.18 0.6 0.1 0.1],'box','off');
%legend show;legend('Location','northeast','box','off');
%set(findall(gcf,'-property','Markersize'),'Markersize',7)
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%xticks([0.2,0.4,0.8,1.2,1.6,2.0]);
%set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%set(findall(gcf,'-property','Markersize'),'Markersize',7)
%screen2jpeg(strcat(outputdir,'cT_centre','.jpeg'));


%%% ---------------- Plot of Integrated and centreline values -------------------------------- %%
%c_eps_avg = c_eps_int./y_nor;c_epsT_avg = c_epsT_int./y_nor;
%
%figure(2),clf
%const=[c_eps_centre',c_epsT_centre',c_eps_int',c_epsT_int',c_eps_avg',c_epsT_avg'];
%const_name=["c_{\epsilon,c}^{-2}","c_{\epsilon{T},c}^{-2}","c_{\epsilon,int}^{-2}","c_{\epsilon{T},int}^{-2}","c_{\epsilon,avg}^{-2}","c_{\epsilon{T},avg}^{-2}"];
%%const_name=["c_{\epsilon,c}","c_{\epsilon{T},c}","c_{\epsilon,int}","c_{\epsilon{T},int}","c_{\epsilon,avg}","c_{\epsilon{T},avg}"];
%for k=1:6
%subplot(2,3,k)
%scatter(Mc,const(:,k),'k','Marker',all_marks{k},'MarkerFaceColor',col{k});
%ax=gca;  ax.YRuler.Exponent = 0;
%grid on;hAx=gca;set(hAx,'xminorgrid','on','yminorgrid','on')
%xlabel('M_c','FontWeight','bold');ylabel(const_name(k),'FontWeight','bold');
%xticks([0.2,0.4,0.8,1.2,1.6,2.0]);
%end
%sgtitle('Dynamic Sigma');
%set(findall(gcf,'-property','FontSize'),'FontSize',6)
%set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%screen2jpeg(strcat(outputdir,'c_subplot_',caseid,'.jpeg'));
