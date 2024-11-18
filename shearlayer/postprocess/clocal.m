clear,clf
inputdir='/home/vsj/Codes/shearlayer_les/runs180/run20/amd_dyn/';
outputdir='/home/vsj/Codes/shearlayer_les/runs180/run20/amd_dyn/postprocess/';
%inputdir='/home/vsj/Codes/check_mgm_dyn/runs128/';
%outputdir='/home/vsj/Codes/check_mgm_dyn/runs128/postprocess/';
mkdir(outputdir)
caseid = "amd_dyn";            %caseid,Mc,domain,time
step=2500:100:7000;
figure(1),clf
for i=1:size(step,2)
fname = strcat(inputdir,'cmodel_',sprintf('%05d',step(i)),'.dat');
%f=readmatrix(fname,'Delimiter', ' ','ConsecutiveDelimitersRule', 'join');
f=importdata(fname);
%c_local(:,i)=f(:,2);c_local_qsgs(:,i)=f(:,3);c_local_tke(:,i)=f(:,4);
%c_local(:,i)=f(:,2)*(1.5^(4/3));c_local_qsgs(:,i)=f(:,3)*(1.5^(4/3));%c_local_tke(:,i)=f(:,4)*(1.5^(4/3));   %Sigma & MGM
c_local(:,i)=f(:,2)*(1.5*1.5*12);c_local_qsgs(:,i)=f(:,3)*(1.5*1.5*12);c_local_tke(:,i)=f(:,4)*(1.5*1.5*12);        %AMD
p=plot(f(:,1),f(:,2),'LineWidth',1,'HandleVisibility','off');p.Color(4)=0.2;hold on
%p=plot(f(:,1),f(:,2),'LineWidth',1,'DisplayName',sprintf('%4d',step(i)));p.Color(4)=0.5;hold on
%p=plot(f(:,1),f(:,2)*(1.5^(4/3)),'LineWidth',1,'HandleVisibility','off');p.Color(4)=0.2;hold on
%p=plot(f(:,1),f(:,2)*(1.5*1.5*12),'LineWidth',1,'HandleVisibility','off');p.Color(4)=0.2;hold on
end
c_local_tavg=sum(c_local,2)/size(step,2);c_local_qsgs_tavg=sum(c_local_qsgs,2)/size(step,2);
c_local_tke_tavg=sum(c_local_tke,2)/size(step,2);
plot(f(:,1),c_local_tavg,'r','LineWidth',2,'DisplayName','c_{\epsilon}' );hold on
plot(f(:,1),c_local_qsgs_tavg,'b','LineWidth',2,'DisplayName','c_{\epsilon T}' );hold on;
plot(f(:,1),c_local_tke_tavg,'g','LineWidth',2,'DisplayName','C_I' );hold on;
c_local_yavg = mean(c_local_tavg);c_local_qsgs_yavg = mean(c_local_qsgs_tavg);
c_local_tke_yavg = mean(c_local_tke_tavg);
txt1=strcat(sprintf('%1.4f',c_local_yavg));
txt2=strcat(sprintf('%1.4f',c_local_qsgs_yavg));
txt3=strcat(sprintf('%1.4f',c_local_tke_yavg));
yline(c_local_yavg,'--r',txt1,'LineWidth',1,'HandleVisibility','off');hold on;
yline(c_local_qsgs_yavg,'--b',txt2,'LineWidth',1,'HandleVisibility','off');
yline(c_local_tke_yavg,'--g',txt3,'LineWidth',1,'HandleVisibility','off');
%ylim([0 0.05])
legend show;legend('Location','Northeast','box','off');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','on','yminorgrid','on')
screen2jpeg(strcat(outputdir,'c_model.jpeg'));
