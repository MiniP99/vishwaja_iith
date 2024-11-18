%%======================= Check self-similar state ===========================%%
clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
caseid = "dns";

c=sqrt(1.4);Mc=1.2;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;
delta_0=1;      %delta_0=Re_theta/(rho_0*Re*du);

inputdir =strcat('/home/vsj/dns_data_csl/mc_1p2/postprocess/');
outputdir=strcat('/home/vsj/dns_data_csl/mc_1p2/postprocess/');
mkdir(outputdir);

%%------------- Read Coordinates x, y, z-----------------------------
Lx = 100; Ly = 100; Lz = 50;
nx = 1024; ny = 1448; nz = 512;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=1.2;du=Mc*2*c;rho_0=1;Re_theta=1000;
delta_0=1;

 
%%----------- Plot of momentum thickness along with linear fit --------
fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
f=readmatrix(fname);
idx1=find(f(:,1)==50);idx2=find(f(:,1)==120);idx3=find(f(:,1)==120);

tau=(f(1:idx3,1)*du)/delta_0;
del=f(1:idx3,2)/delta_0;
plot(tau,del,'-or');
hold on

tau=(f(idx1:idx2,1)*du)/delta_0;
del=f(idx1:idx2,2)/delta_0;
eq=polyfit(tau,del,1);
delta1 = polyval(eq,tau);
plot(tau,delta1,'--b','LineWidth',2);
title('momentum thickness vs time'), 
xlabel('(t*du)/\delta_{\theta}_0');ylabel('\delta_{\theta}/ \delta_{\theta}_0');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
legend({strcat('Mc=',sprintf('%1.1f',Mc)),strcat('Linearfit:y=',sprintf('%1.3f',eq(1,1)),'x+(',sprintf('%1.3f',eq(1,2)),')')},'Location','best','box','off');
grid on;hAx=gca; set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'mom_thick_',sprintf('%1.1f',Mc),'_',caseid,'.png'));

%%-----------------------------------------------------------------------------------------------
fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
delta1= readmatrix(fname);t1=delta1(:,1);d1=delta1(:,2);d2=delta1(:,4);

clf
figure(2)
time=50:10:120;
tau=(time*du)/delta_0;
for i=1:size(time,2)
%f3=zeros(ny,17);
fname=strcat(outputdir,'precorr_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',time(i)),'.dat');
f3=dlmread(fname);  %%  1-ub,2-util, 3-R11,  4-R12,  5-R22
idx=find(t1==time(i));    
del=d1(idx);
xaxis=y/del;
all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};
txt =strcat('\tau=',sprintf('%4.0f',tau(i)));

% plot velocity---------------        'Marker',all_marks{i},'MarkerIndices',1:5:length(y),
%plot(y,f3(:,1),'DisplayName',txt,'LineWidth',1);hold on

p11_int= sum(abs(f3(:,2)))*(dy/del);
p12_int= sum(abs(f3(:,3)))*(dy/del);
p22_int= sum(abs(f3(:,4)))*(dy/del);

fname=strcat(outputdir,'precorr_int_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
fileID=fopen(fname,'a');
fprintf(fileID,'%04d %12.8f %12.8f %12.8f\n',time(i),p11_int,p12_int,p22_int);
fclose(fileID);

%plot reynold stress -------
plot(xaxis,-f3(:,2),'DisplayName',txt,'LineWidth',1);hold on

%plot dissipation----------
%plot(xaxis,f3(:,12),'DisplayName',txt,'LineWidth',1);hold on
end

xlim([-10 10]);
ax=gca;  ax.YRuler.Exponent = 0;  xlabel('y/{\delta_\theta(\tau)}');

%h=ylabel('$\overline{u}/\Delta u$');set(h,'Interpreter','latex')

%h=ylabel('\boldmath $\sqrt{|R_{12}|}/{\Delta{u}}$');set(h,'Interpreter','latex')

ylabel('\pi_{11}{\delta_\theta}/{\Delta{u}^3}');

set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
legend show;legend('box','off');

%screen2jpeg(strcat(outputdir,'ubar_y_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
screen2jpeg(strcat(outputdir,'pi11_y_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
%screen2jpeg(strcat(outputdir,'diss_y_',sprintf('%1.1f',Mc),'_',caseid,'.png'));






