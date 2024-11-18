%%======================= Check self-similar state ===========================%%
clear,clf    %%% mainid,runid,caseid,mc,times(2)
mainid = "runs180";
runid  = "run08";
caseid = "smag_dyn";

c=sqrt(1.4);Mc=0.8;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;

inputdir =strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',runid,'/',caseid,'/');
outputdir=strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',runid,'/',caseid,'/postprocess/');
mkdir(outputdir);

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
%dx = Lx/nx; dy = Ly/ny; dz = Lz/nz;
dx = X(2,1,1)-X(1,1,1); dy = Y(1,2,1)-Y(1,1,1);dz = Z(1,1,2)-Z(1,1,1);

 
%%----------- Plot of momentum thickness along with linear fit --------
fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
f=readmatrix(fname);
idx1=find(f(:,1)==300);idx2=find(f(:,1)==500);idx3=find(f(:,1)==500);

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
time=300:20:500;
tau=(time*du)/delta_0;
for i=1:size(time,2)
f3=zeros(ny,17);
fname=strcat(outputdir,'vel_tke_flux_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',time(i)),'.dat');
f3=dlmread(fname);  %%  1-ub,    2-util,  3-vbar,  4-vtil,  5-R11,  6-R12, 7-R13, 8-R22, 9-R23, 10-R33
idx=find(t1==time(i));    %%  11-prod,12-eps,   13-trans, 14-baro, 15-pre_dil,  16-mass_ru,  17-mass_rv
del=d1(idx);
xaxis=y/del;
all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};
txt =strcat('\tau=',sprintf('%4.0f',tau(i)));

% plot velocity---------------        'Marker',all_marks{i},'MarkerIndices',1:5:length(y),
%plot(y,f3(:,1),'DisplayName',txt,'LineWidth',1);hold on

%plot reynold stress -------
plot(xaxis,f3(:,5),'DisplayName',txt,'LineWidth',1);hold on

%plot dissipation----------
%plot(xaxis,f3(:,12),'DisplayName',txt,'LineWidth',1);hold on
end

xlim([-15 15]);
ax=gca;  ax.YRuler.Exponent = 0;  xlabel('y/{\delta_\theta(\tau)}');

%h=ylabel('$\overline{u}/\Delta u$');set(h,'Interpreter','latex')

h=ylabel('\boldmath $\sqrt{|R_{11}|}/{\Delta{u}}$');set(h,'Interpreter','latex')

%ylabel('\epsilon{\delta_\theta}/{\Delta{u}^3}');

set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
legend show;legend('box','off');

%screen2jpeg(strcat(outputdir,'ubar_y_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
screen2jpeg(strcat(outputdir,'R11_y_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
%screen2jpeg(strcat(outputdir,'diss_y_',sprintf('%1.1f',Mc),'_',caseid,'.png'));



%%%-----------------------------------------------------------------------------------\\

%% Plotting integrated values
clf
figure(3)
fname=strcat(outputdir,'Rij_P12_dis_pi_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
f4=dlmread(fname);  %% 1-time 2-R11 3-R12 4-R13 5-R22 6-R23 7-R33 8-TKE 9-P12 10-diss 11-pi11 12-pi12 13-pi22
tau=(f4(:,1)*du)/delta_0;
yyaxis left
plot(tau(1:idx3),f4(1:idx3,2),'-ob','DisplayName','R_{11,int}');hold on
plot(tau(1:idx3),-f4(1:idx3,3),'-+m','DisplayName','R_{12,int}');hold on
plot(tau(1:idx3),f4(1:idx3,5),'-*g','DisplayName','R_{22,int}');hold on
xline(tau(idx1),'--k',{'self similar'},'Handlevisibility','off');ax=gca;  ax.YRuler.Exponent = 0;  
ylabel('R_{ij,int}')
yyaxis right
ylabel('\epsilon_{int}')
plot(tau(1:idx3),f4(1:idx3,10),'-sr','DisplayName','\epsilon_{int}');hold on;ax=gca;ax.YAxis(2).Exponent = 0;
p=patch([tau(idx1) tau(idx1) tau(idx2) tau(idx2)],[min(ylim) max(ylim) max(ylim) min(ylim)],[0.9 0.9 0.9],'EdgeColor','None','HandleVisibility','off');alpha(p,0.5)
xlabel('\tau=(t*du)/\delta_{\theta}_0');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
legend('Location','southeast','box','off');legend show
screen2jpeg(strcat(outputdir,'Rij_diss_int_',sprintf('%1.1f',Mc),'_',caseid,'.png'));


