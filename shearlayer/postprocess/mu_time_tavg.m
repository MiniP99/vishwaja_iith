clear all,clf

inputdir = '/home/vsj/Codes/diss_time/';
outputdir = '/home/vsj/Codes/diss_time/plots/';
caseid = 'mgm_dyn';
Mc = [0.2,0.4,0.8,1.2,1.6,2.0];
c=sqrt(1.4);du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;
%
%%t_i = [400,280,150,120,0,0]; t_f = [900,600,350,300,0,0];
%%t_i = [400,280,150,160,140,110];t_f=[900,600,350,340,280,250];     %lad
%%t_i = [560,580,200,150,110,100];t_f=[1000,880,400,330,270,240];   %AMD
%%t_i = [560,480,200,150,110,130];t_f=[1000,800,400,330,270,280];   %sigma
%%t_i=[560,0,200,210,140,130];t_f=[1000,0,400,380,300,280];         %Smag
t_i=[460,380,150,120,110,110];t_f=[900,700,350,300,250,250];  %mgm
diss_tavg= zeros(size(Mc,2),4);
%
for i = 1: size(Mc,2)
fname = strcat(inputdir,'diss_time_',caseid,'_',sprintf('%1.1f',Mc(i)),'.dat');
f = readmatrix(fname,'ConsecutiveDelimitersRule','join');
f(:,3) = f(:,3)./f(:,2);f(:,4) = f(:,4)./f(:,2);
f(:,2) = f(:,2)./f(:,2);
display(f)
%
idxi=find(abs(f(:,1)-t_i(i)) < 0.001);idxf=find(abs(f(:,1)-t_f(i)) < 0.001);
diss = mean(f(idxi:idxf,:));
diss_tavg(i,:) = [Mc(i),diss(2),diss(3),diss(4)];
end
%
fname = strcat('/home/vsj/Codes/diss_time/diss_tavg_',caseid,'.dat');
writematrix(diss_tavg,fname,'Delimiter','tab');



fname = strcat(inputdir,caseid,'/','diss_time_',caseid,'_',sprintf('%1.1f',Mc),'.dat');
f = readmatrix(fname,'ConsecutiveDelimitersRule','join');
plot(f(2:end,1)*du,f(2:end,2),'-or','DisplayName','\boldmath$\epsilon^{\mu}_{int}$');hold on
plot(f(2:end,1)*du,f(2:end,3),'>b','DisplayName','\boldmath$\epsilon^{*}_{int}$');hold on
plot(f(2:end,1)*du,f(2:end,4),'^g','DisplayName','\boldmath$\epsilon^{sgs}_{int}$');hold on

%%--3plots:legend:17,ax.=20,label=22;

ax=gca;ax.YAxis.FontSize = 20;ax.XAxis.FontSize = 20;
xlab=xlabel('\boldmath$\tau$','Interpreter','latex');set(xlab,'Fontsize',22);
legend show;leg=legend('Location','northeast','box','off','Fontsize',17,'Interpreter','latex');
leg.ItemTokenSize = [10 10];
ax=gca;  ax.YRuler.Exponent = -3; 
set(findall(gcf,'-property','Markersize'),'Markersize',6)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Linestyle'),'Linestyle','None')
set(findall(gcf,'-property','FontName'),'FontName','Times');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'diss_time_',caseid,'_',sprintf('%1.1f',Mc),'.png'));
