clear,clf
inputdir='/home/vsj/Codes/runs180/';
outputdir='/home/vsj/Codes/plot/myplots/';
outputdir1='/home/vsj/Codes/plot/myplots/runs180/';
mkdir(outputdir)
delta_0=1;
del_the = 0.018;

caseid = ["amd_dyn","sigma_dyn","mgm_dyn","lad","base"];
caseid_name = ["AMD-Dyn","Sigma-Dyn","MGM-Dyn","LAD","Base"];

myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
%col={'r',myb,'g','m','y'};%r,myb
col={myg,'r','b','none','none'};
col1={myg,'r','b','k','k'};
all_marks = {'d','^','v','s','p'};

%delta_base   = [0.2,0.015;0.4,0.012;0.8,0.006];
%delta_lad    = [0.2,0.015;0.4,0.012;0.8,0.007;1.2,0.007;1.6,0.006;2.0,0.004];
%delta_amd    = [0.2,0.019;0.4,0.013;0.8,0.010;1.2,0.008;1.6,0.007;2.0,0.007];
%delta_sigma  = [0.2,0.016;0.4,0.016;0.8,0.009;1.2,0.007;1.6,0.006;2.0,0.006];
%delta_mgm    = [0.2,0.018;0.4,0.014;0.8,0.010;1.2,0.008;1.6,0.006;2.0,0.006];

kris_data    = readmatrix(strcat(outputdir,'matsuno_fit.csv'));
langley_data = readmatrix(strcat(outputdir,'langley_fit.csv'));

figure(1),clf
plot(kris_data(:,1),kris_data(:,2),'Marker','o','MarkerFaceColor','k','MarkerEdgecolor','k','Displayname','DNS');hold on
for i=1:size(caseid,2)
fname = strcat(inputdir,caseid(i),'_data/','delta_tavg_',caseid(i),'.dat');
f = readmatrix(fname,'ConsecutiveDelimitersRule','join'); 
plot(f(:,1),f(:,6)./del_the,'Marker',all_marks{i},'MarkerFaceColor',col{i},'MarkerEdgecolor',col1{i},'DisplayName',caseid_name(i));hold on
end
%plot(kris_data(:,1),kris_data(:,2),'*k','Displayname','DNS');hold on
%plot(langley_data(:,1),langley_data(:,2),'--k','Displayname','Langley fit');
legend show;legend('Location','Northeast','box','off','Fontsize',12,'Fontname','Timesnewroman');
xticks([0.2,0.4,0.8,1.2,1.6,2.0]);
ax=gca;ax.YAxis.FontSize = 13;ax.XAxis.FontSize = 13;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
xlab=xlabel('\boldmath$M_c$','Interpreter','Latex');set(xlab,'Fontsize',17);
ylab=ylabel('\boldmath$\dot{\delta}_{\theta}/\dot{\delta}_{\theta,inc}$', 'Interpreter','latex');set(ylab,'Fontsize',17);
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Linestyle'),'Linestyle','none')
set(findall(gcf,'-property','Fontname'),'Fontname','Timesnewroman')
set(findall(gcf,'-property','Markersize'),'Markersize',8)
screen2jpeg(strcat(outputdir1,'mom_rate.jpeg'));
%%%%%---------------------------------------------------------------------------------------------------------------------------
figure(2),clf
mc=["02","04","08","12","16","20"];
Mc=[0.2,0.4,0.8,1.2,1.6,2.0];
du=2*Mc*sqrt(1.4);
%t_i=[400,280,150,160,140,110];t_f=[900,600,350,340,280,250];  %Lad
%t_i=[560,380,200,150,110,100];t_f=[1000,800,400,330,270,240];  %AMD
t_i=[560,480,200,150,110,130];t_f=[1000,800,400,330,270,280];   %sigma
%t_i=[460,380,150,120,110,110];t_f=[900,700,350,300,250,250];  %mgm

col={'r',myg,'b',myv,myy,'k'};
all_marks = {'o','v','d','^','s','>','x','v','>','<','p','h'};
caseid="sigma_dyn";

for i=1:size(mc,2)
inputdir1=strcat('/home/vsj/Codes/runs180/run',mc(i),'/',caseid,'/postprocess/');
fname = strcat(inputdir1,'delta_',sprintf('%1.1f',Mc(i)),'_',caseid,'.dat');
f=readmatrix(fname);
idx1=find(f(:,1)==t_i(i));idx2=find(f(:,1)==t_f(i));
tau=(f(1:idx2,1)*du(i))/delta_0;
del=f(1:idx2,2)/delta_0;
plot(tau,del,'Marker',all_marks{i},'color',col{i},'LineWidth',1.5,'DisplayName',strcat('M_c=',sprintf('%1.1f',Mc(i))));
hold on

%tau=(f(idx1:idx2,1)*du(i))/delta_0;
%del=f(idx1:idx2,2)/delta_0;
%eq=polyfit(tau,del,1);
%delta1 = polyval(eq,tau);
%fittxt=strcat('Linearfit:y=',sprintf('%1.3f',eq(1,1)),'x+(',sprintf('%1.3f',eq(1,2)),')');
%plot(tau,delta1,'--k','LineWidth',1.5,'HandleVisibility','off');hold on
txt1 = '|';
t = text((f(idx1,1)*du(i))/delta_0,f(idx1,2)/delta_0,txt1);
t = text((f(idx2,1)*du(i))/delta_0,f(idx2,2)/delta_0,txt1);
hold on
end
ylim([1 9.5])
legend show;legend('Location','Northwest','box','off');
xlab=xlabel('\boldmath$(t\Delta{u})/{\delta}_{\theta{,0}}$','Interpreter','latex');set(xlab,'Fontsize',15);
ylab=ylabel('\boldmath$\delta_{\theta}/{\delta}_{\theta{,0}}$','Interpreter','latex');set(ylab,'Fontsize',15);
grid on;hAx=gca; set(hAx,'xminorgrid','off','yminorgrid','off')
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',3)
screen2jpeg(strcat(outputdir1,'mom_thick_',caseid,'.jpeg'));

