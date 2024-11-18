clear all, clf

inputdir = "/home/vsj/Codes/vortexring/";
outputdir = "/home/vsj/Codes/vortexring/axial_plot/";
mkdir(outputdir)

%caseid = ["ctest","cvltest"];
%caseidname = ["cgrid.F90","cvlgrid.F90"];

caseid = ["pr10_285mm_128","pr10_285mm_256","pr10_285mm_512"];
caseidname = ["128\times128\times128","256\times256\times256","512\times512\times512"];
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={myg,'r','b',myy,'c'};
all_marks = {'s','o','+','v','p'};
lin_sty = {'-','--','-.'};

figure(1),clf
faxial = strcat(inputdir,'axial_pr10_285_520.csv');
axial_vel = readmatrix(faxial,'ConsecutiveDelimiters','join');
txt=strcat('$PIV (t=520\mu{s})$');
plot(axial_vel(:,1),axial_vel(:,2),'ok','LineWidth',1,'DisplayName',txt,'Markerfacecolor','k','Markersize',5);hold on

counter = 520;                  %50000-pr10_285mm  %40203-pr03_165mm(25)  %20102-pr03_165mm(50)

for j=1:size(counter,2)
%figure(1),clf
for i=1:size(caseid,2)
fname_coords=strcat(inputdir,caseid(i),'/vortexring_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
xstr=X(:,1,1)';ystr=Y(1,:,1);zstr=squeeze(Z(1,1,:))';
Lx=xstr(end);  Ly=2*ystr(end);    Lz=2*zstr(end);
nx=size(xstr,2);    ny=size(ystr,2);   nz=size(zstr,2);

file=strcat(inputdir,caseid(i),'/vortexring_',sprintf('%04d',counter(j)),'.h5');
time = round(h5readatt(file,'/','Time')*10^6) 
u=h5read(file,strcat('//u'));
u_y0 = u(:,round(ny/2),round(nz/2));
txt=strcat('$',caseidname(i),'$');
plot(xstr,u_y0,'color',col{i},'LineWidth',2, 'Linestyle','None','Displayname',txt,'Markersize',5,'LineStyle',lin_sty{i});hold on
xlim([0 0.17]);
end
ax=gca;ax.YAxis.FontSize = 18;ax.XAxis.FontSize = 18;
set(ax, 'TickLabelInterpreter', 'latex');
legend show;leg=legend('Location','Northeast','box','off','Fontsize',17,'interpreter','latex');
leg.ItemTokenSize = [25 15];
xlab=xlabel('\boldmath$x\,(m)$');set(xlab,'Fontsize',20,'interpreter','latex');
ylab=ylabel('\boldmath$u\,(m/s)$');set(ylab,'Fontsize',20,'interpreter','latex');
ylim([100 750]);yticks([150 300 450 600 750])
file_png= strcat(outputdir,'axial_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)
end
