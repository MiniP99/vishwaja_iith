clear all, clf,clc

inputdir = "/home/vsj/Codes/vortexring/";
outputdir = "/home/vsj/Codes/vortexring/radial_plot/";
mkdir(outputdir)

caseid = ["pr10_285mm_128","pr10_285mm_256","pr10_285mm_512"];
caseidname = ["128\times128\times128","256\times256\times256","512\times512\times512"];
%idx1 = 20;idx2=78; idx3 = 20;idx4=49;  %256x256x256 (uni,noTI), PR=10, DRL=285mm
%idx1 = 37;idx2=156 ; idx3 = 38;idx4=99 ;  %256x256x256 (uni,noTI), PR=10, DRL=285mm
%idx1 = 76; idx2 = 313; idx3 = 77; idx4 = 198; %512x512x512 (uni,noTI), PR=10, DRL=285mm
idx2_mat = [78,156,313];idx4_mat = [49,99,198];
%caseid = ["ctest","cvltest"];
%caseidname = ["cgrid.F90","cvlgrid.F90"];

myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={myg,'r','b',myy,'c'};
all_marks = {'s','o','+','v','p'};
lin_sty = {'-','--','-.'};

figure(1),clf
frad = strcat(inputdir,'radial_pr10_285_520.csv');
radial_vel = readmatrix(frad,'ConsecutiveDelimiters','join');
txt=strcat('$PIV (t=520\mu{s})$');
plot(radial_vel(:,1)-0.0072,radial_vel(:,2),'ok','LineWidth',1,'DisplayName',txt,'Markerfacecolor','k','Markersize',5);hold on


counter = 520;                
for m=1:size(counter,2)
for p=1:size(caseid,2)
fname_coords=strcat(inputdir,caseid(p),'/vortexring_coords.h5');
x=h5read(fname_coords,'//X');y=h5read(fname_coords,'//Y');z=h5read(fname_coords,'//Z');
xstr=x(:,1,1)';ystr=y(1,:,1);zstr=squeeze(z(1,1,:))';
Lx=x(end,1,1);  Ly=2*y(1,end,1);    Lz=2*z(1,1,end);
[nx,ny,nz] = size(x);
dxi = Lx/(nx-1);  deta = Ly/(ny-1);  dzeta= Lz/(nz-1); 
x1 = 0;     xn = Lx;
y1 = -Ly/2; yn = Ly/2;
z1 = -Lz/2; zn = Lz/2;

file=strcat(inputdir,caseid(p),'/vortexring_',sprintf('%04d',counter(m)),'.h5');
time = round(h5readatt(file,'/','Time')*10^6) 
u=h5read(file,strcat('//u'));  v=h5read(file,strcat('//v'));   w=h5read(file,strcat('//w'));
idx2 = idx2_mat(p);
idx4 = idx4_mat(p);

u_y1 = v(:,idx2,round(nz/2));
u_y2 = -v(:,idx4,round(nz/2));
u_y0 = (u_y1+u_y2)/2;

txt=strcat('$',caseidname(p),'$');
plot(xstr,u_y0,'LineStyle',lin_sty{p},'color',col{p},'LineWidth',2,'Displayname',txt,'Markersize',5);hold on
xlim([0 0.17]);
ax=gca;ax.YAxis.FontSize = 18;ax.XAxis.FontSize = 18;
set(ax, 'TickLabelInterpreter', 'latex');
legend show;leg=legend('Location','Southeast','box','off','Fontsize',17,'interpreter','latex');
leg.ItemTokenSize = [25 15];
xlab=xlabel('\boldmath$x\,(m)$');set(xlab,'Fontsize',20,'interpreter','latex');
ylab=ylabel('\boldmath$v\,(m/s)$');set(ylab,'Fontsize',20,'interpreter','latex');
%title(strcat('\boldmath$v(y=',sprintf('%1.3f',ystr(idx2(1))),'m):t=',sprintf('%04d',round((time))),'\mu{s}$'),'Fontsize',15,'interpreter','latex');
end

file_png= strcat(outputdir,'radial_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)
end
