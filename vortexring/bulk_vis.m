clear all,clc
inputdir = "/home/vsj/Codes/vortexring/pr10_285mm_512/";
outputdir = "/home/vsj/Codes/vortexring/pr10_285mm_512/postprocess/bulk/";
%inputdir = "/home/vsj/Codes/vortexring/pr03_165mm_512/";
%outputdir = "/home/vsj/Codes/vortexring/pr03_165mm_512/postprocess/bulk/";
mkdir(outputdir)
 
%Read coordinates from *_coords.h5 file
fname_coords=strcat(inputdir,'vortexring_coords.h5');
x=h5read(fname_coords,'//X');y=h5read(fname_coords,'//Y');z=h5read(fname_coords,'//Z');
Lx=x(end,1,1);  Ly=2*y(1,end,1);    Lz=2*z(1,1,end);
[nx,ny,nz] = size(x);
dxi = Lx/(nx-1);  deta = Ly/(ny-1);  dzeta= Lz/(nz-1); 
x1 = 0;     xn = Lx;
y1 = -Ly/2; yn = Ly/2;
z1 = -Lz/2; zn = Lz/2;

counter = 10;
for i=1:size(counter,2)
file=strcat(inputdir,'vortexring_',sprintf('%04d',counter(i)),'.h5');
bulk=h5read(file,strcat('//bulk'));
time = round(h5readatt(file,'/','Time')*10^6) 

u_draw=squeeze(bulk(:,:,round(nz/2)));

n=256;
t1=linspace(min(min(u_draw)),max(max(u_draw)),n+1);
t1=t1(1:n);

figure(1),clf
contourf(x(:,1,1)', y(1,:,1), u_draw','LineColor','None','LevelList',t1);hold on
ax=gca;ax.YAxis.FontSize = 12;ax.XAxis.FontSize = 12;
xlab=xlabel('\boldmath$x(m)$');set(xlab,'Fontsize',15,'interpreter','latex');
ylab=ylabel('\boldmath$y(m)$');set(ylab,'Fontsize',15,'interpreter','latex');
title(strcat('\boldmath$\beta^{*}',':t=',sprintf('%04d',round((time))),'\mu{s}$'),'Fontsize',15,'interpreter','latex');
xmin = min(u_draw,[],'all')
xmax = max(u_draw,[],'all')
c=colorbar;c.FontSize=10;
%caxis([0 0.2])
xlim([0 0.6]);ylim([-0.4 0.4])
%daspect([0.7 0.7 0.7]); % Set data aspect ratio
colormap(jet(n));
file_png= strcat(outputdir,'bulk_contxy_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)

end
