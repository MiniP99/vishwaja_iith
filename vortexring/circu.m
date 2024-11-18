clear all,clc
inputdir = "/home/vsj/Codes/vortexring/check/";
outputdir = "/home/vsj/Codes/vortexring/check/postprocess/circu/";
%inputdir = "/home/vsj/Codes/vortexring/pr03_165mm_512_TI10/";
%outputdir = "/home/vsj/Codes/vortexring/pr03_165mm_512_TI10/postprocess/circu/";
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


counter = 1;
for i=1:size(counter,2)
file=strcat(inputdir,'vortexring_',sprintf('%04d',counter(i)),'.h5');
rho=h5read(file,strcat('//rho'));
u=h5read(file,strcat('//u'));
v=h5read(file,strcat('//v'));
w=h5read(file,strcat('//w'));
time = round(h5readatt(file,'/','Time')*10^6) 

u_draw =  squeeze(u(:,:,round(nz/2)));
v_draw =  squeeze(v(:,:,round(nz/2)));
xstr=x(:,1,1)';ystr=y(1,:,1);zstr=squeeze(z(1,1,:))';
[X,Y]=meshgrid(xstr,ystr);
figure(1),clf
quiver(X,Y,u_draw',v_draw','r');
xlim([0 0.1]);ylim([0.05 0.1]);
ax=gca;ax.YAxis.FontSize = 12;ax.XAxis.FontSize = 12;
xlab=xlabel('\boldmath$x(m)$');set(xlab,'Fontsize',15,'interpreter','latex');
ylab=ylabel('\boldmath$y(m)$');set(ylab,'Fontsize',15,'interpreter','latex');
file_png= strcat(outputdir,'uv_contxy_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)

u_draw =  squeeze(u(1,:,:));
n=256;
t1=linspace(min(min(u_draw)),max(max(u_draw)),n+1);
t1=t1(1:n);
xstr=x(:,1,1)';ystr=y(1,:,1);zstr=squeeze(z(1,1,:))';
contourf(ystr, zstr, u_draw,'LineColor','None','LevelList',t1);hold on
ax=gca;ax.YAxis.FontSize = 12;ax.XAxis.FontSize = 12;
xlim([-0.05 0.05]);ylim([-0.05 0.05])
c=colorbar;c.FontSize=10;%c.Ticks=[0 10000 20000];
colormap(jet(n));
file_png= strcat(outputdir,'u_contyz_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)

figure(1),clf
plot(ystr,u_draw(:,round(nz/2)-1),'-or');
%plot(zstr,u_draw(round(ny/2)-1,:),'-or');
file_png= strcat(outputdir,'u_zline_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)





end
