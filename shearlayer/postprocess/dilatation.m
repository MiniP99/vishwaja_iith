%%============ Code for finding dilatation and vorticity  ====================%%
clear
inputdir='/home/vsj/Codes/shearlayer_les/check/';
dir=strcat('/home/vsj/Codes/shearlayer_les/check/postprocess/');
dir1=strcat('/home/vsj/Codes/shearlayer_les/runs180/mgm_dyn_data/');
outputdir='/home/vsj/Codes/shearlayer_les/check/contour_dil/';
mkdir(outputdir)
caseid="mgm_dyn";          %% change mc, domain

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);     Ly=2*y(end);  Lz=z(end);
nx=size(x,2);  ny=size(y,2); nz=size(z,2);
dx = Lx/nx;    dy= Ly/ny;    dz= Lz/nz; 
kx = fftshift((-nx/2:nx/2-1)*(2*pi/Lx));
ky = fftshift((-ny/2:ny/2-1)*(2*pi/Ly));
kz = fftshift((-nz/2:nz/2-1)*(2*pi/Lz));
[kx3D,ky3D,kz3D] = ndgrid(kx,ky,kz);

% inputs
c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;
%delta_0=Re_theta/(rho_0*Re*du);
delta_0=1;

clf
time=0:250;

file_video='/home/vsj/Codes/shearlayer_les/check/contour_dil/dil.mp4';
obj=VideoWriter(file_video);
obj.Quality=100;
obj.FrameRate=12;
open(obj);


fname = strcat(dir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
delta1= readmatrix(fname);t1=delta1(:,1);d1=delta1(:,2);d2=delta1(:,4);

for i=1:size(time,2)
counter=time(i)/1;
tau = time(i)*du/delta_0;

%%--------start loop-------------------------------------- 
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');

%%------------- Read fields -----------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');

%---------------------------calculating derivatives-----------------------------
dudx=ddx_compact(u,dx,nx,ny,nz);dudy=ddy_compact(u,dy,nx,ny,nz);
dudz=ddz_compact(u,dz,nx,ny,nz);
dvdx=ddx_compact(v,dx,nx,ny,nz);dvdy=ddy_compact(v,dy,nx,ny,nz);
dvdz=ddz_compact(v,dz,nx,ny,nz);
dwdx=ddx_compact(w,dx,nx,ny,nz);dwdy=ddy_compact(w,dy,nx,ny,nz);
dwdz=ddz_compact(w,dz,nx,ny,nz);

%--------------------------Calculating dilatation,vorticity-----------------------------------
dil     = dudx+dvdy+dwdz;
ome_x   = 0.5*(dwdy-dvdz);ome_y=0.5*(dudz-dwdx);ome_z=0.5*(dvdx-dudy);
ome_mag = sqrt(ome_x.^2+ome_y.^2+ome_z.^2);

idx=find(t1==time(i));   
del=d1(idx)  %To find momentum thickness at particular time

%--------------------------Normalize components and plot------------------------------------------------

dil = (dil*del)/du;  ome_x = (ome_x*del)/du; ome_y = (ome_y*del)/du; ome_z = (ome_z*del)/du; ome_mag = (ome_mag*del)/du;

dil_draw = squeeze(dil(:,:,nz/2))';
n=256;
lev=linspace(min(min(dil_draw)),max(max(dil_draw)),n+1);
lev=lev(1:n);

figure(1),clf
contourf(x, y, dil_draw,'LineColor','None','LevelList',lev)
hold on
contour(x, y, squeeze(ome_mag(:,:,nz/2))',[0.01 0.01],'Showtext','off','LineColor','r','Linewidth',1);
pbaspect([2 2 1])
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
tit=title(strcat('\boldmath$\tau=',sprintf('%4.0f',tau),'$'));set(tit,'Fontsize',20,'interpreter','latex');
xlab=xlabel('\boldmath$x/\delta_{\theta,0}$');set(xlab,'Fontsize',20,'interpreter','latex');
ylab=ylabel('\boldmath$y/\delta_{\theta,0}$');set(ylab,'Fontsize',20,'interpreter','latex');
yticks([-40 -20 0 20 40])
ylim([-35 35])
c=colorbar,c.FontSize=15;colormap(gray(256));caxis([-1 0.8]);

set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',1)
file_png= strcat(outputdir,'dil_contxy_',sprintf('%04d',time(i)),'_',sprintf('%1.1f',Mc),'_',caseid,'.png');
screen2jpeg(file_png)
%%-----create video-----
f= imread(file_png);
writeVideo(obj,f);
end
obj.close();
