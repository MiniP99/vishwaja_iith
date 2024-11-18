clear all,clc
%inputdir = "/home/vsj/Codes/vortexring/pr10_285mm_512/";
%outputdir = "/home/vsj/Codes/vortexring/pr10_285mm_512/postprocess/gradrho_yz/";
inputdir = "/home/vsj/Codes/vortexring/pr03_165mm_512/";
outputdir = "/home/vsj/Codes/vortexring/pr03_165mm_512/postprocess/gradrho_yz/";
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
for i=1:nx
    for j=1:ny
        for k=1:nz
            xi(i,j,k)    = x1 + ((i-1)*dxi); 
            eta(i,j,k)   = y1 + ((j-1)*deta); 
            zeta(i,j,k)  = z1 + ((k-1)*dzeta); 
        end
    end
end

%Find metric multipliers to find derivatives in physical space
dxidxij = get_metrics(x,y,z,dxi,deta,dzeta);

counter = 16;
for i=1:size(counter,2)
file=strcat(inputdir,'vortexring_',sprintf('%04d',counter(i)),'.h5');
rho=h5read(file,strcat('//rho'));
u=h5read(file,strcat('//u'));
v=h5read(file,strcat('//v'));
w=h5read(file,strcat('//w'));
time = round(h5readatt(file,'/','Time')*10^6) 

drho_dxi   = ddxi_compact10_gen(rho,dxi);
drho_deta  = ddeta_compact10_sym(rho,deta);
drho_dzeta = ddzeta_compact10_sym(rho,dzeta);
drho_dx = dxidxij(:,:,:,1).*drho_dxi + dxidxij(:,:,:,4).*drho_deta + dxidxij(:,:,:,7).*drho_dzeta;
drho_dy = dxidxij(:,:,:,2).*drho_dxi + dxidxij(:,:,:,5).*drho_deta + dxidxij(:,:,:,8).*drho_dzeta;
drho_dz = dxidxij(:,:,:,3).*drho_dxi + dxidxij(:,:,:,6).*drho_deta + dxidxij(:,:,:,9).*drho_dzeta;
grad_rho_mag = sqrt( drho_dx.^2 + drho_dy.^2 + drho_dz.^2 );

dudxi   = ddxi_compact10_gen(u,dxi);
dudeta  = ddeta_compact10_sym(u,deta);
dudzeta = ddzeta_compact10_sym(u,dzeta);
dudx = dxidxij(:,:,:,1).*dudxi + dxidxij(:,:,:,4).*dudeta + dxidxij(:,:,:,7).*dudzeta;
dudy = dxidxij(:,:,:,2).*dudxi + dxidxij(:,:,:,5).*dudeta + dxidxij(:,:,:,8).*dudzeta;
dudz = dxidxij(:,:,:,3).*dudxi + dxidxij(:,:,:,6).*dudeta + dxidxij(:,:,:,9).*dudzeta;

dvdxi   = ddxi_compact10_gen(v,dxi);
dvdeta  = ddeta_compact10_sym(v,deta);
dvdzeta = ddzeta_compact10_sym(v,dzeta);
dvdx = dxidxij(:,:,:,1).*dvdxi + dxidxij(:,:,:,4).*dvdeta + dxidxij(:,:,:,7).*dvdzeta;
dvdy = dxidxij(:,:,:,2).*dvdxi + dxidxij(:,:,:,5).*dvdeta + dxidxij(:,:,:,8).*dvdzeta;
dvdz = dxidxij(:,:,:,3).*dvdxi + dxidxij(:,:,:,6).*dvdeta + dxidxij(:,:,:,9).*dvdzeta;

dwdxi   = ddxi_compact10_gen(w,dxi);
dwdeta  = ddeta_compact10_sym(w,deta);
dwdzeta = ddzeta_compact10_sym(w,dzeta);
dwdx = dxidxij(:,:,:,1).*dwdxi + dxidxij(:,:,:,4).*dwdeta + dxidxij(:,:,:,7).*dwdzeta;
dwdy = dxidxij(:,:,:,2).*dwdxi + dxidxij(:,:,:,5).*dwdeta + dxidxij(:,:,:,8).*dwdzeta;
dwdz = dxidxij(:,:,:,3).*dwdxi + dxidxij(:,:,:,6).*dwdeta + dxidxij(:,:,:,9).*dwdzeta;

ome_x=  0.5*(dwdy-dvdz);
ome_y=  0.5*(dudz-dwdx);
ome_z=  0.5*(dvdx-dudy);

mag_omega = sqrt(ome_x.^2 + ome_y.^2 + ome_z.^2);
u_draw_xy =  squeeze(mag_omega(:,:,round(nz/2)));
%xmax1 = max(u_draw_xy(:,152:end),[],'all'); %pr10_285_256
xmax1 = max(u_draw_xy,[],'all')
[idx1,idx2]=find(abs(u_draw_xy-xmax1)<10^-6);
u_draw = squeeze(grad_rho_mag(idx1(1),:,:));

%idx1 = 76; idx2 = 313; idx3 = 77; idx4 = 198; %512x512x512 (uni,noTI), PR=03, DRL=285mm
%u_draw = squeeze(mag_omega(idx1,:,:));

n=256;
t1=linspace(min(min(u_draw)),max(max(u_draw)),n+1);
t1=t1(1:n);

figure(1),clf
xstr=x(:,1,1)';ystr=y(1,:,1);zstr=squeeze(z(1,1,:))';
contourf(zstr, ystr, u_draw,'LineColor','None','LevelList',t1);hold on
ax=gca;ax.YAxis.FontSize = 12;ax.XAxis.FontSize = 12;
xlab=xlabel('\boldmath$z(m)$');set(xlab,'Fontsize',15,'interpreter','latex');
ylab=ylabel('\boldmath$y(m)$');set(ylab,'Fontsize',15,'interpreter','latex');
title(strcat('\boldmath$\nabla{\rho}(x=',sprintf('%1.3f',xstr(idx1(1))),'m):t=',sprintf('%04d',round((time))),'\mu{s}$'),'Fontsize',15,'interpreter','latex');
xmin = min(u_draw,[],'all')
xmax = max(u_draw,[],'all')
%caxis([0 40000])
c=colorbar;c.FontSize=10;%c.Ticks=[0 10000 20000];
colormap(flipud(gray(n)));
caxis([0 30])   % omega_z  %pr03_165mm
xlim([-0.1 0.1]);ylim([-0.1 0.1])
file_png= strcat(outputdir,'gradrho_contyz_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)

end
