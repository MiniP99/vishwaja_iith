clear all,clc
inputdir = "/home/vsj/Codes/vortexring/pr10_285mm_256/";
outputdir = "/home/vsj/Codes/vortexring/pr10_285mm_256/postprocess/iso_ome/";
%inputdir = "/home/vsj/Codes/vortexring/pr03_165mm_256/";
%outputdir = "/home/vsj/Codes/vortexring/pr03_165mm_256/postprocess/iso_ome/";
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

counter = 6;
for i=1:size(counter,2)
file=strcat(inputdir,'vortexring_',sprintf('%04d',counter(i)),'.h5');
rho=h5read(file,strcat('//rho'));
u=h5read(file,strcat('//u'));
v=h5read(file,strcat('//v'));
w=h5read(file,strcat('//w'));
time = round(h5readatt(file,'/','Time')*10^6) 

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
xmin = min(mag_omega,[],'all')
xmax = max(mag_omega,[],'all')

figure(1),clf
xstr=x(:,1,1)';ystr=y(1,:,1);zstr=squeeze(z(1,1,:))';
[x_new y_new z_new]=meshgrid(ystr,xstr,zstr);
p=patch(isosurface(x_new,y_new,z_new,mag_omega, 30000)) ;
set(p,'FaceColor','red','EdgeColor','none');
%view([134.5,17.14])
view([110,30])
camlight 
lighting gouraud
xlab=xlabel('\boldmath$y(m)$');set(xlab,'Fontsize',20,'interpreter','latex');
ylab=ylabel('\boldmath$x(m)$');set(ylab,'Fontsize',20,'interpreter','latex');
zlab=zlabel('\boldmath$z(m)$');set(zlab,'Fontsize',20,'interpreter','latex');
ax=gca;ax.YAxis.FontSize = 18;ax.XAxis.FontSize = 18;ax.ZAxis.FontSize = 18;
set(ax, 'TickLabelInterpreter', 'latex');
%tit=title(strcat('\boldmath$Iso surface:|\omega|=30000 s^{-1}, time=',sprintf('%04d',round((time))),'\mu{s}$'));set(tit,'Fontsize',20,'interpreter','latex');
ylim([xstr(3) 0.2]);xlim([-0.08 0.08]);zlim([-0.08 0.08]);
grid off;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
file_png= strcat(outputdir,'isosurface_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png);
end
