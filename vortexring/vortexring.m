clear all,clc
inputdir = "/home/vsj/Codes/vortexring/pr10_285mm_512/";
outputdir = "/home/vsj/Codes/vortexring/pr10_285mm_512/postprocess/grad_rho/";
%inputdir = "/home/vsj/Codes/vortexring/pr03_165mm_512/";
%outputdir = "/home/vsj/Codes/vortexring/pr03_165mm_512/postprocess/grad_rho/";
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

counter = 8;
for i=1:size(counter,2)
file=strcat(inputdir,'vortexring_',sprintf('%04d',counter(i)),'.h5');
rho=h5read(file,strcat('//rho'));
bulk=h5read(file,strcat('//bulk'));
time = round(h5readatt(file,'/','Time')*10^6) 

drho_dxi   = ddxi_compact10_gen(rho,dxi);
drho_deta  = ddeta_compact10_sym(rho,deta);
drho_dzeta = ddzeta_compact10_sym(rho,dzeta);
%dfdx = dxidx*dfdxi + detadx*dfdeta + dzetadx*dfdzeta
%dfdy = dxidy*dfdxi + detady*dfdeta + dzetady*dfdzeta
%dfdz = dxidz*dfdxi + detadz*dfdeta + dzetadz*dfdzeta
drho_dx = dxidxij(:,:,:,1).*drho_dxi + dxidxij(:,:,:,4).*drho_deta + dxidxij(:,:,:,7).*drho_dzeta;
drho_dy = dxidxij(:,:,:,2).*drho_dxi + dxidxij(:,:,:,5).*drho_deta + dxidxij(:,:,:,8).*drho_dzeta;
drho_dz = dxidxij(:,:,:,3).*drho_dxi + dxidxij(:,:,:,6).*drho_deta + dxidxij(:,:,:,9).*drho_dzeta;

grad_rho_mag = sqrt( drho_dx.^2 + drho_dy.^2 + drho_dz.^2 );
u_draw=squeeze(grad_rho_mag(:,:,round(nz/2)));

n=256;
t1=linspace(min(min(u_draw)),max(max(u_draw)),n+1);
t1=t1(1:n);

figure(1),clf
contourf(x(:,1,1)', y(1,:,1), u_draw','LineColor','None','LevelList',t1);hold on
ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
set(ax, 'TickLabelInterpreter', 'latex');
xlab=xlabel('\boldmath$x\,(m)$');set(xlab,'Fontsize',16,'interpreter','latex');
ylab=ylabel('\boldmath$y\,(m)$');set(ylab,'Fontsize',16,'interpreter','latex');
%title(strcat('\boldmath$\nabla{\rho}',':t=',sprintf('%04d',round((time))),'\mu{s}$'),'Fontsize',15,'interpreter','latex');
xmin = min(u_draw,[],'all')
xmax = max(u_draw,[],'all')
caxis([0 150])   % omega_z  %pr10_285
%caxis([0 40])   % omega_z  %pr03_165mm
c=colorbar;c.FontSize=15;
set(c,'TickLabelInterpreter','latex')
colormap(flipud(gray(n)));
%colormap(jet(n));
yticks([-0.1 0 0.1])
%yticks([-0.5 -0.2 0 0.2 0.5])
xlim([0 0.4]);ylim([-0.1 0.1])
daspect([0.7 0.7 0.7]); % Set data aspect ratio
%legend show;legend('Location','northeast','box','off','FontSize',10,'Numcolumns',2)
file_png= strcat(outputdir,'grad_rho_contxy_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)

end
