clear
inputdir='/home/vsj/Codes/shearlayer_les/check/';
dir=strcat('/home/vsj/Codes/shearlayer_les/check/postprocess/');
dir1=strcat('/home/vsj/Codes/shearlayer_les/runs180/mgm_dyn_data/');
outputdir='/home/vsj/Codes/shearlayer_les/check/contour_lad/';
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

%file_video='/home/vsj/Codes/shearlayer_les/check/contours/uvel.mp4';
%obj=VideoWriter(file_video);
%obj.Quality=100;
%obj.FrameRate=6;
%open(obj);

for i=1:size(time,2)
counter=time(i);
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5')
u=h5read(file,strcat('//u'));

u_draw=squeeze(u(:,:,nz/2))';
n=256;
t1=linspace(min(min(u_draw)),max(max(u_draw)),n+1);
t1=t1(1:n);


figure(1),clf
contourf(x, y, u_draw,'LineColor','None','LevelList',t1);
ax=gca;ax.YAxis.FontSize = 13;ax.XAxis.FontSize = 13;
xlab=xlabel('\boldmath$x$');set(xlab,'Fontsize',17,'interpreter','latex');
ylab=ylabel('\boldmath$y$');set(ylab,'Fontsize',17,'interpreter','latex');
c=colorbar;c.FontSize=13;colormap(jet(n));%caxis([-3 3])
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
file_png= strcat(outputdir,var(1,j),'_contxy_',sprintf('%04d',time(i)),'_',sprintf('%1.1f',Mc),'_',caseid,'.png');
screen2jpeg(file_png)


%%-----create video-----
%f= imread(file_png);
%writeVideo(obj,f);
end
%obj.close();
