%%============ Code for finding dilatation and vorticity  ====================%%

clear,clf    %%% mainid,runid,caseid,mc,time_step,time,file_video
mainid = "runs180";
runid  = "run20";
caseid = "mgm_dyn";

c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;
delta_0=1;      %delta_0=Re_theta/(rho_0*Re*du);

inputdir =strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid,'/');
outputdir=strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid,'/postprocess/');
outputdir1=strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid,'/contours/');
mkdir(outputdir1);

%file_video='/home/vsj/Codes/runs180/run20/mgm_dyn/contours/dilatation.mp4';
%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
dx = Lx/nx; dy = Ly/ny; dz = Lz/nz;

fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
delta1= readmatrix(fname);t1=delta1(:,1);d1=delta1(:,2);d2=delta1(:,4);
%t1=time,t2=mom_thick,t3=mom_thick_rate,t4=vor_thick,t5=ratio,t6=del_99

%obj=VideoWriter(file_video);
%obj.Quality=100;
%obj.FrameRate=2;
%open(obj);
%delta_0=Re_theta/(rho_0*Re*du);
delta_0=1;

for i=0:5:50
counter=i;
time=counter*5;
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');

%%------------- Read fields -----------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');

%---------------------------calculating derivatives-----------------------------
dudx=ddx_compact(u,Lx,nx,ny,nz);dudy=ddy_compact(u,Ly,nx,ny,nz);dudz=ddz_compact(u,Lz,nx,ny,nz);
dvdx=ddx_compact(v,Lx,nx,ny,nz);dvdy=ddy_compact(v,Ly,nx,ny,nz);dvdz=ddz_compact(v,Lz,nx,ny,nz);
dwdx=ddx_compact(w,Lx,nx,ny,nz);dwdy=ddy_compact(w,Ly,nx,ny,nz);dwdz=ddz_compact(w,Lz,nx,ny,nz);

%--------------------------Calculating dilatation,vorticity-----------------------------------
dil     = dudx+dvdy+dwdz;
ome_x   = 0.5*(dwdy-dvdz);ome_y=0.5*(dudz-dwdx);ome_z=0.5*(dvdx-dudy);
ome_mag = sqrt(ome_x.^2+ome_y.^2+ome_z.^2);

idx=find(t1==time);   
del=d1(idx)  %To find momentum thickness at particular time

dil = (dil*del)/du;  ome_x = (ome_x*del)/du; ome_y = (ome_y*del)/du; ome_z = (ome_z*del)/du; ome_mag = (ome_mag*del)/du;
%--------------------------Normalize components and plot------------------------------------------------
figure(1),clf
[Y,X,Z]=meshgrid(y,x,z);
hold on
slice(Y,X,Z,dil,min(min(min(Y))),min(min(min(X))),min(min(min(Z))))
slice(Y,X,Z,dil,max(max(max(Y))),max(max(max(X))),max(max(max(Z))))
view(2);
grid on;
colormap(gray(256));
colorbar('vertical');
shading interp;
caxis([-1 1]);
camup([1 0 0])
%contour(y, x, squeeze(ome_mag(:,:,nz/2)),[0.01 0.01],'Showtext','off','LineColor','r','Linewidth',1);
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
file_jpeg=strcat(outputdir1,'dil_new_','contxy_',sprintf('%04d',time),'_',sprintf('%1.1f',Mc),'_',caseid,'.jpeg');
screen2jpeg(file_jpeg);
end
