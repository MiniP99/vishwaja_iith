
%%============ Code for finding kernals ====================%%

clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
mainid = "runs180";
runid  = "run20";
caseid = "amd_dyn";

c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;
delta_0=1;      %delta_0=Re_theta/(rho_0*Re*du);

inputdir =strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid,'/');
outputdir=strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid,'/postprocess_new/');
mkdir(outputdir);

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
dx = Lx/nx; dy = Ly/ny; dz = Lz/nz; 

%%-------------Start loop-------------------------------------------
for i=20:48
k=i-19;
counter=i;
time=counter*5;
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');

%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');
T=h5read(file,'//T');
t11sgs=h5read(file,'//t11');t12sgs=h5read(file,'//t12');t13sgs=h5read(file,'//t13');
q1sgs=h5read(file,'//q1');q2sgs=h5read(file,'//q2');q3sgs=h5read(file,'//q3');

%calculating derivatives-----------------------------
dudx=ddx_compact(u,dx,nx,ny,nz);dudy=ddy_compact(u,dy,nx,ny,nz);dudz=ddz_compact(u,dz,nx,ny,nz);
dvdx=ddx_compact(v,dx,nx,ny,nz);dvdy=ddy_compact(v,dy,nx,ny,nz);dvdz=ddz_compact(v,dz,nx,ny,nz);
dwdx=ddx_compact(w,dx,nx,ny,nz);dwdy=ddy_compact(w,dy,nx,ny,nz);dwdz=ddz_compact(w,dz,nx,ny,nz);
dTdx=ddx_compact(T,dx,nx,ny,nz);dTdy=ddy_compact(T,dy,nx,ny,nz);dTdz=ddz_compact(T,dz,nx,ny,nz);

S11 = dudx; S22 = dvdy ; S33=dwdz;
S12 = 0.5*(dudy+dvdx);
S13 = 0.5*(dudz+dwdx);
S23 = 0.5*(dvdz+dwdy);

num1 = (dudx.*(dudx.*S11 + dvdx.*S12 + dwdx.*S13) + dvdx.*(dudx.*S12 + dvdx.*S22 + dwdx.*S23) + dwdx.*(dudx.*S13 + dvdx.*S23 + dwdx.*S33)).*(dx^2);
num2 = (dudy.*(dudy.*S11 + dvdy.*S12 + dwdy.*S13) + dvdy.*(dudy.*S12 + dvdy.*S22 + dwdy.*S23) + dwdy.*(dudy.*S13 + dvdy.*S23 + dwdy.*S33)).*(dy^2);
num3 = (dudz.*(dudz.*S11 + dvdz.*S12 + dwdz.*S13) + dvdz.*(dudz.*S12 + dvdz.*S22 + dwdz.*S23) + dwdz.*(dudz.*S13 + dvdz.*S23 + dwdz.*S33)).*(dz^2);
numer = -(num1+num2+num3);

denom = (dudx.^2 + dudy.^2 + dudz.^2) + (dvdx.^2 + dvdy.^2 + dvdz.^2) + (dwdx.^2 + dwdy.^2 + dwdz.^2) + (1e-32);

nu_sgs=max(numer./denom,0);
nu_sgs_avg(k,:) = avgxz(nu_sgs,nx,nz);

n1 = (dudx.*dTdx.*dTdx + dvdx.*dTdx.*dTdy + dwdx.*dTdx.*dTdz).*(dx^2) ;
n2 = (dudy.*dTdy.*dTdx + dvdy.*dTdy.*dTdy + dwdy.*dTdy.*dTdz).*(dy^2) ;
n3 = (dudz.*dTdz.*dTdx + dvdz.*dTdz.*dTdy + dwdz.*dTdz.*dTdz).*(dz^2) ;
numer2 = -(n1+n2+n3);

denom2 = (dTdx.^2 + dTdy.^2 + dTdz.^2) + (1e-32);

kappa_sgs=max(numer2./denom2,0);
kappa_sgs_avg(k,:) = avgxz(kappa_sgs,nx,nz);
end

nu_tavg = mean(nu_sgs_avg);kappa_tavg = mean(kappa_sgs_avg);

nu = mean(nu_tavg);kappa = mean(kappa_tavg);

fname = strcat('/home/vsj/Codes/kernalsgs_',caseid,'_',sprintf('%1.1f',Mc),'.dat');
fileID = fopen(fname,'w');
for i=1:ny
fprintf(fileID,'%1.1f %2.8f %2.8f %2.8f \r\n',Mc,y(i),nu_tavg(i),kappa_tavg(i));
end
fclose(fileID);

fname = strcat('/home/vsj/Codes/kernalsgs_',caseid,'.dat');
fileID = fopen(fname,'a');
fprintf(fileID,'%1.1f %2.8f %2.8f\r\n',Mc,nu,kappa);
fclose(fileID);

plot(y,nu_tavg);
screen2jpeg('nusgs.png');
