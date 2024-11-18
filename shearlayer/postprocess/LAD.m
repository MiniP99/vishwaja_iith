%%============ Code for finding artificial mu ====================%%

clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
mainid = "runs180";
runid  = "run20";
caseid = "amd_dyn";

c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;T_ref=1;
delta_0=1;      %delta_0=Re_theta/(rho_0*Re*du);

inputdir =strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',runid,'/',caseid,'/');
outputdir=strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',runid,'/',caseid,'/postprocess_new/');
mkdir(outputdir);

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
%dx = Lx/nx; dy = Ly/ny; dz = Lz/nz; 
dx = X(2,1,1)-X(1,1,1); dy = Y(1,2,1)-Y(1,1,1);dz = Z(1,1,2)-Z(1,1,1);

%%-------------Start loop-------------------------------------------
for i=20:48
k=i-19;
counter=i;
time=counter*5;
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');

%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');T=h5read(file,'//T');
mu=h5read(file,'//mu');
t11sgs=h5read(file,'//t11');t12sgs=h5read(file,'//t12');t13sgs=h5read(file,'//t13');
q1sgs=h5read(file,'//q1');q2sgs=h5read(file,'//q2');q3sgs=h5read(file,'//q3');

%calculating derivatives-----------------------------
dudx=ddx_compact(u,dx,nx,ny,nz);dudy=ddy_compact(u,dy,nx,ny,nz);dudz=ddz_compact(u,dz,nx,ny,nz);
dvdx=ddx_compact(v,dx,nx,ny,nz);dvdy=ddy_compact(v,dy,nx,ny,nz);dvdz=ddz_compact(v,dz,nx,ny,nz);
dwdx=ddx_compact(w,dx,nx,ny,nz);dwdy=ddy_compact(w,dy,nx,ny,nz);dwdz=ddz_compact(w,dz,nx,ny,nz);
dTdx=ddx_compact(T,dx,nx,ny,nz);dTdy=ddy_compact(T,dy,nx,ny,nz);dTdz=ddz_compact(T,dz,nx,ny,nz);

% mu physical
muphy_ref = (1/Re_theta);
T0        = 101325/(287*1.2);
Sconst    = 110.4/T0;
Skconst   = 194/T0;  %For Pr
muphy     = muphy_ref .*((T./T_ref).^1.5).*((T_ref+Sconst)./(T+Sconst));

% LAD
S11 = dudx; S22 = dvdy ; S33=dwdz;
S12 = 0.5*(dudy+dvdx);
S13 = 0.5*(dudz+dwdx);
S23 = 0.5*(dvdz+dwdy);

func = sqrt (S11.^2 + S22.^2 + S33.^2 + 2*(S12.^2 + S13.^2 + S23.^2)); %Magnitude of Strain rate tensor

der1 = ddx_compact(func,dx,nx,ny,nz);der2 = ddx_compact(der1,dx,nx,ny,nz);
der3 = ddx_compact(der2,dx,nx,ny,nz);der4 = ddx_compact(der3,dx,nx,ny,nz);
termx = der4 * (dx^6);

der1 = ddy_compact(func,dy,nx,ny,nz);der2 = ddy_compact(der1,dy,nx,ny,nz);
der3 = ddy_compact(der2,dy,nx,ny,nz);der4 = ddy_compact(der3,dy,nx,ny,nz);
termy = der4 * (dy^6);

der1 = ddz_compact(func,dz,nx,ny,nz);der2 = ddz_compact(der1,dz,nx,ny,nz);
der3 = ddz_compact(der2,dz,nx,ny,nz);der4 = ddz_compact(der3,dz,nx,ny,nz);
termz = der4 * (dz^6);

term = abs (termx + termy + termz);

cmu = 0.002 ; 
mustar = cmu .* rho .* term; 

% Filter twice
mustar = gaussian_filter(mustar,nx,ny,nz);
mustar = gaussian_filter(mustar,nx,ny,nz);

mustar_avg(k,:) = avgxz(mustar,nx,nz);
muphy_avg(k,:)  = avgxz(muphy,nx,nz);
mu_avg(k,:)  = avgxz(mu,nx,nz);

end

mustar_tavg = mean(mustar_avg); muphy_tavg = mean(muphy_avg); mu_tavg = mean(mu_avg);

mustar_one = mean(mustar_tavg); muphy_one = mean(muphy_tavg); mu_one = mean(mu_tavg);

%fname = strcat('/home/vsj/Codes/mustar_',caseid,'_',sprintf('%1.1f',Mc),'.dat');
%fileID = fopen(fname,'w');
%for i=1:ny
%fprintf(fileID,'%1.1f %2.8f %2.8f\r\n',Mc,y(i),mustar_tavg(i));
%end
%fclose(fileID);

fname = strcat('/home/vsj/Codes/mustar_',caseid,'.dat');
fileID = fopen(fname,'a');
fprintf(fileID,'%1.1f %2.8f %2.8f %2.8f %2.8f\r\n',Mc,mu_one,muphy_one,mustar_one,muphy_one+mustar_one);
fclose(fileID);

plot(y,mu_tavg,'r',y,muphy_tavg,'b',y,mustar_tavg,'g',y,mustar_tavg+muphy_tavg,'--r');
screen2jpeg('/home/vsj/Codes/mustar.png');
