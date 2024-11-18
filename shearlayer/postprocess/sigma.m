%%============ Code for finding kernals ====================%%

clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
mainid = "runs180";
runid  = "run20";
caseid = "sigma_dyn";

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
%dx = Lx/nx; dy = Ly/ny; dz = Lz/nz; 
dx = X(2,1,1)-X(1,1,1); dy = Y(1,2,1)-Y(1,1,1);dz = Z(1,1,2)-Z(1,1,1);

%%-------------Start loop-------------------------------------------
for i=26:56
k=i-25;
counter=i;
time=counter*5;
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');

%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');

%calculating derivatives-----------------------------
dudx=ddx_compact(u,dx,nx,ny,nz);dudy=ddy_compact(u,dy,nx,ny,nz);dudz=ddz_compact(u,dz,nx,ny,nz);
dvdx=ddx_compact(v,dx,nx,ny,nz);dvdy=ddy_compact(v,dy,nx,ny,nz);dvdz=ddz_compact(v,dz,nx,ny,nz);
dwdx=ddx_compact(w,dx,nx,ny,nz);dwdy=ddy_compact(w,dy,nx,ny,nz);dwdz=ddz_compact(w,dz,nx,ny,nz);


G11 = dudx.^2 + dvdx.^2 + dwdx.^2;
G12 = dudx.*dudy + dvdx.*dvdy + dwdx.*dwdy;
G13 = dudx.*dudz + dvdx.*dvdz + dwdx.*dwdz;
G22 = dudy.^2 + dvdy.^2 + dwdy.^2;
G23 = dudy.*dudz + dvdy.*dvdz + dwdy.*dwdz;
G33 = dudz.^2 + dvdz.^2 + dwdz.^2;


I1   = G11 + G22 + G33;
I1sq = I1.*I1;
I1cu = I1sq.*I1;

I2 = -G11.*G11 - G22.*G22 - G33.*G33;
I2 = I2 - 2.*G12.*G12 - 2.*G13.*G13;
I2 = I2 - 2.*G23.*G23;
I2 = I2 + I1sq;
I2 = 0.5*I2;

I3 = G11.*(G22.*G33 - G23.*G23);
I3 = I3 + G12.*(G13.*G23 - G12.*G33);
I3 = I3 + G13.*(G12.*G23 - G22.*G13);

alpha1 = I1sq/9 - I2/3;
alpha1 = max(alpha1,0);

alpha2 = I1cu/27 - I1.*I2/6 + I3/2;
alpha1sqrt = sqrt(alpha1);
alpha1tmp = alpha1.*alpha1sqrt;
alpha1tmp = alpha2./(alpha1tmp + 1e-13);
alpha1tmp = min(alpha1tmp,1);
alpha1tmp = max(alpha1tmp,-1);
alpha1tmp = acos(alpha1tmp);
alpha3 = (1/3)*(alpha1tmp);

sigma1sq = I1/3 + 2.*alpha1sqrt.*cos(alpha3);
sigma1sq = max(sigma1sq,0);
sigma1 = sqrt(sigma1sq);

sigma2 = pi/3 + alpha3;
sigma2 = (-2).*alpha1sqrt.*cos(sigma2);
sigma2 = sigma2 + I1/3;
sigma2 = max(sigma2,0);
sigma2 = sqrt(sigma2);
            
sigma3 = pi/3 - alpha3;
sigma3 = (-2).*alpha1sqrt.*cos(sigma3);
sigma3 = sigma3 + I1/3;
sigma3 = max(sigma3,0);
sigma3 = sqrt(sigma3);

nusgs_sigma = sigma3.*(sigma1 - sigma2).*(sigma2 - sigma3)./(sigma1sq + 1e-15);
nu_sgs_avg(k,:) = avgxz(nusgs_sigma,nx,nz);
end

nu_tavg = mean(nu_sgs_avg);

nu = mean(nu_tavg);

fname = strcat('/home/vsj/Codes/kernalsgs_sigma_',caseid,'_',sprintf('%1.1f',Mc),'.dat');
fileID = fopen(fname,'w');
for i=1:ny
fprintf(fileID,'%1.1f %2.8f %2.8f \r\n',Mc,y(i),nu_tavg(i));
end
fclose(fileID);

fname = strcat('/home/vsj/Codes/kernalsgs_sigma_',caseid,'.dat');
fileID = fopen(fname,'a');
fprintf(fileID,'%1.1f %2.8f\r\n',Mc,nu);
fclose(fileID);

plot(y,nu_tavg);
screen2jpeg('nusgs.png');
