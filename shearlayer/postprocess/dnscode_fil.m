%%============ Code for finding momentum thickness, del_99, TKE Budgets ====================%%

clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
caseid = "dns";

inputdir =strcat('/home/vsj/dns_data_csl/mc2p0/my_fil/');
outputdir=strcat('/home/vsj/dns_data_csl/mc2p0/postprocess_fil/');
mkdir(outputdir);

%%-------------Start loop-------------------------------------------
for i=5:11
counter=i
time=counter*10;
file=strcat(inputdir,'sl_fil',sprintf('%04d',counter),'.h5');
%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');p=h5read(file,'//p');
mu=h5read(file,'//mu');
bulk=h5read(file,'//bulk');
[nx ny nz]=size(u);
%%------------- Read Coordinates x, y, z-----------------------------
Lx = 80; Ly = 80; Lz = 40;
%nx = 1024; ny = 1448; nz = 512;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;
delta_0=1;
t11sgs=zeros(nx,ny,nz);t12sgs=zeros(nx,ny,nz);t13sgs=zeros(nx,ny,nz);
t22sgs=zeros(nx,ny,nz);t23sgs=zeros(nx,ny,nz);t33sgs=zeros(nx,ny,nz);
%calculating derivatives-----------------------------
dudx=ddx_compact(u,dx,nx,ny,nz);dudy=ddy_compact(u,dy,nx,ny,nz);dudz=ddz_compact(u,dz,nx,ny,nz);
dvdx=ddx_compact(v,dx,nx,ny,nz);dvdy=ddy_compact(v,dy,nx,ny,nz);dvdz=ddz_compact(v,dz,nx,ny,nz);
dwdx=ddx_compact(w,dx,nx,ny,nz);dwdy=ddy_compact(w,dy,nx,ny,nz);dwdz=ddz_compact(w,dz,nx,ny,nz);
% calculating averages,fluctuations------------------
rho_bar=avgxz(rho,nx,nz);
u_bar=avgxz(u,nx,nz);[u_tilde,u_pp]=favre_avg_fluct(rho,u,nx,ny,nz);
v_bar=avgxz(v,nx,nz);[v_tilde,v_pp]=favre_avg_fluct(rho,v,nx,ny,nz);
w_bar=avgxz(w,nx,nz);[w_tilde,w_pp]=favre_avg_fluct(rho,w,nx,ny,nz);
% tausgs average
t11sgs_bar=avgxz(t11sgs,nx,nz);t12sgs_bar=avgxz(t12sgs,nx,nz);t13sgs_bar=avgxz(t13sgs,nx,nz);
t22sgs_bar=avgxz(t22sgs,nx,nz);t23sgs_bar=avgxz(t23sgs,nx,nz);t33sgs_bar=avgxz(t33sgs,nx,nz);
% calculating derivatives of favre_avg
du_tildx=zeros(1,ny);du_tildy=ddy1D_compact(u_tilde,dy,nx,ny,nz);du_tildz=zeros(1,ny);
dv_tildx=zeros(1,ny);dv_tildy=ddy1D_compact(v_tilde,dy,nx,ny,nz);dv_tildz=zeros(1,ny);
dw_tildx=zeros(1,ny);dw_tildy=ddy1D_compact(w_tilde,dy,nx,ny,nz);dw_tildz=zeros(1,ny);
% calculating  derivatives of favre_fluct--------------
du_ppdx=ddx_compact(u_pp,dx,nx,ny,nz);du_ppdy=ddy_compact(u_pp,dy,nx,ny,nz);du_ppdz=ddz_compact(u_pp,dz,nx,ny,nz);
dv_ppdx=ddx_compact(v_pp,dx,nx,ny,nz);dv_ppdy=ddy_compact(v_pp,dy,nx,ny,nz);dv_ppdz=ddz_compact(v_pp,dz,nx,ny,nz);
dw_ppdx=ddx_compact(w_pp,dx,nx,ny,nz);dw_ppdy=ddy_compact(w_pp,dy,nx,ny,nz);dw_ppdz=ddz_compact(w_pp,dz,nx,ny,nz);
%--------------------------------- calculating momentum thickness------------------------------------
[R11,R12,R13,R22,R23,R33]=reynold_stress(rho,u_pp,v_pp,w_pp,nx,nz);
R11 = rho_bar.*R11;R12 = rho_bar.*R12;R13 = rho_bar.*R13;R22 = rho_bar.*R22;R23 = rho_bar.*R23;R33 = rho_bar.*R33;

R11=R11-t11sgs_bar;R12=R12-t12sgs_bar;R13=R13-t13sgs_bar;R22=R22-t22sgs_bar;R23=R23-t23sgs_bar;R33=R33-t33sgs_bar;

func1=rho_bar.*(0.5*du-u_tilde).*(0.5*du+u_tilde);
delta_theta=(1/(rho_0*(du^2)))*sum(func1)*dy;

%-------------------------------- calculating momentum thickness rate----------------------------------
%func2=rho_bar.*R12.*du_tildy;
func2=R12.*du_tildy;         %%%rho_bar included in R12
delta_dot=(-2/(rho_0*(du^3)))*sum(func2)*dy;

%----------------------------------Calculating vorticity thickness-------------
delta_omega=du/max(abs(ddy1D_compact(u_bar,dy,nx,ny,nz)));
D_w=delta_omega/delta_theta;

%----------------------------------calculating delta_99-------------(Matsuno)
u1=0.99*du/2;u2=-0.99*du/2;
[val_1,idx_1]=min(abs(u_bar(1,:)-u1));[val_2,idx_2]=min(abs(u_bar(1,:)-u2));
delta_1=y(idx_1);delta_2=y(idx_2);
delta_99=abs(delta_1)+abs(delta_2);

%----------------------------------calculating Mach numbers-------------(Matsuno)
S = ddy1D_compact(u_tilde,dy,nx,ny,nz);
C = avgxz(sqrt(1.4*(p./rho)),nx,nz);

[value_yc,idx_yc] = max(abs(S)); 
y_c  = y(idx_yc);    %[val_yc,idx_yc]=min(abs(y_c(i)-y));
S_yc = S(idx_yc);
C_yc = C(idx_yc);

%%Turbulent Mach number---
num1 = sqrt(abs((R11+R22+R33)./rho_bar));  num1_yc=num1(idx_yc);
num2 = sqrt(abs(R22./rho_bar))          ;  num2_yc=num2(idx_yc);
num3 = sqrt(abs(R12./rho_bar))          ;  num3_yc=num3(idx_yc);
M_t  = num1_yc/C_yc;  M_tv  = num2_yc/C_yc;  M_tau = num3_yc/C_yc;

%% Gradient Mach Number--
[v_bar,v_rey_fluct]    = reynolds_avg_fluct(v,nx,ny,nz);
RHS = 0.1* avgxz(v_rey_fluct(:,idx_yc,:).*v_rey_fluct(:,idx_yc,:),nx,nz);
LHS = zeros(1,ny/2);
for k=1:idx_yc-1%:size(delta_y,2)
index1=idx_yc-k;index2=idx_yc+k;
if(index2<=ny)
LHS(k)=avgxz(v_rey_fluct(:,index1,:).*v_rey_fluct(:,index2,:),nx,nz);
end
end
[variable,index] = min(abs(RHS-LHS));
delta_y  = index*2*dy;
M_g      = S_yc*delta_y/C_yc;

%u_delta-------
y1 = y_c+delta_y/2; y2 = y_c-delta_y/2; [val_y1,idx_y1] = min(abs(y1-y)); [val_y2,idx_y2] = min(abs(y2-y));
U_del= u_tilde(idx_y1)-u_tilde(idx_y2);


%Integrated Reynold stresses
R11_int=(sum(R11)*dy)/(rho_0*(du^2)*delta_theta);R12_int=(sum(R12)*dy)/(rho_0*(du^2)*delta_theta);
R13_int=(sum(R13)*dy)/(rho_0*(du^2)*delta_theta);R22_int=(sum(R22)*dy)/(rho_0*(du^2)*delta_theta);
R23_int=(sum(R23)*dy)/(rho_0*(du^2)*delta_theta);R33_int=(sum(R33)*dy)/(rho_0*(du^2)*delta_theta);
%Production-----------------------------
Pij = -(R11.*du_tildx+R12.*du_tildy+R13.*du_tildz+R12.*dv_tildx+R22.*dv_tildy+R23.*dv_tildz+R13.*dw_tildx+R23.*dw_tildy+R33.*dw_tildz);
%Pij = rho_bar.*Pij;         %%tausgs included in Rij and also rho_bar

%P12 from above
P12 = R12.*du_tildy;           %%In K.Matsuno

%Dissipation------------------------------
[tau11,tau12,tau13,tau22,tau23,tau33]=shearstress(mu,bulk,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,t11sgs,t12sgs,t13sgs,t22sgs,t23sgs,t33sgs);               %%tausgs included
[tau11_bar,tau11_rey_fluct,tau12_bar,tau12_rey_fluct,tau13_bar,tau13_rey_fluct,tau22_bar,tau22_rey_fluct,tau23_bar,tau23_rey_fluct,tau33_bar,tau33_rey_fluct]=tau_avg_fluct(tau11,tau12,tau13,tau22,tau23,tau33,nx,ny,nz);

func_tau=tau11.*du_ppdx+tau12.*du_ppdy+tau13.*du_ppdz+tau12.*dv_ppdx+tau22.*dv_ppdy+tau23.*dv_ppdz+tau13.*dw_ppdx+tau23.*dw_ppdy+tau33.*dw_ppdz;
%func_tau=tau11_rey_fluct.*du_ppdx+tau12_rey_fluct.*du_ppdy+tau13_rey_fluct.*du_ppdz+tau12_rey_fluct.*dv_ppdx+tau22_rey_fluct.*dv_ppdy+tau23_rey_fluct.*dv_ppdz+tau13_rey_fluct.*dw_ppdx+tau23_rey_fluct.*dw_ppdy+tau33_rey_fluct.*dw_ppdz;

diss=avgxz(func_tau,nx,nz);


% Transport------------------------

[p_bar,p_reynold_fluct]    = reynolds_avg_fluct(p,nx,ny,nz);
[rho_bar,rho_reynold_fluct]= reynolds_avg_fluct(rho,nx,ny,nz);
[u_bar,u_reynold_fluct]    = reynolds_avg_fluct(u,nx,ny,nz);
[v_bar,v_reynold_fluct]    = reynolds_avg_fluct(v,nx,ny,nz);
[w_bar,w_reynold_fluct]    = reynolds_avg_fluct(w,nx,ny,nz);

term1=(avgxz(rho.*u_pp.*u_pp.*v_pp,nx,nz)+avgxz(rho.*v_pp.*v_pp.*v_pp,nx,nz)+avgxz(rho.*w_pp.*w_pp.*v_pp,nx,nz))/2;
term2=avgxz(p_reynold_fluct.*v_reynold_fluct,nx,nz);
term3=avgxz(tau12_rey_fluct.*u_pp,nx,nz)+avgxz(tau22_rey_fluct.*v_pp,nx,nz)+avgxz(tau23_rey_fluct.*w_pp,nx,nz);
func_trans=term1+term2-term3;
trans=-ddy1D_compact(func_trans,dy,nx,ny,nz);

%R11 = R11/rho_bar;R12 = R12/rho_bar;R13 = R13/rho_bar;R22 = R22/rho_bar;R23 = R23/rho_bar;R33 = R33/rho_bar;

% Baropycnal----------------
baro = - avgxz(v_pp,nx,nz) .* ddy1D_compact(p_bar,dy,nx,ny,nz);

%Pressure dilation-----------------------------
dil_func =  p_reynold_fluct .* (du_ppdx+dv_ppdy+dw_ppdz);
pre_dil = avgxz(dil_func,nx,ny);


%Pressure-strain terms----------------------
pi_11 = avgxz(p_reynold_fluct.*2.*du_ppdx,nx,ny);
pi_12 = avgxz(p_reynold_fluct.*(du_ppdy+dv_ppdx),nx,ny);
pi_22 = avgxz(p_reynold_fluct.*2.*dv_ppdy,nx,ny);

%mass flux-------------------------
mass_flux_ru = avgxz(rho_reynold_fluct.*u_reynold_fluct,nx,nz)/(rho_0*du);
mass_flux_rv = avgxz(rho_reynold_fluct.*v_reynold_fluct,nx,nz)/(rho_0*du);

%Anisotropy--------------------------------
Rkk = R11 + R22 + R33; KE = Rkk/2; K1 = (2/3)*KE;
KE_int = (sum(KE)*dy)/(rho_0*(du^2)*delta_theta);
%b11 = (R11-K1)./Rkk;
%b12 = R12./Rkk;
%b22 = (R22-K1)./Rkk;

% Reynolds number,kolmogorov's scale------
nu = mu./rho; 
%Re_theta_t = (du*delta_theta)./nu; 
Re_theta_t = (delta_theta)./nu;
Re_theta_t_max = max(Re_theta_t,[],'all');
Re_theta_t_min = min(Re_theta_t,[],'all');
Re_theta_t_avg = mean(Re_theta_t,'all');
eeta = min((diss.^(-1/4)).*(nu.^(3/4)),[],'all'); 

%------------------------------Write delta to a file-------------------------------
fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
fileID = fopen(fname,'a');   %%% U_del/du
fprintf(fileID,'%04d %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f\r\n',time,delta_theta,delta_dot,delta_omega,D_w,delta_99,delta_y,y_c,U_del/du);
fclose(fileID);
%----------------------------Write mach numbers to a file-------------------------
fname = strcat(outputdir,'mach_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
fileID = fopen(fname,'a');   
fprintf(fileID,'%04d %2.8f %2.8f %2.8f %2.8f \r\n',time,M_t,M_tv,M_tau,M_g);
fclose(fileID);
%----------------------------Write reynold numbers and kolmogorov's scale to a file-------------------------
fname = strcat(outputdir,'Rey_kol_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
fileID = fopen(fname,'a');   %%% eeta/dx
fprintf(fileID,'%04d %6.2f %6.2f %6.2f %2.8f \r\n',time,Re_theta_t_max,Re_theta_t_min,Re_theta_t_avg,eeta/dx);
fclose(fileID);

%%%---------------------------Non-dimensional terms-------------------------------

P12_int = (abs(sum(P12)*dy))/(du^3); diss_int = (sum(diss)*dy)/(du^3);
pi_11_int = (sum(pi_11)*dy)/(du^3) ; pi_12_int = (sum(pi_12)*dy)/(du^3) ; pi_22_int = (sum(pi_22)*dy)/(du^3);

u_bar=u_bar/du;u_tilde=u_tilde/du;v_bar=v_bar/du;v_tilde=v_tilde/du;
R11=sqrt(abs(R11))/du;R12=sqrt(abs(R12))/du;R13=sqrt(abs(R13))/du;R22=sqrt(abs(R22))/du;R23=sqrt(abs(R23))/du;R33=sqrt(abs(R33))/du;
Pij=(Pij*delta_theta)/(du^3);diss=(diss*delta_theta)/(du^3);trans=(trans*delta_theta)/(du^3);baro=(baro*delta_theta)/(du^3);pre_dil=(pre_dil*delta_theta)/(du^3);

fname=strcat(outputdir,'vel_tke_flux_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',time),'.dat');
fileID=fopen(fname,'w');
for j=1:ny
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n',u_bar(1,j),u_tilde(1,j),v_bar(1,j),v_tilde(1,j),R11(1,j),R12(1,j),R13(1,j),R22(1,j),R23(1,j),R33(1,j),Pij(1,j),diss(1,j),trans(1,j),baro(1,j),pre_dil(1,j),mass_flux_ru(1,j),mass_flux_rv(1,j));
end
fclose(fileID);

fname = strcat(outputdir,'Rij_P12_dis_pi_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
fileID = fopen(fname,'a');
fprintf(fileID,'%04d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\r\n',time,R11_int,R12_int,R13_int,R22_int,R23_int,R33_int,KE_int,P12_int,diss_int,pi_11_int,pi_12_int,pi_22_int);
fclose(fileID);                %% In non-dimensional form



end
