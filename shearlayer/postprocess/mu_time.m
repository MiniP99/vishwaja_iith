%%============ Code for finding artificial mu ====================%%
clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
mainid = "runs180";
runid  = "run20";
caseid = "mgm_dyn";

sgsid = 0; %1=amd,2=sigma,3-Smag

c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;T_ref=1;
delta_0=1;      %delta_0=Re_theta/(rho_0*Re*du);

inputdir =strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',runid,'/',caseid,'/');
%outputdir=strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',runid,'/',caseid,'/postprocess_new/');
%mkdir(outputdir);

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
dx = Lx/nx; dy = Ly/ny; dz = Lz/nz; 
deltales = (dx*dy*dz)^(1/3);

t_i = 0; t_f = 250; time_step = 5;
time_total = t_i:time_step:t_f

diss_phy_int = zeros(size(time_total,2),1);
diss_star_int = zeros(size(time_total,2),1);
diss_sgs_int = zeros(size(time_total,2),1);

k = 1;
%%-------------Start loop-------------------------------------------
for i = t_i : time_step : t_f
disp(k)
counter = i/time_step
time = time_total(k)
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');

%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');T=h5read(file,'//T');
mu=h5read(file,'//mu');bulk=h5read(file,'//bulk');
t11sgs=h5read(file,'//t11');t12sgs=h5read(file,'//t12');t13sgs=h5read(file,'//t13');
t22sgs=h5read(file,'//t22');t23sgs=h5read(file,'//t23');t33sgs=h5read(file,'//t33');
%t11sgs=zeros(nx,ny,nz);t12sgs=zeros(nx,ny,nz);t13sgs=zeros(nx,ny,nz);
%t22sgs=zeros(nx,ny,nz);t23sgs=zeros(nx,ny,nz);t33sgs=zeros(nx,ny,nz);

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
funcmom=rho_bar.*(0.5*du-u_tilde).*(0.5*du+u_tilde);
delta_theta=(1/(rho_0*(du^2)))*sum(funcmom)*dy;

% mu physical
muphy_ref = (1/Re_theta);
T0        = 101325/(287*1.2);
Sconst    = 110.4/T0;
Skconst   = 194/T0;  %For Pr
muphy     = muphy_ref .*((T./T_ref).^1.5).*((T_ref+Sconst)./(T+Sconst));
%muphy = mu - mustar;
%muphy_avg(k,:) = avgxz(muphy,nx,nz)./(du * delta_theta * rho_0);

% LAD
S11 = dudx; S22 = dvdy ; S33=dwdz;
S12 = 0.5*(dudy+dvdx);
S13 = 0.5*(dudz+dwdx);
S23 = 0.5*(dvdz+dwdy);

% ------------- mu*(LAD) -------------------%
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
%cmu = 0.00; 
%mustar = cmu .* rho .* term; 

% Filter twice
%mustar = gaussian_filter(mustar,nx,ny,nz);
%mustar = gaussian_filter(mustar,nx,ny,nz);

mustar = mu - muphy;
%mustar_avg(k,:) = avgxz(mustar,nx,nz)./(du * delta_theta * rho_0);


% ------------------- SGS --------------------%
if(sgsid==1)
   num1 = (dudx.*(dudx.*S11 + dvdx.*S12 + dwdx.*S13) + dvdx.*(dudx.*S12 + dvdx.*S22 + dwdx.*S23) + dwdx.*(dudx.*S13 + dvdx.*S23 + dwdx.*S33)).*((dx^2) / 12);
   num2 = (dudy.*(dudy.*S11 + dvdy.*S12 + dwdy.*S13) + dvdy.*(dudy.*S12 + dvdy.*S22 + dwdy.*S23) + dwdy.*(dudy.*S13 + dvdy.*S23 + dwdy.*S33)).*((dy^2) / 12);
   num3 = (dudz.*(dudz.*S11 + dvdz.*S12 + dwdz.*S13) + dvdz.*(dudz.*S12 + dvdz.*S22 + dwdz.*S23) + dwdz.*(dudz.*S13 + dvdz.*S23 + dwdz.*S33)).*((dz^2) / 12);
   numer = -(num1+num2+num3);

   denom = (dudx.^2 + dudy.^2 + dudz.^2) + (dvdx.^2 + dvdy.^2 + dvdz.^2) + (dwdx.^2 + dwdy.^2 + dwdz.^2) + (1e-32);

   musgs = rho .* max(numer./denom,0)  ;

elseif(sgsid==2)
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

   musgs = (sigma3.*(sigma1 - sigma2).*(sigma2 - sigma3)./(sigma1sq + 1e-15)) * (deltales^2) .* rho;
   
else
   %display('I am here in Smag')
   %modS_sq =  2*(S11.^2 + S22.^2 + S33.^2 + 2*(S12.^2 + S13.^2 + S23.^2));
   %musgs = sqrt(modS_sq) * (deltales^2) .* rho;

end

%Physical Dissipation------------------------------
bulk_phy = zeros(nx,ny,nz);
lambda=bulk_phy-(2/3)*muphy; bambda=(4/3)*muphy + bulk_phy;
tau11=bambda.*dudx+lambda.*(dvdy+dwdz); tau12=muphy.*(dudy+dvdx); tau13=muphy.*(dudz+dwdx);
tau22=bambda.*dvdy+lambda.*(dudx+dwdz); tau23=muphy.*(dvdz+dwdy); tau33=bambda.*dwdz+lambda.*(dvdy+dudx);

[tau11_bar,tau11_rey_fluct,tau12_bar,tau12_rey_fluct,tau13_bar,tau13_rey_fluct,tau22_bar,tau22_rey_fluct,tau23_bar,tau23_rey_fluct,tau33_bar,tau33_rey_fluct]=tau_avg_fluct(tau11,tau12,tau13,tau22,tau23,tau33,nx,ny,nz);

func_tau=tau11_rey_fluct.*du_ppdx+tau12_rey_fluct.*du_ppdy+tau13_rey_fluct.*du_ppdz+tau12_rey_fluct.*dv_ppdx+tau22_rey_fluct.*dv_ppdy+tau23_rey_fluct.*dv_ppdz+tau13_rey_fluct.*dw_ppdx+tau23_rey_fluct.*dw_ppdy+tau33_rey_fluct.*dw_ppdz;

diss_phy=avgxz(func_tau,nx,nz);
diss_phy_int(k,1) = (sum(diss_phy)*dy)/(du^3);

%Artificial Dissipation -------------------------------
lambda=bulk-(2/3)*mustar; bambda=(4/3)*mustar + bulk;
tau11=bambda.*dudx+lambda.*(dvdy+dwdz); tau12=mustar.*(dudy+dvdx); tau13=mustar.*(dudz+dwdx);
tau22=bambda.*dvdy+lambda.*(dudx+dwdz); tau23=mustar.*(dvdz+dwdy); tau33=bambda.*dwdz+lambda.*(dvdy+dudx);

[tau11_bar,tau11_rey_fluct,tau12_bar,tau12_rey_fluct,tau13_bar,tau13_rey_fluct,tau22_bar,tau22_rey_fluct,tau23_bar,tau23_rey_fluct,tau33_bar,tau33_rey_fluct]=tau_avg_fluct(tau11,tau12,tau13,tau22,tau23,tau33,nx,ny,nz);

func_tau=tau11_rey_fluct.*du_ppdx+tau12_rey_fluct.*du_ppdy+tau13_rey_fluct.*du_ppdz+tau12_rey_fluct.*dv_ppdx+tau22_rey_fluct.*dv_ppdy+tau23_rey_fluct.*dv_ppdz+tau13_rey_fluct.*dw_ppdx+tau23_rey_fluct.*dw_ppdy+tau33_rey_fluct.*dw_ppdz;

diss_star=avgxz(func_tau,nx,nz);
diss_star_int(k,1) = (sum(diss_star)*dy)/(du^3);

%Subgrid Dissipation -------------------------------
tau11=-t11sgs; tau12=-t12sgs; tau13=-t13sgs; tau22=-t22sgs; tau23=-t23sgs; tau33=-t33sgs;
[tau11_bar,tau11_rey_fluct,tau12_bar,tau12_rey_fluct,tau13_bar,tau13_rey_fluct,tau22_bar,tau22_rey_fluct,tau23_bar,tau23_rey_fluct,tau33_bar,tau33_rey_fluct]=tau_avg_fluct(tau11,tau12,tau13,tau22,tau23,tau33,nx,ny,nz);

func_tau=tau11_rey_fluct.*du_ppdx+tau12_rey_fluct.*du_ppdy+tau13_rey_fluct.*du_ppdz+tau12_rey_fluct.*dv_ppdx+tau22_rey_fluct.*dv_ppdy+tau23_rey_fluct.*dv_ppdz+tau13_rey_fluct.*dw_ppdx+tau23_rey_fluct.*dw_ppdy+tau33_rey_fluct.*dw_ppdz;

diss_sgs=avgxz(func_tau,nx,nz);
diss_sgs_int(k,1) = (sum(diss_sgs)*dy)/(du^3);

k=k+1;
end


fname = strcat('/home/vsj/Codes/diss_time/diss_time_',caseid,'_',sprintf('%1.1f',Mc),'.dat');
fileID = fopen(fname,'a');
for p=1:size(time_total,2)
fprintf(fileID,'%04d %2.8f %2.8f %2.8f\r\n',time_total(p),diss_phy_int(p),diss_star_int(p),diss_sgs_int(p));
end
fclose(fileID);

figure(1),clf
plot(time_total,diss_phy_int,'or','LineWidth',2,'DisplayName','\epsilon_{phy}');hold on
plot(time_total,diss_star_int,'>b','LineWidth',2,'DisplayName','\epsilon_{*}');hold on
plot(time_total,diss_sgs_int,'^g','LineWidth',2,'DisplayName','\epsilon_{sgs}');hold on
xlab=xlabel('\boldmath$Time$','Interpreter','latex');set(xlab,'Fontsize',20);
%h=ylabel('\boldmath$\sqrt{|R_{22}|}/{\Delta{u}}$');set(h,'Interpreter','latex');set(h,'Fontsize',20);

legend show;legend('Location','northeast','box','off','FontSize',10)

fnamepng = strcat('/home/vsj/Codes/diss_time/diss_time_',caseid,'_',sprintf('%1.1f',Mc),'.png');
screen2jpeg(fnamepng)

